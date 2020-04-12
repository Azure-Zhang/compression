// ------------------------------------------------------------------
//   sections.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "sections.h"
#include "buffer.h"
#include "file.h"
#include "vblock.h"
#include "endianness.h"
#include "random_access_vcf.h"

// ZIP only: create section list that goes into the genozip header, as we are creating the sections
void sections_add_to_list (VBlock *vb, const SectionHeader *header)
{
    DictIdType dict_id = { 0 };
    bool is_dict = (section_type_is_dictionary (header->section_type));

    if (is_dict)                                          dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    else if (section_type_is_b250 (header->section_type)) dict_id = ((SectionHeaderBase250    *)header)->dict_id;

    // 1. if this is a vcf_header, random_access or genozip_header - it goes directly into the z_file by the I/O thread
    //    before or after all the compute threads are operational
    // 2. if this is a dictionary - it goes directly into z_file by the Compute thread while merge holds the mutex:
    //    mtf_merge_in_vb_ctx_one_dict_id -> zfile_compress_dictionary_data
    // 3. if we this section is part of a VB (other than a dictionary), we store the entry within the VB and merge it to
    //    the z_file in the correct order of VBs by the I/O thread after the compute thread is finished.
    //
    // offsets in case 2 and 3 are relative to their buffer at this point, and will be updated later

    uint64_t offset;
    Buffer *buf;
    char *name;
    VBlockP alc_vb;
    if (!vb->vblock_i) {  // case 1 - vcf header, random access, genotype header
        buf    = &z_file->section_list_buf;
        alc_vb = evb; // z_file buffer goes to evb
        name   = "z_file->section_list_buf";
        offset = z_file->disk_so_far + vb->z_data.len;
    }
    else if (is_dict) {          // case 2 - dictionaries
        buf    = &z_file->section_list_dict_buf;
        alc_vb = evb; // z_file buffer goes to evb
        name   = "z_file->section_list_dict_buf";
        offset = z_file->dict_data.len;
    }
    else {                       // case 3 - VB content
        buf    = &vb->section_list_buf;
        alc_vb = vb;
        name   = "section_list_buf";
        offset = vb->z_data.len;
    }

    buf_alloc (alc_vb, buf, MAX (buf->len + 1, 50) * sizeof(SectionListEntry), 2, name, vb->vblock_i);
    
    SectionListEntry *ent = &NEXTENT (SectionListEntry, *buf);
    ent->section_type     = header->section_type;
    ent->vblock_i         = BGEN32 (header->vblock_i); // big endian in header - convert back to native
    ent->dict_id          = dict_id;
    ent->offset           = offset;  // this is a partial offset (within d) - we will correct it later
}

// Called by ZIP I/O thread. concatenates a vb or dictionary section list to the z_file sectinon list - just before 
// writing those sections to the disk. we use the current disk position to update the offset
void sections_list_concat (VBlock *vb, BufferP section_list_buf)
{
    buf_alloc (evb, &z_file->section_list_buf, 
              (z_file->section_list_buf.len + section_list_buf->len) * sizeof(SectionListEntry), 2, 
              "z_file->section_list_buf", 0);
  
    SectionListEntry *dst = AFTERENT (SectionListEntry, z_file->section_list_buf);
    SectionListEntry *src = FIRSTENT (SectionListEntry, *section_list_buf);

    // update the offset
    for (unsigned i=0; i < section_list_buf->len; i++)
        src[i].offset += z_file->disk_so_far;

    // copy all entries
    memcpy (dst, src, section_list_buf->len * sizeof(SectionListEntry));
    z_file->section_list_buf.len += section_list_buf->len;

    buf_free (section_list_buf);
}

// called by PIZ I/O thread: zfile_read_on_vb
uint8_t sections_count_info_b250s (unsigned vb_i)
{
    SectionListEntry *sl = ARRAY (SectionListEntry, z_file->section_list_buf);    

    // skip to the first SEC_VCF_INFO_SF_B250 (if there is one...)
    while (z_file->sl_cursor < z_file->section_list_buf.len &&
           sl[z_file->sl_cursor].vblock_i == vb_i &&
           sl[z_file->sl_cursor].section_type != SEC_VCF_INFO_SF_B250) 
        z_file->sl_cursor++;

    // count the SEC_VCF_INFO_SF_B250 sections
    uint32_t start = z_file->sl_cursor;
    while (sl[z_file->sl_cursor].section_type == SEC_VCF_INFO_SF_B250) z_file->sl_cursor++;

    return (uint8_t)(z_file->sl_cursor - start);
}

// called by PIZ I/O to know if next up is a VB Header or VCF Header or EOF
SectionType sections_get_next_header_type (SectionListEntry **sl_ent, 
                                           bool *skipped_vb,   // out (VB only) - true if this vb should be skipped
                                           Buffer *region_ra_intersection_matrix) // out (VB only) - a bytemap - rows are ra's of this VB, columns are regions, a cell is 1 if there's an intersection
{
    // find the next VB or VCF header section
    if (skipped_vb) *skipped_vb = false;

    while (z_file->sl_cursor < z_file->section_list_buf.len) {
        *sl_ent = ENT (SectionListEntry, z_file->section_list_buf, z_file->sl_cursor++);
 
        SectionType sec_type = (*sl_ent)->section_type;
        if (sec_type == SEC_TXT_HEADER) 
            return SEC_TXT_HEADER;

        if (sec_type == SEC_VCF_VB_HEADER) {
            if (random_access_is_vb_included ((*sl_ent)->vblock_i, region_ra_intersection_matrix))
                return SEC_VCF_VB_HEADER;
            
            else if (skipped_vb) *skipped_vb = true;
        }
    }

    return SEC_EOF; // no more headers
}

// dictionary section iterator. returns true if another dictionary was found.
bool sections_get_next_dictionary (SectionListEntry **sl_ent) // if *sl_ent==NULL - initialize cursor
{
    if (! *sl_ent) z_file->sl_dir_cursor = 0;
    
    // find the next VB or VCF header section
    while (z_file->sl_dir_cursor < z_file->section_list_buf.len) {
        *sl_ent = ENT (SectionListEntry, z_file->section_list_buf, z_file->sl_dir_cursor++);
        if (section_type_is_dictionary((*sl_ent)->section_type)) 
            return true;
    }

    return false; // no more dictionaries
}

// called by PIZ I/O : zfile_vcf_read_one_vb. Sets *sl_ent to the first section of this vb_i, and returns its offset
uint64_t sections_vb_first (uint32_t vb_i, SectionListEntry **sl_ent)
{
    unsigned i=0; for (; i < z_file->section_list_buf.len; i++) {
        *sl_ent = ENT (SectionListEntry, z_file->section_list_buf, i);
        if ((*sl_ent)->vblock_i == vb_i) break; // found!
    }

    ASSERT (i < z_file->section_list_buf.len, "Error in sections_get_next_vb_section: cannot find any section for vb_i=%u", vb_i);

    return (*sl_ent)->offset;
}

// called by PIZ I/O : zfile_vcf_read_one_vb. Returns the offset of the next section with this vb_i
uint64_t sections_vb_next (SectionListEntry **sl_ent /* in / out */)
{
    ASSERT0 (*sl_ent, "Error in sections_vb_next: *sl_ent is NULL");

    uint32_t vb_i = (*sl_ent)->vblock_i;

    // get next section - just verify that it belongs to this vb_i
    ASSERT (1 + *sl_ent <= LASTENT (SectionListEntry, z_file->section_list_buf), 
            "Error in sections_get_next_vb_section: there are no more sections (called for vb_i=%u)", vb_i);

    (*sl_ent)++;

    ASSERT ((*sl_ent)->vblock_i == vb_i, "Error in sections_get_next_vb_section: vb_i=%u has no more sections", vb_i);

    return (*sl_ent)->offset;
}

// called by PIZ I/O when splitting a concatenated file - to know if there are any more VCF components remaining
bool sections_has_more_vcf_components()
{
    return z_file->sl_cursor==0 || 
           ENT (SectionListEntry, z_file->section_list_buf, z_file->sl_cursor-1)->section_type == SEC_TXT_HEADER;
}

void BGEN_sections_list()
{
    SectionListEntry *ent = ARRAY (SectionListEntry, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) {
        ent[i].vblock_i = BGEN32 (ent[i].vblock_i);
        ent[i].offset          = BGEN64 (ent[i].offset);
    }
}

void sections_show_gheader (SectionHeaderGenozipHeader *header)
{
    unsigned num_sections = BGEN32 (header->num_sections);
    char size_str[50];

    fprintf (stderr, "Below are the contents of the genozip header. This is the output of --show-gheader:\n");
    fprintf (stderr, "  genozip_version: %u\n",               header->genozip_version);
    fprintf (stderr, "  data_type: %s\n",                     dt_name (BGEN16 (header->data_type)));
    fprintf (stderr, "  encryption_type: %s\n",               encryption_name (header->encryption_type)); 
    fprintf (stderr, "  num_samples: %u\n",                   BGEN32 (header->num_samples));
    fprintf (stderr, "  uncompressed_data_size: %s\n",        buf_display_uint (BGEN64 (header->uncompressed_data_size), size_str));
    fprintf (stderr, "  num_items_concat: %"PRIu64"\n",       BGEN64 (header->num_items_concat));
    fprintf (stderr, "  num_sections: %u\n",                  num_sections);
    fprintf (stderr, "  num_vcf_components: %u\n",            BGEN32 (header->num_vcf_components));
    fprintf (stderr, "  md5_hash_concat: %s\n",               md5_display (&header->md5_hash_concat, false));
    fprintf (stderr, "  created: %*s\n",                      -FILE_METADATA_LEN, header->created);

    fprintf (stderr, "  sections:\n");

    SectionListEntry *ents = ARRAY (SectionListEntry, z_file->section_list_buf);

    for (unsigned i=0; i < num_sections; i++) {
     
        uint64_t this_offset = ents[i].offset;
        uint64_t next_offset = (i < num_sections-1) ? ents[i+1].offset : z_file->disk_so_far;

        fprintf (stderr, "    %3u. %-24.24s %*.*s vb_i=%u offset=%"PRIu64" size=%"PRId64"\n", 
                 i, st_name(ents[i].section_type), 
                 -DICT_ID_LEN, DICT_ID_LEN, ents[i].dict_id.num ? dict_id_printable (ents[i].dict_id).id : ents[i].dict_id.id, 
                 ents[i].vblock_i, this_offset, next_offset - this_offset);
    }
}

static const char *sections_type_name (unsigned item, const char **names, unsigned num_names)
{
    if (item > num_names) {
        static char str[50];
        sprintf (str, "%u (out of range)", item);
        return str;
    }
    
    return names[item];    
}

const char *st_name(SectionType sec_type)
{
    static const char *names[NUM_SEC_TYPES] = SECTIONTYPE_NAMES;
    
    if (sec_type == SEC_EOF) return "SEC_EOF";
    
    return sections_type_name (sec_type, names, sizeof(names)/sizeof(names[0]));
}

const char *dt_name (unsigned data_type)
{
    static const char *names[NUM_DATATYPES] = DATATYPE_NAMES;
    return sections_type_name (data_type, names, sizeof(names)/sizeof(names[0]));
}

const char *encryption_name (unsigned encryption_type)
{
    static const char *names[NUM_ENCRYPTION_TYPES] = ENCRYPTION_TYPE_NAMES;
    return sections_type_name (encryption_type, names, sizeof(names)/sizeof(names[0]));
}

// called by PIZ I/O
uint64_t sections_get_offset_first_section_of_type (SectionType st)
{
    SectionListEntry *sl = ARRAY (SectionListEntry, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++)
        if (sl[i].section_type == st) return sl[i].offset;

    ABORT ("Error in sections_get_offset_first_section_of_type: Cannot find section_type=%s in z_file", st_name (st));
    return 0; // never reaches here - squash compiler warning
}
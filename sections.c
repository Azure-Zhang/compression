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
#include "random_access.h"
#include "strings.h"
#include "crypt.h"
#include "dict_id.h"
#include "zfile.h"

// ZIP only: create section list that goes into the genozip header, as we are creating the sections
void sections_add_to_list (VBlock *vb, const SectionHeader *header)
{
    DictId dict_id = DICT_ID_NONE;

    if      (header->section_type == SEC_DICT ) dict_id = ((SectionHeaderDictionary *)header)->dict_id;
    else if (header->section_type == SEC_B250 ) dict_id = ((SectionHeaderCtx        *)header)->dict_id;
    else if (header->section_type == SEC_LOCAL) dict_id = ((SectionHeaderCtx        *)header)->dict_id;

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
    if (header->section_type == SEC_DICT) {          // case 1 - dictionaries
        buf    = &z_file->section_list_dict_buf;
        alc_vb = evb; // z_file buffer goes to evb
        name   = "z_file->section_list_dict_buf";
        offset = z_file->dict_data.len;
    }
    else if (!vb->vblock_i) {  // case 2 - vcf header, random access, genotype header
        buf    = &z_file->section_list_buf;
        alc_vb = evb; // z_file buffer goes to evb
        name   = "z_file->section_list_buf";
        offset = z_file->disk_so_far + vb->z_data.len;
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
void sections_list_bind (VBlock *vb, BufferP section_list_buf)
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

// called by PIZ I/O to know if next up is a VB Header or VCF Header or EOF
SectionType sections_get_next_header_type (SectionListEntry **sl_ent, 
                                           bool *skipped_vb,   // out (VB only) - true if this vb should be skipped
                                           Buffer *region_ra_intersection_matrix) // out (VB only) - a bytemap - rows are ra's of this VB, columns are regions, a cell is 1 if there's an intersection
{
    // find the next VB or TXT header section
    if (skipped_vb) *skipped_vb = false;

    while (z_file->sl_cursor < z_file->section_list_buf.len) {
        *sl_ent = ENT (SectionListEntry, z_file->section_list_buf, z_file->sl_cursor++);
 
        SectionType sec_type = (*sl_ent)->section_type;
        if (sec_type == SEC_TXT_HEADER) 
            return sec_type;

        if (sec_type == SEC_VB_HEADER) {
            if (random_access_is_vb_included ((*sl_ent)->vblock_i, region_ra_intersection_matrix))
                return SEC_VB_HEADER;
            
            else if (skipped_vb) *skipped_vb = true;
        }
    }

    return SEC_NONE; // no more headers
}

// section iterator. returns true if another section of this type was found.
bool sections_get_next_section_of_type (SectionListEntry **sl_ent, uint32_t *cursor, SectionType st1, SectionType st2) // if *sl_ent==NULL - initialize cursor
{
    // case: first time
    if (! *sl_ent) {
        *cursor = 0;
        while (*cursor < (uint32_t)z_file->section_list_buf.len) {
            *sl_ent = ENT (SectionListEntry, z_file->section_list_buf, (*cursor)++);
            if ((*sl_ent)->section_type == st1 || (*sl_ent)->section_type == st2) 
                return true;
        }
    }
        
    if (*cursor == z_file->section_list_buf.len) return false; 

    *sl_ent = ENT (SectionListEntry, z_file->section_list_buf, (*cursor)++);
    return (*sl_ent)->section_type == st1 || (*sl_ent)->section_type == st2;
}

// returns the next section's type without moving the cursor
SectionType sections_peek (uint32_t cursor)
{
    if (cursor == z_file->section_list_buf.len) return SEC_NONE;

    return ENT (SectionListEntry, z_file->section_list_buf, cursor)->section_type;
}

bool sections_seek_to (SectionType st)
{
    SectionListEntry *sl=NULL;
    unsigned i=0; for (; i < z_file->section_list_buf.len; i++) {
        sl = ENT (SectionListEntry, z_file->section_list_buf, i);
        if (sl->section_type == st) {
            file_seek (z_file, sl->offset, SEEK_SET, false);
            return true; // section found
        }
    }

    return false; // section not found
}

// count how many sections we have of a certain type
uint32_t sections_count_sections (SectionType st)
{
    uint32_t count=0;
    for (uint32_t i=0; i < z_file->section_list_buf.len; i++) 
        if (ENT (SectionListEntry, z_file->section_list_buf, i)->section_type == st)
            count++;

    return count;
}

// called by PIZ I/O : vcf_zfile_read_one_vb. Sets *sl_ent to the first section of this vb_i, and returns its offset
SectionListEntry *sections_vb_first (uint32_t vb_i)
{
    SectionListEntry *sl=NULL;
    unsigned i=0; for (; i < z_file->section_list_buf.len; i++) {
        sl = ENT (SectionListEntry, z_file->section_list_buf, i);
        if (sl->vblock_i == vb_i) break; // found!
    }

    ASSERT (i < z_file->section_list_buf.len, "Error in sections_get_next_vb_section: cannot find any section for vb_i=%u", vb_i);

    return sl;
}

bool sections_has_reference(void)
{
    for (int i=z_file->section_list_buf.len-1; i >= 0; i--) { // search backwards as the reference sections are near the end
        SectionType st = ENT (SectionListEntry, z_file->section_list_buf, i)->section_type;
        if (st == SEC_REFERENCE)
            return true; // found
        else if (st == SEC_DICT)
            return false; // we arrived at a SEC_DICT without seeing any reference, so there isn't any
    }
    
    return false;
}

// get genome size - this is determined by gpos+num_bases of the last SEC_REFERENCE or SEC_REF_IS_SET section
// NOT THREAD SAFE - can only be called by I/O thread
int64_t sections_get_genome_size (void)
{
    for (int i=z_file->section_list_buf.len-1; i >= 0; i--) { // search backwards as the reference sections are near the end
        const SectionListEntry *sl = ENT (const SectionListEntry, z_file->section_list_buf, i);
        
        if (sl->section_type == SEC_REFERENCE) {
            // check if we have a SEC_REF_IS_SET section - in that case, take num_bases from SEC_REF_IS_SET
            // as the num_bases in SEC_REFERENCE indicates a compacted section
            if (i > 0 && (sl-1)->section_type == SEC_REF_IS_SET) sl--;

            SectionHeaderReference *header = zfile_read_section_header (sl->offset, sizeof (SectionHeaderReference));

            return (int64_t)(BGEN64(header->gpos) + BGEN32 (header->num_bases));
        }
    }

    ABORT0 ("Error in sections_get_genome_size: cannot find a SEC_REFERENCE section");
    return 0;
}

// called by refhash_initialize - get details of the refhash ahead of loading it from the reference file 
// NOT THREAD SAFE - can only be called by I/O thread
void sections_get_refhash_details (uint32_t *num_layers, uint32_t *base_layer_bits) // optional outs
{
    ASSERT0 (flag_reading_reference, "Error in sections_get_refhash_details: can only be called while reading reference");

    for (int i=z_file->section_list_buf.len-1; i >= 0; i--) { // search backwards as the refhash sections are near the end
        SectionListEntry *sl = ENT (SectionListEntry, z_file->section_list_buf, i);
        if (sl->section_type == SEC_REF_HASH) {

            SectionHeaderRefHash *header = zfile_read_section_header (sl->offset, sizeof (SectionHeaderRefHash));
            if (num_layers) *num_layers = header->num_layers;
            if (base_layer_bits) *base_layer_bits = header->layer_bits + header->layer_i; // layer_i=0 is the base layer, layer_i=1 has 1 bit less etc
            return;
        }
        else if (sl->section_type == SEC_REFERENCE)
            break; // we arrived at a SEC_REFERENCE - there won't be any more SEC_REF_HASH sections
    }

    ABORT ("Error in sections_get_refhash_details: can't find SEC_REF_HASH sections in %s", z_name);
}


// called by PIZ I/O when splitting a bound file - to know if there are any more VCF components remaining
bool sections_has_more_components()
{
    return z_file->sl_cursor==0 || 
           ENT (SectionListEntry, z_file->section_list_buf, z_file->sl_cursor-1)->section_type == SEC_TXT_HEADER;
}

void BGEN_sections_list()
{
    ARRAY (SectionListEntry, ent, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++) {
        ent[i].vblock_i = BGEN32 (ent[i].vblock_i);
        ent[i].offset   = BGEN64 (ent[i].offset);
    }
}

void sections_show_gheader (SectionHeaderGenozipHeader *header)
{
    if (flag_reading_reference) return; // don't show gheaders of reference file
    
    unsigned num_sections = BGEN32 (header->num_sections);
    char size_str[50];

    fprintf (stderr, "Contents of the genozip header (output of --show-gheader) of %s:\n", z_name);
    fprintf (stderr, "  genozip_version: %u\n",         header->genozip_version);
    fprintf (stderr, "  data_type: %s\n",               dt_name (BGEN16 (header->data_type)));
    fprintf (stderr, "  encryption_type: %s\n",         encryption_name (header->encryption_type)); 
    fprintf (stderr, "  num_samples: %u\n",             BGEN32 (header->num_samples));
    fprintf (stderr, "  uncompressed_data_size: %s\n",  str_uint_commas (BGEN64 (header->uncompressed_data_size), size_str));
    fprintf (stderr, "  num_items_bound: %"PRIu64"\n", BGEN64 (header->num_items_bound));
    fprintf (stderr, "  num_sections: %u\n",            num_sections);
    fprintf (stderr, "  num_components: %u\n",          BGEN32 (header->num_components));
    fprintf (stderr, "  md5_hash_bound: %s\n",          md5_display (header->md5_hash_bound));
    fprintf (stderr, "  created: %*s\n",                -FILE_METADATA_LEN, header->created);
    fprintf (stderr, "  license_hash: %s\n",            md5_display (header->license_hash));
    fprintf (stderr, "  reference filename: %*s\n",     -REF_FILENAME_LEN, header->ref_filename);
    fprintf (stderr, "  reference file hash: %s\n",     md5_display (header->ref_file_md5));

    fprintf (stderr, "  sections:\n");

    ARRAY (SectionListEntry, ents, z_file->section_list_buf);

    for (unsigned i=0; i < num_sections; i++) {
     
        uint64_t this_offset = ents[i].offset;
        uint64_t next_offset;
        
        if (i < num_sections-1) 
            next_offset = ents[i+1].offset;

        else // we're at the last section genozip header+footer
            next_offset = this_offset + BGEN32 (header->h.data_compressed_len) + BGEN32 (header->h.compressed_offset) + sizeof (SectionFooterGenozipHeader);

        fprintf (stderr, "    %3u. %-24.24s %*.*s vb_i=%u offset=%"PRIu64" size=%"PRId64"\n", 
                 i, st_name(ents[i].section_type), 
                 -DICT_ID_LEN, DICT_ID_LEN, ents[i].dict_id.num ? dict_id_printable (ents[i].dict_id).id : ents[i].dict_id.id, 
                 ents[i].vblock_i, this_offset, next_offset - this_offset);
    }
}

const char *st_name(SectionType sec_type)
{
    static const struct {const char *name; } abouts[NUM_SEC_TYPES] = SECTIONTYPE_ABOUT;
    
    if (sec_type == SEC_NONE) return "SEC_NONE";
    
    return type_name (sec_type, &abouts[sec_type].name , sizeof(abouts)/sizeof(abouts[0]));
}

// called by PIZ I/O
SectionListEntry *sections_get_offset_first_section_of_type (SectionType st)
{
    ARRAY (SectionListEntry, sl, z_file->section_list_buf);

    for (unsigned i=0; i < z_file->section_list_buf.len; i++)
        if (sl[i].section_type == st) return &sl[i];

    ABORT ("Error in sections_get_offset_first_section_of_type: Cannot find section_type=%s in z_file", st_name (st));
    return 0; // never reaches here - squash compiler warning
}
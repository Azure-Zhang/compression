// ------------------------------------------------------------------
//   vblock.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// vb stands for VBlock - it started its life as VBlockVCF when genozip could only compress VCFs, but now
// it means a block of lines from the text file. 

#include "genozip.h"
#include "move_to_front.h"
#include "vblock.h"
#include "file.h"

// pool of VBs allocated based on number of threads
static VBlockPool *pool = NULL;

// one VBlock outside of pool
VBlock *evb = NULL;

// cleanup vb and get it ready for another usage (without freeing memory held in the Buffers)
static void vb_vcf_release_vb (VBlockVCF *vb); 
static void vb_sam_release_vb (VBlockSAM *vb);
void vb_release_vb (VBlock **vb_p) 
{
    if (! *vb_p) return; // nothing to release

    VBlock *vb = *vb_p;
    *vb_p = NULL;

    vb->num_lines = vb->first_line = vb->vblock_i = vb->txt_data_next_offset = 0;
    vb->vb_data_size = vb->vb_data_read_size = 0;
    vb->ready_to_dispatch = vb->is_processed = false;
    vb->z_next_header_i = 0;
    vb->num_dict_ids = 0;
    
    memset(vb->txt_section_bytes, 0, sizeof(vb->txt_section_bytes));
    memset(vb->z_section_bytes, 0, sizeof(vb->z_section_bytes));
    memset(vb->z_num_sections, 0, sizeof(vb->z_num_sections));
    memset(vb->z_section_entries, 0, sizeof(vb->z_section_entries));

    memset (&vb->profile, 0, sizeof (vb->profile));

    buf_free(&vb->compressed);
    buf_free(&vb->txt_data);
    buf_free(&vb->txt_data_spillover);
    buf_free(&vb->z_data);
    buf_free(&vb->z_section_headers);
    buf_free(&vb->spiced_pw);
    buf_free(&vb->show_headers_buf);
    buf_free(&vb->show_b250_buf);
    buf_free(&vb->section_list_buf);

    for (unsigned i=0; i < MAX_DICTS; i++) 
        if (vb->mtf_ctx[i].dict_id.num)
            mtf_free_context (&vb->mtf_ctx[i]);

    for (unsigned i=0; i < NUM_COMPRESS_BUFS; i++)
        buf_free (&vb->compress_bufs[i]);
        
    vb->in_use = false; // released the VB back into the pool - it may now be reused

    if      (vb->data_type == DATA_TYPE_VCF) vb_vcf_release_vb ((VBlockVCF*)vb);
    else if (vb->data_type == DATA_TYPE_SAM) vb_sam_release_vb ((VBlockSAM*)vb);

    // STUFF THAT PERSISTS BETWEEN VBs (i.e. we don't free / reset):
    // vb->num_lines_alloced
    // vb->buffer_list : we DON'T free this because the buffers listed are still available and going to be re-used/
    //                   we have logic in vb_get_vb() to update its vb_i
    // vb->num_sample_blocks : we keep this value as it is needed by vb_cleanup_memory, and it doesn't change
    //                         between VBs of a file or concatenated files.
    // vb->data_type : type of this vb 
}

void vb_create_pool (unsigned num_vbs)
{
    ASSERT (!pool || num_vbs==pool->num_vbs, 
            "Error: vb pool already exists, but with the wrong number of vbs - expected %u but it has %u", num_vbs, pool->num_vbs);

    if (!pool)  {
        // allocation includes array of pointers (initialized to NULL)
        pool = (VBlockPool *)calloc (1, sizeof (VBlockPool) + num_vbs * sizeof (VBlockVCF *)); // note we can't use Buffer yet, because we don't have VBs yet...
        ASSERT0 (pool, "Error: failed to calloc pool");

        pool->num_vbs = num_vbs; 
    }
}

VBlockPool *vb_get_pool (void)
{
    return pool;
}

void vb_external_vb_initialize(void)
{
    ASSERT0 (!evb, "Error: evb already initialized");

    evb = calloc (1, sizeof (VBlock));
    ASSERT0 (evb, "Error: failed to calloc evb");
    evb->data_type = DATA_TYPE_NONE;
    evb->id = -1;
}

// allocate an unused vb from the pool. seperate pools for zip and unzip
VBlock *vb_get_vb (unsigned vblock_i)
{
    // see if there's a VB avaiable for recycling
    unsigned vb_i; for (vb_i=0; vb_i < pool->num_vbs; vb_i++) {
    
        // free if this is a VB allocated by a previous file, with a different data type
        if (pool->vb[vb_i] && pool->vb[vb_i]->data_type != z_file->data_type) {
            FREE (pool->vb[vb_i]);
            pool->vb[vb_i] = NULL; 
            pool->num_allocated_vbs--; // we will immediately allocate and increase this back
        }
        
        if (!pool->vb[vb_i]) { // VB is not allocated - allocate it
            switch (z_file->data_type) {
                case DATA_TYPE_VCF : pool->vb[vb_i] = calloc (sizeof (VBlockVCF), 1); break;
                case DATA_TYPE_SAM : pool->vb[vb_i] = calloc (sizeof (VBlockSAM), 1); break;
                default            : ABORT ("Error in vb_get_vb: Invalid data_type=%d", z_file->data_type);
            }
            pool->num_allocated_vbs++;
            pool->vb[vb_i]->data_type = z_file->data_type;
        }

        if (!pool->vb[vb_i]->in_use) {
            pool->vb[vb_i]->id = vb_i;
            break;
        }
    }

    ASSERT (vb_i < pool->num_vbs, "Error: VB pool is full - it already has %u VBs", pool->num_vbs)

    pool->vb[vb_i]->in_use           = true;
    pool->vb[vb_i]->vblock_i         = vblock_i;
    pool->vb[vb_i]->buffer_list.vb_i = vblock_i;

    return pool->vb[vb_i];
}

// free memory allocations that assume subsequent files will have the same number of samples.
static void vb_vcf_cleanup_memory (VBlockVCF *vb);
static void vb_sam_cleanup_memory (VBlockSAM *vb);
void vb_cleanup_memory (void)
{

    for (unsigned vb_i=0; vb_i < pool->num_vbs; vb_i++) {
        VBlock *vb = pool->vb[vb_i];
        if      (vb && vb->data_type == DATA_TYPE_VCF) vb_vcf_cleanup_memory ((VBlockVCF *)vb);
        else if (vb && vb->data_type == DATA_TYPE_SAM) vb_sam_cleanup_memory ((VBlockSAM *)vb);
    }

    global_vcf_num_samples = 0;
}

unsigned vb_vcf_num_samples_in_sb (const VBlockVCF *vb, unsigned sb_i)
{
    // case: last block has less than a full block of samples
    if (sb_i == vb->num_sample_blocks-1 && global_vcf_num_samples % vb->num_samples_per_block)
        return global_vcf_num_samples % vb->num_samples_per_block;

    else
        return vb->num_samples_per_block;
} 

unsigned vb_vcf_num_sections(VBlockVCF *vb) 
{
    return 1 + vb->has_genotype_data + (vb->phase_type == PHASE_MIXED_PHASED) + (vb->num_haplotypes_per_line > 0);
}


//--------------------------------
// SAM stuff
//--------------------------------

// cleanup vb (except common) and get it ready for another usage (without freeing memory held in the Buffers)
static void vb_vcf_release_vb (VBlockVCF *vb) 
{
    // note: vb->data_line is not freed but rather used by subsequent vbs
    if (command == ZIP && vb->data_lines)
        memset (vb->data_lines, 0, sizeof(ZipDataLineVCF) * vb->num_lines_alloced);

    else if (command != ZIP && vb->data_lines) {
        for (unsigned i=0; i < vb->num_lines_alloced; i++) {
            PizDataLineVCF *dl = &((PizDataLineVCF *)vb->data_lines)[i];
            
            dl->has_haplotype_data = dl->has_genotype_data = 0;
            dl->format_mtf_i = 0;

            buf_free(&dl->line);
            buf_free(&dl->v1_variant_data);
        }
    }

    for (unsigned i=0; i < vb->num_sample_blocks; i++) {
        if (vb->haplotype_sections_data) buf_free(&vb->haplotype_sections_data[i]);
        if (vb->genotype_sections_data)  buf_free(&vb->genotype_sections_data[i]);
        if (vb->phase_sections_data)     buf_free(&vb->phase_sections_data[i]);
    }

    vb->ploidy = vb->num_haplotypes_per_line = 0;
    vb->has_genotype_data = vb->has_haplotype_data = false;
    vb->phase_type = PHASE_UNKNOWN;
    vb->max_gt_line_len = vb->last_pos = 0;
    vb->curr_ra_ent = NULL; 
    vb->curr_ra_ent_is_initialized = false;

    buf_free(&vb->line_variant_data);
    buf_free(&vb->line_gt_data);
    buf_free(&vb->line_ht_data);
    buf_free(&vb->line_phase_data);

    buf_free(&vb->sample_iterator);
    buf_free(&vb->genotype_one_section_data);
    buf_free(&vb->is_sb_included);
    buf_free(&vb->genotype_section_lens_buf);

    buf_free (&vb->format_mapper_buf);
    buf_free (&vb->iname_mapper_buf);

    vb->num_info_subfields = vb->num_format_subfields = 0;

    buf_free(&vb->optimized_gl_dict);
    buf_free(&vb->haplotype_permutation_index);
    buf_free(&vb->haplotype_permutation_index_squeezed);
    
    buf_free(&vb->gt_sb_line_starts_buf);
    buf_free(&vb->gt_sb_line_lengths_buf);
    buf_free(&vb->helper_index_buf);
    buf_free(&vb->ht_columns_data);
    buf_free(&vb->format_info_buf);
    buf_free(&vb->ra_buf);
    buf_free(&vb->region_ra_intersection_matrix);
    buf_free(&vb->column_of_zeros);

    buf_free(&vb->gtshark_db_db_data);
    buf_free(&vb->gtshark_db_gt_data);
    buf_free(&vb->gtshark_exceptions_line_i);
    buf_free(&vb->gtshark_exceptions_ht_i);
    buf_free(&vb->gtshark_exceptions_allele);
    buf_free(&vb->gtshark_vcf_data);

    // backward compatibility with genozip v1
    buf_free(&vb->v1_subfields_start_buf);        
    buf_free(&vb->v1_subfields_len_buf);
    buf_free(&vb->v1_num_subfields_buf);
    buf_free(&vb->v1_variant_data_section_data);
}

// freeing arrays of buffers allocated by calloc 
static void vb_vcf_free_buffer_array (Buffer **buf_array, unsigned buf_array_len)
{
    if (! (*buf_array)) return; // array not allocated - nothing to do

    for (unsigned i=0; i < buf_array_len; i++) 
        buf_destroy (&(*buf_array)[i]);

    FREE (*buf_array);

    *buf_array = NULL;
}

// free memory allocations that assume subsequent files will have the same number of samples.
static void vb_vcf_cleanup_memory (VBlockVCF *vb)
{
    vb_vcf_free_buffer_array (&vb->genotype_sections_data, vb->num_sample_blocks);
    vb_vcf_free_buffer_array (&vb->haplotype_sections_data, vb->num_sample_blocks);
    vb_vcf_free_buffer_array (&vb->phase_sections_data, vb->num_sample_blocks);
    vb->num_sample_blocks = 0;
}

//--------------------------------
// SAM stuff
//--------------------------------

static void vb_sam_release_vb (VBlockSAM *vb)
{
    vb->num_optional_subfield_b250s = vb->next_seq = vb->next_qual = vb->next_random_pos = 0;
    vb->nm_did_i = vb->strand_did_i = vb->mc_did_i = 0;
    vb->last_pos = vb->last_rname_node_index = 0;

    memset (&vb->qname_mapper, 0, sizeof (vb->qname_mapper));
    
    // note: vb->data_line is not freed but rather used by subsequent vbs
    if (command == ZIP && vb->data_lines)
        memset (vb->data_lines, 0, sizeof(ZipDataLineSAM) * vb->num_lines_alloced);

    buf_free (&vb->random_pos_data);
    buf_free (&vb->optional_mapper_buf);
    buf_free (&vb->seq_data);
    buf_free (&vb->qual_data);
}

static void vb_sam_cleanup_memory (VBlockSAM *vb)
{
  // nothing to do
}

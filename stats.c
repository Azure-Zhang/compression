// ------------------------------------------------------------------
//   stats.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "header.h"
#include "dict_id.h"
#include "strings.h"
#include "stats.h"
#include "sections.h"
#include "file.h"

#define stats_verify_all_covered(covered) stats_verify_all_covered_do (covered, __FUNCTION__)
static void stats_verify_all_covered_do (const unsigned *covered, const char *caller)
{
    for (unsigned sec_i=0; sec_i < NUM_SEC_TYPES; sec_i++) {

        ASSERTW (covered[sec_i]==1 || !txt_file->section_bytes[sec_i], 
                 "Warning: section %s has bytes in txt_file, but covered=%u in %s", st_name (sec_i), covered[sec_i], caller);

        ASSERTW (covered[sec_i]==1 || !z_file->section_bytes[sec_i], 
                 "Warning: section %s has bytes in z_file, but covered=%u in %s", st_name (sec_i), covered[sec_i], caller);
    }
}

static void stats_show_file_metadata (void)
{
    fprintf (stderr, "\n\n");
    if (txt_file->name) fprintf (stderr, "%s file name: %s\n", dt_name (z_file->data_type), txt_file->name);
    
    char ls[30];
    if (z_file->data_type == DT_VCF) 
        fprintf (stderr, "Individuals: %u   %s: %s   Dictionaries: %u   Vblocks: %u\n", 
                 global_vcf_num_samples, DTPZ (show_sections_line_name),  
                 str_uint_commas (z_file->num_lines, ls), z_file->num_dict_ids, z_file->num_vbs);
    else
        fprintf (stderr, "%s: %s   Dictionaries: %u   Vblocks: %u\n", 
                 DTPZ (show_sections_line_name), str_uint_commas (z_file->num_lines, ls), z_file->num_dict_ids, z_file->num_vbs);
}

void stats_show_sections (void)
{
    // the order in which we want them displayed
    const SectionType secs[] = {
        SEC_GENOZIP_HEADER,     SEC_RANDOM_ACCESS,      SEC_TXT_HEADER,        SEC_VB_HEADER,
        SEC_CHROM_B250,         SEC_CHROM_DICT,         SEC_POS_B250,          SEC_POS_DICT, 
        SEC_ID_B250,            SEC_ID_DICT,            SEC_VCF_REFALT_B250,   SEC_VCF_REFALT_DICT, 
        SEC_VCF_QUAL_B250,      SEC_VCF_QUAL_DICT,      SEC_VCF_FILTER_B250,   SEC_VCF_FILTER_DICT, 
        SEC_VCF_INFO_B250,      SEC_VCF_INFO_DICT,      SEC_VCF_INFO_SF_B250,  SEC_VCF_INFO_SF_DICT, 
        SEC_VCF_FORMAT_B250,    SEC_VCF_FORMAT_DICT,    SEC_VCF_GT_DATA,       SEC_VCF_FRMT_SF_DICT,
        SEC_HT_DATA,            SEC_STATS_HT_SEPERATOR, SEC_VCF_PHASE_DATA,

        SEC_SEQ_DATA,           SEC_QUAL_DATA,          SEC_SAM_MD_DATA,
        SEC_SAM_BD_DATA,        SEC_SAM_BI_DATA,        SEC_RANDOM_POS_DATA,  
        SEC_NUMERIC_ID_DATA,    SEC_ENST_DATA,          SEC_FASTA_COMMENT_DATA,
        SEC_SAM_QNAME_B250,     SEC_SAM_QNAME_DICT,     SEC_SAM_QNAME_SF_B250,  SEC_SAM_QNAME_SF_DICT,
        SEC_SAM_FLAG_B250,      SEC_SAM_FLAG_DICT,      SEC_SAM_RNAME_B250,     SEC_SAM_RNAME_DICT, 
        SEC_SAM_POS_B250,       SEC_SAM_POS_DICT,
        SEC_SAM_MAPQ_B250,      SEC_SAM_MAPQ_DICT,      SEC_SAM_CIGAR_B250,     SEC_SAM_CIGAR_DICT,     
        SEC_SAM_TLEN_B250,      SEC_SAM_TLEN_DICT,      SEC_SAM_PNEXT_B250,     SEC_SAM_PNEXT_DICT,      
        SEC_SAM_OPTIONAL_B250,  SEC_SAM_OPTIONAL_DICT,  SEC_SAM_OPTNL_SF_B250,  SEC_SAM_OPTNL_SF_DICT,

        SEC_GFF3_SEQID_B250,    SEC_GFF3_SEQID_DICT,    SEC_GFF3_SOURCE_B250,   SEC_GFF3_SOURCE_DICT,
        SEC_GFF3_TYPE_B250,     SEC_GFF3_TYPE_DICT,     SEC_GFF3_START_B250,    SEC_GFF3_START_DICT,
        SEC_GFF3_END_B250,      SEC_GFF3_END_DICT,      SEC_GFF3_SCORE_B250,    SEC_GFF3_SCORE_DICT,
        SEC_GFF3_STRAND_B250,   SEC_GFF3_STRAND_DICT,   SEC_GFF3_PHASE_B250,    SEC_GFF3_PHASE_DICT,
        SEC_GFF3_ATTRS_B250,    SEC_GFF3_ATTRS_DICT,    SEC_GFF3_ATTRS_SF_B250, SEC_GFF3_ATTRS_SF_DICT, 
        
        SEC_FAST_LINEMETA_B250, SEC_FAST_LINEMETA_DICT,
        SEC_FAST_DESC_B250,     SEC_FAST_DESC_DICT,     SEC_FAST_DESC_SF_B250,  SEC_FAST_DESC_SF_DICT
    };

    static const char *categories[] = {
        "Genozip header", "Random access index", "Original file header", "VBlock header", 
        "CHROM b250", "CHROM dict", "POS b250", "POS dict", 
        "ID b250", "ID dict", "REF+ALT b250", "REF+ALT dict", 
        "QUAL b250", "QUAL dict", "FILTER b250", "FILTER dict",
        "INFO names b250", "INFO names dict", "INFO values b250", "INFO values dict", 
        "FORMAT b250", "FORMAT dict", "Other sample info", "FORMAT subfields dict",
        "Haplotype data", "HT separator char", "Phasing char",

        "SEQ data", "QUAL data", "MD data", "BD data", "BI data", "RAND_POS data", 
        "NUMERIC_ID data", "ENST numeric ids", "Comment data",
        "QNAME b250", "QNAME dict", "QNAME subfields b250", "QNAME subfields dict", 
        "FLAG b250", "FLAG dict", "RNAME b250", "RNAME dict", 
        "POS b250 (delta)", "POS dict (delta)", 
        "MAPQ b250", "MAPQ dict", "CIGAR b250", "CIGAR dict", 
        "TLEN b250", "TLEN dict", "PNEXT b250 (delta)", "PNEXT dict (delta)",
        "OPTIONAL names b250", "OPTIONAL names dict", "OPTIONAL values b250", "OPTIONAL values dict",

        "SEQID b250", "SEQID dict", "SOURCE b250", "SOURCE dict", "TYPE b250", "TYPE dict", 
        "START b250", "START dict", "END b250", "END dict", "SCORE b250", "SCORE dict", "STRAND b250", "STRAND dict", 
        "PHASE b250", "PHASE dict",
        "ATTRS names b250", "ATTRS names dict", "ATTRS values b250", "ATTRS values dict",  

        "Line metadata b250", "Line metadata dict",
        "DESC b250", "DESC dict", "DESC subfields b250", "DESC subfields dict"
    };

    unsigned num_secs       = sizeof(secs)/sizeof(secs[0]);
    unsigned num_categories = sizeof(categories)/sizeof(categories[0]);
    ASSERT (num_categories == num_secs, "Error in stats_show_sections: num_categories=%u but num_secs=%u, expecting them to be equal", num_categories, num_secs);

    stats_show_file_metadata();

    char vsize[30], zsize[30], zentries_str[30];

    fprintf (stderr, "Sections stats:\n");
    fprintf (stderr, "                            #Sec          #Entries   %8s          %%        GENOZIP     %%   Ratio\n",
             dt_name (z_file->data_type));
    const char *format = "%22s    %6u  %16s  %9s      %5.1f      %9s %5.1f  %6.1f%s\n";

    int64_t total_txt=0, total_z=0, total_entries=0;
    uint32_t total_sections=0;
    
    unsigned covered[NUM_SEC_TYPES]; memset (covered, 0, sizeof(covered));

    for (unsigned sec_i=0; sec_i < num_secs; sec_i++) {
        int64_t txtbytes  = txt_file->section_bytes[secs[sec_i]];
        int64_t zbytes    = z_file->section_bytes[secs[sec_i]];
        int64_t zentries  = z_file->section_entries[secs[sec_i]];
        int32_t zsections = z_file->num_sections[secs[sec_i]];
        covered[secs[sec_i]]++;;

        if (!txtbytes && !zbytes) continue;

        bool is_dict = section_type_is_dictionary (secs[sec_i]);
        bool is_b250 = section_type_is_b250 (secs[sec_i]);

        int64_t ratio_zbytes;
        if (secs[sec_i] == SEC_VCF_INFO_B250  || 
            secs[sec_i] == SEC_SAM_QNAME_B250 || secs[sec_i] == SEC_SAM_OPTIONAL_B250 ||
            secs[sec_i] == SEC_FAST_DESC_B250 ||
            secs[sec_i] == SEC_GFF3_ATTRS_B250) 
            ratio_zbytes = zbytes + z_file->section_bytes[secs[sec_i+1]] + z_file->section_bytes[secs[sec_i+2]] 
                                  + z_file->section_bytes[secs[sec_i+3]];
        
        else if (is_dict || 
                 secs[sec_i] == SEC_VCF_INFO_SF_B250  || 
                 secs[sec_i] == SEC_SAM_QNAME_SF_B250 || secs[sec_i] == SEC_SAM_OPTNL_SF_B250 ||
                 secs[sec_i] == SEC_FAST_DESC_SF_B250 ||
                 secs[sec_i] == SEC_GFF3_ATTRS_SF_B250)
            ratio_zbytes = 0;
        
        else if (is_b250) 
            ratio_zbytes = zbytes + z_file->section_bytes[secs[sec_i+1]]; // b250 and dict combined
        
        else              
            ratio_zbytes = zbytes;

        int64_t ratio_txtbytes;
        if (secs[sec_i] == SEC_VCF_INFO_B250  || 
            secs[sec_i] == SEC_SAM_QNAME_B250 || secs[sec_i] == SEC_SAM_OPTIONAL_B250 ||
            secs[sec_i] == SEC_FAST_DESC_B250 ||
            secs[sec_i] == SEC_GFF3_ATTRS_B250) 
            ratio_txtbytes = txtbytes + txt_file->section_bytes[secs[sec_i+1]] + txt_file->section_bytes[secs[sec_i+2]] 
                                      + txt_file->section_bytes[secs[sec_i+3]];
        else
            ratio_txtbytes = txtbytes;

        fprintf (stderr, format, categories[sec_i], zsections, 
                 str_uint_commas (zentries, zentries_str),
                 txtbytes ? str_size(txtbytes, vsize) : "       ", 
                 100.0 * (double)txtbytes / (double)txt_file->txt_data_size_single,
                 str_size(zbytes, zsize), 
                 100.0 * (double)zbytes / (double)z_file->disk_size,
                 ratio_zbytes ? (double)ratio_txtbytes / (double)ratio_zbytes : 0,
                 !ratio_zbytes ? (ratio_txtbytes && secs[sec_i] != SEC_FAST_DESC_SF_B250  && 
                                                    secs[sec_i] != SEC_VCF_INFO_SF_B250   && 
                                                    secs[sec_i] != SEC_GFF3_ATTRS_SF_B250 && 
                                                    secs[sec_i] != SEC_SAM_QNAME_SF_B250  && 
                                                    secs[sec_i] != SEC_SAM_OPTNL_SF_B250 ? "\b\b\bInf" : "\b\b\b---") : "");

        total_sections += zsections;
        total_entries  += zentries;
        total_txt      += txtbytes;
        total_z        += zbytes;
    }

    fprintf (stderr, format, "TOTAL", total_sections, str_uint_commas (total_entries, zentries_str),
             str_size(total_txt, vsize), 100.0 * (double)total_txt / (double)txt_file->txt_data_size_single,
             str_size(total_z, zsize),   100.0 * (double)total_z   / (double)z_file->disk_size,
             (double)total_txt / (double)total_z, "");

    uint64_t all_dicts_comp_len=0, all_dicts_uncomp_len=0;
    char s1[20], s2[20], s3[20], s4[20], s5[20], s6[20];

    fprintf (stderr, "\nDictionaries:\n");
    fprintf (stderr, "did_i Name     Type            #Words        #Uniq         Hash    uncomp      comp      comp     comp     %% of \n");
    fprintf (stderr, "                                                                   dict        dict      b250     TOTAL    file \n");
    for (uint32_t i=0; i < z_file->num_dict_ids; i++) { // don't show CHROM-FORMAT as they are already showed above
        uint32_t dict_compressed_size, b250_compressed_size;

        const MtfContext *ctx = &z_file->mtf_ctx[i];

        if (!ctx->mtf_i.len) continue;
        
        sections_get_sizes (ctx->dict_id, &dict_compressed_size, &b250_compressed_size);

        all_dicts_uncomp_len += ctx->dict.len;
        all_dicts_comp_len   += (uint64_t)dict_compressed_size;

        fprintf (stderr, "%-2u    %*.*s %-6.6s %15s %12s %12s %9s %9s %9s %9s %5.1f\n", i, -DICT_ID_LEN, DICT_ID_LEN, err_dict_id (ctx->dict_id), 
                 dict_id_display_type (z_file->data_type, ctx->dict_id), str_uint_commas (ctx->mtf_i.len, s1), str_uint_commas (ctx->mtf.len, s2), 
                 str_uint_commas (ctx->global_hash_prime, s3), str_size (ctx->dict.len, vsize),
                 str_size (dict_compressed_size, s4), str_size (b250_compressed_size, s5),
                 str_size (dict_compressed_size + b250_compressed_size, s6),
                 100.0 * (double)(dict_compressed_size + b250_compressed_size) / (double)total_z);
    }

    fprintf (stderr, "\nTotal dictionaries size: On disk: %s ; In memory: %s\n", 
             str_size (all_dicts_comp_len, s1), str_size (all_dicts_uncomp_len, s2));

    stats_verify_all_covered (covered);

    ASSERTW (total_z == z_file->disk_size, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%s but file size is %s (diff=%d)", 
             str_uint_commas (total_z, s1), str_uint_commas (z_file->disk_size, s2), (int32_t)(z_file->disk_size - total_z));

    // note: we use txt_data_so_far_single and not txt_data_size_single, because the latter has estimated size if disk_size is 
    // missing, while txt_data_so_far_single is what was actually processed
    ASSERTW (total_txt == txt_file->txt_data_so_far_single || flag_optimize, "Hmm... incorrect calculation for %s sizes: total section sizes=%s but file size is %s (diff=%d)", 
             dt_name (z_file->data_type), str_uint_commas (total_txt, s1), str_uint_commas (txt_file->txt_data_size_single, s2), 
             (int32_t)(txt_file->txt_data_so_far_single - total_txt)); 
}

void stats_show_content (void)
{
    stats_show_file_metadata();

    char vsize[30], zsize[30];
    int64_t total_txt=0, total_z=0;

    fprintf (stderr, "Compression stats:\n");
    fprintf (stderr, "                              Original %%       GENOZIP     %%  Ratio\n");
    const char *format = "%22s   %8s %5.1f      %9s %5.1f  %5.1f%s\n";

    const char *categories[] = {"Haplotype data", "Other sample data", 
                                "SEQ data", "QUAL data",
                                "Header & other fields", "Index" };

#define NUM_CATEGORIES 6

    int sections_per_category[NUM_CATEGORIES][100] = { 
        { SEC_HT_DATA ,  NIL },
        { SEC_VCF_PHASE_DATA, SEC_VCF_GT_DATA, SEC_VCF_FRMT_SF_DICT, SEC_STATS_HT_SEPERATOR, NIL},
        { SEC_SEQ_DATA,  NIL },
        { SEC_QUAL_DATA, NIL },
        { SEC_TXT_HEADER, SEC_VB_HEADER, 
          
          SEC_CHROM_B250, SEC_POS_B250, SEC_ID_B250, SEC_VCF_REFALT_B250, 
          SEC_VCF_QUAL_B250, SEC_VCF_FILTER_B250, SEC_VCF_INFO_B250, SEC_VCF_FORMAT_B250, SEC_VCF_INFO_SF_B250, 
          SEC_CHROM_DICT, SEC_POS_DICT, SEC_ID_DICT, SEC_VCF_REFALT_DICT, SEC_VCF_QUAL_DICT,
          SEC_VCF_FILTER_DICT, SEC_VCF_INFO_DICT, SEC_VCF_INFO_SF_DICT, SEC_VCF_FORMAT_DICT,       

          SEC_NUMERIC_ID_DATA,    SEC_ENST_DATA,
          
          SEC_RANDOM_POS_DATA,  SEC_SAM_MD_DATA,        SEC_SAM_BD_DATA,       SEC_SAM_BI_DATA, 
          SEC_SAM_QNAME_B250,     SEC_SAM_QNAME_DICT,     SEC_SAM_QNAME_SF_B250, SEC_SAM_QNAME_SF_DICT,
          SEC_SAM_FLAG_B250,      SEC_SAM_FLAG_DICT,      SEC_SAM_RNAME_B250,    SEC_SAM_RNAME_DICT, 
          SEC_SAM_POS_B250,       SEC_SAM_POS_DICT,       SEC_SAM_MAPQ_B250,     SEC_SAM_MAPQ_DICT,      
          SEC_SAM_CIGAR_B250,     SEC_SAM_CIGAR_DICT,     SEC_SAM_TLEN_B250,     SEC_SAM_TLEN_DICT,      
          SEC_SAM_PNEXT_B250,     SEC_SAM_PNEXT_DICT, 
          SEC_SAM_OPTIONAL_B250,  SEC_SAM_OPTIONAL_DICT,  SEC_SAM_OPTNL_SF_B250, SEC_SAM_OPTNL_SF_DICT,
          
          SEC_FAST_DESC_B250,     SEC_FAST_DESC_DICT,     SEC_FAST_DESC_SF_B250,  SEC_FAST_DESC_SF_DICT,  
          SEC_FAST_LINEMETA_B250, SEC_FAST_LINEMETA_DICT, SEC_FASTA_COMMENT_DATA,

          SEC_GFF3_SEQID_B250,    SEC_GFF3_SEQID_DICT,    SEC_GFF3_SOURCE_B250,   SEC_GFF3_SOURCE_DICT,
          SEC_GFF3_TYPE_B250,     SEC_GFF3_TYPE_DICT,     SEC_GFF3_START_B250,    SEC_GFF3_START_DICT,
          SEC_GFF3_END_B250,      SEC_GFF3_END_DICT,      SEC_GFF3_SCORE_B250,    SEC_GFF3_SCORE_DICT,
          SEC_GFF3_STRAND_B250,   SEC_GFF3_STRAND_DICT,   SEC_GFF3_PHASE_B250,    SEC_GFF3_PHASE_DICT,
          SEC_GFF3_ATTRS_B250,    SEC_GFF3_ATTRS_DICT,    SEC_GFF3_ATTRS_SF_B250, SEC_GFF3_ATTRS_SF_DICT, 

          NIL },
        { SEC_RANDOM_ACCESS, SEC_GENOZIP_HEADER, NIL }
    };

    unsigned covered[NUM_SEC_TYPES]; memset (covered, 0, sizeof(covered));

    for (unsigned i=0; i < NUM_CATEGORIES; i++) {

        int64_t txtbytes=0, zbytes=0;
        for (int *sec_i = sections_per_category[i]; *sec_i != NIL; sec_i++) {
            txtbytes += txt_file->section_bytes[*sec_i];
            zbytes += z_file->section_bytes[*sec_i];
            covered[*sec_i]++;
        }

        if (!txtbytes && !zbytes) continue;

        fprintf (stderr, format, categories[i], 
                 str_size(txtbytes, vsize), 100.0 * (double)txtbytes / (double)txt_file->txt_data_size_single,
                 str_size(zbytes, zsize), 100.0 * (double)zbytes / (double)z_file->disk_size,
                 zbytes ? (double)txtbytes / (double)zbytes : 0,
                 !zbytes ? (txtbytes ? "\b\b\bInf" : "\b\b\b---") : "");

        total_txt      += txtbytes;
        total_z        += zbytes;
    }

    fprintf (stderr, format, "TOTAL", 
             str_size(total_txt, vsize), 100.0 * (double)total_txt / (double)txt_file->txt_data_size_single,
             str_size(total_z, zsize),   100.0 * (double)total_z   / (double)z_file->disk_size,
             (double)total_txt / (double)total_z, "");

    stats_verify_all_covered (covered);

    ASSERTW (total_z == z_file->disk_size, "Hmm... incorrect calculation for GENOZIP sizes: total section sizes=%"PRId64" but file size is %"PRId64" (diff=%"PRId64")", 
             total_z, z_file->disk_size, z_file->disk_size - total_z);

    ASSERTW (total_txt == txt_file->txt_data_size_single || flag_optimize, "Hmm... incorrect calculation for %s sizes: total section sizes=%"PRId64" but file size is %"PRId64" (diff=%"PRId64")", 
             dt_name (z_file->data_type), total_txt, txt_file->txt_data_size_single, txt_file->txt_data_size_single - total_txt);

}

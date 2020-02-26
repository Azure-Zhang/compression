// ------------------------------------------------------------------
//   sections.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef SECTIONS_INCLUDED
#define SECTIONS_INCLUDED

#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <stdbool.h>
#else
#include "compatability/visual_c_stdint.h"
#include "compatability/visual_c_stdbool.h"
#endif

#include "md5.h"
#include "dict_id.h"

// note: the numbering of the sections cannot be modified, for backward compatability
typedef enum {
    // data sections - statring in v1
    SEC_VCF_HEADER         = 0,  SEC_VB_HEADER           = 1, 
    SEC_GENOTYPE_DICT      = 2,  SEC_GENOTYPE_DATA       = 3, 
    SEC_PHASE_DATA         = 4,  SEC_HAPLOTYPE_DATA      = 5,

    // data sections added in v2
    SEC_GENOZIP_HEADER     = 6,   SEC_RANDOM_ACCESS      = 7,

    SEC_CHROM_DICT         = 8,   SEC_CHROM_B250         = 9,
    SEC_POS_DICT           = 10,  SEC_POS_B250           = 11,
    SEC_ID_DICT            = 12,  SEC_ID_B250            = 13,  
    SEC_REFALT_DICT        = 14,  SEC_REFALT_B250        = 15,  
    SEC_QUAL_DICT          = 16,  SEC_QUAL_B250          = 17, 
    SEC_FILTER_DICT        = 18,  SEC_FILTER_B250        = 19, 
    SEC_INFO_DICT          = 20,  SEC_INFO_B250          = 21, 
    SEC_FORMAT_DICT        = 22,  SEC_FORMAT_B250        = 23,
    SEC_INFO_SUBFIELD_DICT = 24,  SEC_INFO_SUBFIELD_B250 = 25,
    
    // These sections are not real sections - they don't appear in the genozip file - just for stats. They can be changed if needed.
    SEC_STATS_HT_SEPERATOR, 
} SectionType;

// we put the names here in a #define so we can eyeball their identicality to SectionType
#define SECTIONTYPE_NAMES { \
    "SEC_VCF_HEADER"        ,  "SEC_VB_HEADER",\
    "SEC_GENOTYPE_DICT"     ,  "SEC_GENOTYPE_DATA",\
    "SEC_PHASE_DATA"        ,  "SEC_HAPLOTYPE_DATA",\
    \
    "SEC_GENOZIP_HEADER"    ,  "SEC_RANDOM_ACCESS",\
    \
    "SEC_CHROM_DICT"        ,  "SEC_CHROM_B250",\
    "SEC_POS_DICT"          ,  "SEC_POS_B250",\
    "SEC_ID_DICT"           ,  "SEC_ID_B250",\
    "SEC_REFALT_DICT"       ,  "SEC_REFALT_B250",\
    "SEC_QUAL_DICT"         ,  "SEC_QUAL_B250",\
    "SEC_FILTER_DICT"       ,  "SEC_FILTER_B250",\
    "SEC_INFO_DICT"         ,  "SEC_INFO_B250",\
    "SEC_FORMAT_DICT"       ,  "SEC_FORMAT_B250",\
    "SEC_INFO_SUBFIELD_DICT",  "SEC_INFO_SUBFIELD_B250",\
    \
    "SEC_STATS_HT_SEPERATOR" \
}

#define NUM_SEC_TYPES (SEC_STATS_HT_SEPERATOR+1) // put this here and not in sections.h as its used in vb.h that is widely used

#define section_type_is_dictionary(s) (((s) >= SEC_CHROM_DICT && (s) <= SEC_FORMAT_DICT && (s) % 2 == SEC_CHROM_DICT % 2) ||       \
                                        (s) == SEC_INFO_SUBFIELD_DICT || \
                                        (s) == SEC_GENOTYPE_DICT)

#define section_type_is_b250(s)       (((s) >= SEC_CHROM_B250 && (s) <= SEC_FORMAT_B250 && (s) % 2 == SEC_CHROM_B250 % 2) ||       \
                                        (s) == SEC_INFO_SUBFIELD_B250)

#define section_type_is_vb(s)         (((s) >= SEC_CHROM_DICT && (s) <= SEC_INFO_SUBFIELD_B250) || \
                                       ((s) >= SEC_VB_HEADER && (s) <= SEC_HAPLOTYPE_DATA))

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

#pragma pack(push, 1) // structures that are part of the genozip format are packed.

// section headers are encoded in Big Endian (see https://en.wikipedia.org/wiki/Endianness)
// the reason for selecting big endian is that I am developing on little endian CPU (Intel) so
// endianity bugs will be discovered more readily this way

typedef struct {
    uint32_t magic; 
    uint32_t compressed_offset;     // number of bytes from the start of the header that is the start of compressed data (sizeof header + header encryption padding)
    uint32_t data_encrypted_len;    // = data_compressed_len + padding if encrypted, 0 if not
    uint32_t data_compressed_len;
    uint32_t data_uncompressed_len;
    uint32_t variant_block_i;       // VB with in file starting from 1 ; 0 for VCF Header
    uint16_t section_i;             // section within VB - 0 for Variant Data
    uint8_t  section_type;          
    uint8_t  flags;                 // section-type specific flags
} SectionHeader; 

// data types genozip can compress
#define DATA_TYPE_VCF 0

#define ENC_TYPE_NONE   0
#define ENC_TYPE_AES256 1

// Note: see v1 sections in v1.c

typedef struct {
    SectionHeader h;
    uint8_t genozip_version;
    uint8_t encryption_type;          // one of ENC_TYPE_*
    uint16_t data_type;               // one of DATA_TYPE_*
    uint32_t num_samples;             // number of samples. "samples" is data_type-dependent. 
    uint64_t uncompressed_data_size;  // data size of uncompressed file, if uncompressed as a single file
    uint64_t num_items_concat;        // number of items in a concatenated file. "item" is data_type-dependent. For VCF, it is lines.
    uint32_t num_sections;            // number sections in this file (including this one)
    uint32_t num_vcf_components;      // number of vcf concatenated components in this file (1 if no concatenation)
    uint32_t num_info_dictionary_sections;  // DO WE NEED THIS?
    uint32_t num_gt_dictionary_sections;    // AND THIS?

    Md5Hash md5_hash_concat;   // md5 of original VCF file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.

    uint8_t password_test[16]; // short encrypted block - used to test the validy of a password
#define FILE_METADATA_LEN 72
    char created[FILE_METADATA_LEN];    
} 
SectionHeaderGenozipHeader;

// this footer appears AFTER the genozip header data, facilitating reading the genozip header in reverse from the end of the file
typedef struct {
    uint64_t genozip_header_offset;
    uint32_t magic;
} SectionFooterGenozipHeader;

#define COMPRESSION_TYPE_NONE  0
#define COMPRESSION_TYPE_GZIP  1
#define COMPRESSION_TYPE_BZIP2 2
#define COMPRESSION_TYPE_BGZIP 3

// The VCF header section appears once in the file (or multiple times in case of concatenation), and includes the VCF file header 
typedef struct {
    SectionHeader h;
    uint64_t vcf_data_size;    // number of bytes in the original VCF file. Concat mode: entire file for first SectionHeaderVCFHeader, and only for that VCF if not first
#define NUM_LINES_UNKNOWN ((uint64_t)-1) 
    uint64_t num_lines;        // number of variants (data lines) in the original VCF file. Concat mode: entire file for first SectionHeaderVCFHeader, and only for that VCF if not first
    uint32_t num_samples;      // number of samples in the original VCF file
    uint32_t max_lines_per_vb; // log2 of the upper bound on how many variants (data lines) a VB can have in this file
    uint8_t compression_type;  // compression type of original file, one of COMPRESSION_TYPE_*
    Md5Hash md5_hash_single;   // non-0 only if this genozip file is a result of concatenatation with --md5. md5 of original single VCF file.

#define VCF_FILENAME_LEN 256
    char vcf_filename[VCF_FILENAME_LEN];    // filename of this single component. without path, 0-terminated.

} SectionHeaderVCFHeader; 

// The variant data section appears for each variant block

// Note: the reason we index the haplotypes across all samples rather than for each sample group, despite more bits in haplotype_index entry
// is that the better global sorting results in overall smaller file. 
// 
// Tested with the first 1024 variant block of chr22 of 1000-genomes (5096 haplotypes):
//   - Sorting all haplotypes together:                                        62.8KB Index 13bit x 5096 = 8.3K Total: 71.1K <--- better sort together despite bigger index
//   - Sorting each 1024 haplotypes (5 sample groups) and sorting separately - 68.1KB Index 10bit x 5096 = 6.4K Total: 74.5K
//
// Note: this doesn't affect retrieval time of a specific sample, because the sample block is just looked up via the index
 
typedef struct {
    SectionHeader h;
    uint32_t first_line;               // line (starting from 1) of this variant block in the VCF file
                                       // new in v2: if this value is 0, then this is the terminating section of the file. after it is either EOF or a VCF Header section of the next concatenated file
    uint32_t num_lines;                // number of variants in this block
    uint8_t phase_type;
    
    // flags
    uint8_t has_genotype_data : 1;     // 1 if there is at least one variant in the block that has FORMAT with have anything except for GT 
    uint8_t for_future_use    : 7;
    uint16_t unused1;                  // new in v2: padding / ffu

    // features of the data
    uint32_t num_samples;
    uint32_t num_haplotypes_per_line;  // 0 if no haplotypes
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block;
    uint16_t ploidy;
    uint8_t  unused2;                  // new in v2: padding / ffu
/*    uint8_t  field_dictionary_sections_bitmap; // new in v2: bitmap on the existance of field dictionaries: from LSb=CHROM to MSb=FORMAT
                                       // note: each section corresponds to a dict_id, BUT not all dict_ids have sections - 
                                       // they won't, if all their snips have already been entered to the dictionary by previous VBs
    uint32_t num_info_dictionary_sections;  // new in v2
    uint32_t num_gt_dictionary_sections;    // rename from num_dictionary_sections in v1
*/    uint32_t num_dict_ids;             // v2 change: uint16_t->uint32_t. number of dict_ids used in this VB. the actual fields are deduced from the FORMAT column on each line
    uint32_t num_info_subfields;       // v2 addition: number INFO subfields present in this VB. each subfield has a dictionary. the dictionary for each subfield is deduced from the INFO names, in the order the appear in the VB. Eg. the first INFO name, is subfield=0, and its dictionary is looked up by the name.
    uint32_t max_gt_line_len;

    uint32_t vb_data_size;             // size of variant block as it appears in the source file
    uint32_t z_data_bytes;             // total bytes of this variant block in the genozip file including all sections and their headers
    uint16_t haplotype_index_checksum;
    uint16_t unused3;                  // new in v2: padding / ffu

    // uint8_t haplotype_index[];      // removed in V2 - haplotype_index is now the data of this section rather than part of the header
} SectionHeaderVbHeader; 

typedef struct {
    SectionHeader h;                   // in v1, the section_type was SEC_DICTIONARY, in v2 is is the specific SEC_*_DICT
    uint32_t num_snips;                // number of items in dictionary
    DictIdType dict_id;           
} SectionHeaderDictionary; 

// b250 data encoded according to a dictionary
typedef struct {
    SectionHeader h;
    uint32_t num_b250_items;           // number of items in b250 items
    uint8_t encoding;                  // one of Base250Encoding
    uint8_t unused[3];
    DictIdType dict_id;           
} SectionHeaderBase250;     

// the data of SEC_SECTION_LIST is an array of the following type, as is the z_file->section_list_buf
typedef struct {
    uint64_t offset;              // offset of this section in the file
    DictIdType dict_id;           // used if this section is a DICT or a B250 section
    uint32_t variant_block_i;
    uint8_t section_type;
    uint8_t include         : 1;  // include this variant block in zcat
    uint8_t for_future_use  : 7;
} SectionListEntry;

extern void sections_add_to_list (VariantBlockP pseudo_vb, const SectionHeader *header, uint64_t offset);
extern void sections_update_list (VariantBlockP vb, bool is_z_data);
extern void BGEN_sections_list (BufferP sections_list_buf);
extern const char *st_name (unsigned sec_type);

// ------------------------------------------------------------------------------------------------------
// GENOZIP_FILE_FORMAT_VERSION==1 historical version - we support uncomrpessing old version files

typedef struct {
    SectionHeader h;
    uint8_t  genozip_version;
    uint32_t num_samples;    // number of samples in the original VCF file
    uint64_t vcf_data_size;  // number of bytes in the original VCF file. Concat mode: entire file for first SectionHeaderVCFHeader, and only for that VCF if not first
    uint64_t num_lines;      // number of variants (data lines) in the original VCF file. Concat mode: entire file for first SectionHeaderVCFHeader, and only for that VCF if not first
    Md5Hash md5_hash_concat; // md5 of original VCF file, or 0s if no hash was calculated. if this is a concatenation - this is the md5 of the entire concatenation.
    Md5Hash md5_hash_single; // non-0 only if this genozip file is a result of concatenatation with --md5. md5 of original single VCF file.

#define v1_VCF_FILENAME_LEN 256
    char vcf_filename[v1_VCF_FILENAME_LEN];    // filename of this single component. without path, 0-terminated.
    
#define v1_FILE_METADATA_LEN 72
    char created[v1_FILE_METADATA_LEN];    
} v1_SectionHeaderVCFHeader; 

typedef struct {
    SectionHeader h;
    uint32_t first_line;               // line (starting from 1) of this variant block in the VCF file
    uint32_t num_lines;                // number of variants in this block
    uint8_t phase_type;
    
    // flags
    uint8_t has_genotype_data : 1;     // 1 if there is at least one variant in the block that has FORMAT with have anything except for GT 
    uint8_t is_sorted_by_pos  : 1;     // 1 if variant block is sorted by POS
    uint8_t for_future_use    : 6;

    // features of the data
    uint32_t num_samples;
    uint32_t num_haplotypes_per_line;  // 0 if no haplotypes
    uint32_t num_sample_blocks;
    uint32_t num_samples_per_block;
    uint16_t ploidy;
    uint16_t num_dict_ids;            // number of subfields used in this VB. the actual fields are deduced from the FORMAT column on each line
    uint16_t num_dictionary_sections;  // we have a dictionary for each dict_id, in which new words were introduced in this VB that weren't already in the dictionary. so num_dictionary_sections <= num_dict_ids
    uint32_t max_gt_line_len;

    char chrom[MAX_CHROM_LEN];         // a null-terminated ID of the chromosome
    int64_t min_pos, max_pos;          // minimum and maximum POS values in this VB. -1 if unknown. Note: our format support 64bit POS, but VCF specification as well as the POS dictionary supports 32bit (values 0 to 2M-1)

    uint32_t vb_data_size;             // size of variant block as it appears in the source file
    uint32_t z_data_bytes;             // total bytes of this variant block in the genozip file including all sections and their headers
    uint16_t haplotype_index_checksum;
    uint8_t haplotype_index[];         // length is num_haplotypes. e.g. the first entry shows for the first haplotype in the original file, its index into the permuted block. # of bits per entry is roundup(log2(num_samples*ploidy)).
} v1_SectionHeaderVariantData; 

// ------------------------------------------------------------------------------------------------------


#pragma pack(pop)

#endif

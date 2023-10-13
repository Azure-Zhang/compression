// ------------------------------------------------------------------
//   sections.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "digest.h"
#include "local_type.h"

// this data must be perfectly aligned with SectionType.
#define SECTIONTYPE_ABOUT {    \
    {"SEC_RANDOM_ACCESS",   sizeof (SectionHeader)              }, \
    {"SEC_REFERENCE",       sizeof (SectionHeaderReference)     }, \
    {"SEC_REF_IS_SET",      sizeof (SectionHeaderReference)     }, \
    {"SEC_REF_HASH",        sizeof (SectionHeaderRefHash)       }, \
    {"SEC_REF_RAND_ACC",    sizeof (SectionHeader)              }, \
    {"SEC_REF_CONTIGS",     sizeof (SectionHeader)              }, \
    {"SEC_GENOZIP_HEADER",  sizeof (SectionHeaderGenozipHeader) }, \
    {"SEC_DICT_ID_ALIASES", sizeof (SectionHeader)              }, \
    {"SEC_TXT_HEADER",      sizeof (SectionHeaderTxtHeader)     }, \
    {"SEC_VB_HEADER",       sizeof (SectionHeaderVbHeader)      }, \
    {"SEC_DICT",            sizeof (SectionHeaderDictionary)    }, \
    {"SEC_B250",            sizeof (SectionHeaderCtx)           }, \
    {"SEC_LOCAL",           sizeof (SectionHeaderCtx)           }, \
    {"SEC_CHROM2REF_MAP",   sizeof (SectionHeader)              }, \
    {"SEC_STATS",           sizeof (SectionHeader)              }, \
    {"SEC_BGZF",            sizeof (SectionHeader)              }, \
    {"SEC_RECON_PLAN",      sizeof (SectionHeaderReconPlan)     }, \
    {"SEC_COUNTS",          sizeof (SectionHeaderCounts)        }, \
    {"SEC_REF_IUPACS",      sizeof (SectionHeader)              }, \
}

#define HEADER_IS(_st) (header->section_type == SEC_##_st)
#define ST(_st) (st == SEC_##_st)

// Section headers - big endian

#define GENOZIP_MAGIC 0x27052012

#pragma pack(1) // structures that are part of the genozip format are packed.

// section headers are encoded in Big Endian (see https://en.wikipedia.org/wiki/Endianness)
// the reason for selecting big endian is that I am developing on little endian CPU (Intel) so
// endianity bugs will be discovered more readily this way

typedef enum __attribute__ ((__packed__)) { STORE_NONE, STORE_INT, STORE_FLOAT, STORE_INDEX        } StoreType; // values for SectionFlags.ctx.store
typedef enum __attribute__ ((__packed__)) { B250_BYTES_4, B250_BYTES_3, B250_BYTES_2, B250_BYTES_1 } B250Size; // part of the file format - goes into SectionHeaderCtx.b250_size

typedef enum __attribute__ ((__packed__)) { SAG_NONE, SAG_BY_SA, SAG_BY_NH, SAG_BY_SOLO, SAG_BY_CC, SAG_BY_FLAG, NUM_SAG_TYPES } SagType;
#define SAM_SAG_TYPE_NAMES                { "NONE",   "BY_SA",   "BY_NH",   "BY_SOLO",   "BY_CC",   "BY_FLAG" }
#define IS_SAG_SA   (segconf.sag_type == SAG_BY_SA)
#define IS_SAG_NH   (segconf.sag_type == SAG_BY_NH)
#define IS_SAG_SOLO (segconf.sag_type == SAG_BY_SOLO)
#define IS_SAG_CC   (segconf.sag_type == SAG_BY_CC)
#define IS_SAG_FLAG (segconf.sag_type == SAG_BY_FLAG)

// goes into SectionHeader.flags and also SectionEnt.flags
typedef union SectionFlags {  
    uint8_t flags;

    struct FlagsGenozipHeader {
        // note: if updating dts_* flags, update in zfile_compress_genozip_header, sections_show_header too
        #define dts_ref_internal dt_specific // SAM: REF_INTERNAL was used for compressing (i.e. SAM file without reference) (introduced v6)
        #define v14_dts_paired   dt_specific // FASTQ: This z_file contains one or more pairs of FASTQs compressed with --pair (introduced v9.0.13 until v14, since v15 moved to TxtHeader)
        #define dts_mismatch     dt_specific // CHAIN: This chain file's contigs mismatch its references, so it cannot be used with --chain (introduced v12.0.35)
        uint8_t dt_specific      : 1;  // this flag has a different meaning depending on the data_type, may be one of the above ^ 
        uint8_t aligner          : 1;  // SAM, FASTQ: our aligner may have used to align sequences to the reference (always with FASTQ if compressed with a reference, sometimes with SAM)
        uint8_t txt_is_bin       : 1;  // BAM: Source file is binary 
        uint8_t has_digest       : 1;  // true if compressed with digest (i.e. not --optimize etc) (v15) 
        #define v14_bgzf has_digest    // Up to v14: Reconstruct as BGZF (user may override) (determined by the last component) (since v15, determined by SectionHeaderTxtHeader.src_codec)
        uint8_t adler            : 1;  // true if Adler32 is used, false if MD5 is used 
        uint8_t has_gencomp      : 1;  // VCF: file supports dual coordinates - last two components are the "liftover rejects" data (v12)
                                       // SAM/BAM: PRIM and/or DEPN components exist (v14)
        uint8_t has_taxid        : 1;  // each line in the file has Taxonomic ID information (v12)
        #define dts2_deep        dt_specific2 // SAM: this file is compressed with --deep (v15)
        uint8_t dt_specific2     : 1;  // this flag has a different meaning depending on the data_type, may be one of the above ^ (v15, "unused" up to v14)
    } genozip_header;

    struct FlagsTxtHeader {
        uint8_t pair             : 2;  // FASTQ component (inc. in a Deep SAM): PAIR_R1 or PAIR_R2 if this component is a paired FASTQ component (v15)
        #define v13_dvcf_comp_i pair   // v12-13: DVCF: 0=Main 1=Primary-only rejects 2=Luft-only rejects (in v14, this moved to SectionEnt.comp_i)
        uint8_t is_txt_luft      : 1;  // VCF: true if original source file was a dual-coordinates file in Luft rendition (v12)
        uint8_t unused           : 5;
    } txt_header;

    union FlagsVbHeader {
        struct FlagsVbHeaderVcf {
            uint8_t coords             : 2; // DC_PRIMARY if it contains TOPLEVEL container, DC_LUFT if LUFT toplevel container, or DC_BOTH if both (DC_NONE prior to v12)
            uint8_t use_null_DP_method : 1; // "null_DP" method is used (13.0.5)
            uint8_t unused             : 5;
        } vcf;
        struct FlagsVbHeaderSam {
            uint8_t unused             : 2; // for now, reserved for coords, should we have a dual-coord SAM in the future
            uint8_t v13_is_sorted      : 1; // Likely sorted (copied from segconf.is_sorted)     - introduced 13.0.3 and canceled v14
            uint8_t v13_is_collated    : 1; // Likely collated (copied from segconf.is_collated) - introduced 13.0.3 and canceled v14
            uint8_t unused2            : 4;
        } sam;
        struct FlagsVbHeaderGff {
            uint8_t embedded_fasta     : 1; // this VB consists embedded FASTA and has data_type=DT_FASTA (v15)
            uint8_t unused             : 7;
        } gff;
    } vb_header;

    struct FlagsBgzf {
        uint8_t has_eof_block    : 1;
        uint8_t level            : 4;  // 0-12 for libdeflate or 0-9 for zlib level: 15 means unknown
        BgzfLibraryType library  : 3;  // ignored if level=15 (introduced 9.0.16)
    } bgzf;

    struct FlagsCtx {
        StoreType store          : 2;  // after reconstruction of a snip, store it in ctx.last_value
        uint8_t paired           : 1;  // FASTQ VB (inc. Deep): R1 context: pair-identical: this section should be used for both R1 and R2 if the corresponding R2 section missing (v15)
                                       //                       R2 context: pair-assisted: reconstruction of this context requires loading the corresponding R1 context to ctx->b250R1/localR1
        uint8_t store_delta      : 1;  // introduced v12.0.41: after reconstruction of a snip, store last_delta. notes: 1. last_delta also stored in case of a delta snip. 2. if using this, store=STORE_INT must be set.
        #define v13_copy_local_param spl_custom   // up to v13: copy ctx.b250/local.param from SectionHeaderCtx.param. since v14, piz always copies, except for LT_BITMAP
        uint8_t spl_custom       : 1;  // introduced v14: similar to store_per_line, but storing is done by the context's SPECIAL function, instead of in reconstruct_store_history
        uint8_t all_the_same     : 1;  // SEC_B250: the b250 data contains only one element, and should be used to reconstruct any number of snips from this context

        #define same_line        ctx_specific_flag // v13.0.5: Valid for contexts that use SNIP_OTHER_DELTA and SNIP_DIFF: if true, reconstructs gets the value in the line (whether before or after). if false, it gets the last value.
        #define no_textual_seq   ctx_specific_flag // v14.0.0: SAM_SQBITMAP: indicates that sam_piz_sam2bam_SEQ doesn't need to store textual_seq. 
        #define depn_clip_hard   ctx_specific_flag // v14.0.0: OPTION_SA_Z in SAM_COMP_MAIN: if true: depn lines, if their CIGAR has a clipping, it is hard clipping (H)
        #define lookback0_ok     ctx_specific_flag // v14.0.0: contexts that are items of a container with lookback. indicates that a SNIP_LOOKBACK when lookback=0 is not an error.
        #define trailing_zero_ok ctx_specific_flag // v14.0.17: INFO_QD: whether reconstruct with or without trailing zeros. 
        #define acgt_no_x        ctx_specific_flag // v15.0.13: contexts compressed with CODEC_ACGT: (NONREF and others): this context has no _X companion context
        uint8_t ctx_specific_flag: 1;  // v10.0.3: flag specific a context 

        uint8_t store_per_line   : 1;  // v12.0.41: store value or text for each line - in context->history        
    } ctx;

    struct FlagsDict {
        uint8_t deep_sam         : 1;  // v15: Deep: dictionary used by a SAM   
        uint8_t deep_fastq       : 1;  // v15: Deep: dictionary used by a FASTQ 
        uint8_t unused           : 6;
    } dictionary;

    struct FlagsRandomAccess {
        uint8_t luft             : 1;  // this section is for reconstructing the file in the Luft coordinates (introduced v12)
        uint8_t unused           : 7;
    } random_access;
    
    struct FlagsReconPlan {
        uint8_t luft             : 1;  // this section is for reconstructing the file in the Luft coordinates (introduced v12)
        uint8_t frag_len_bits    : 4;  // 2^(frag_len_bits+17) is the maximum fragment length (v14) 
        uint8_t unused           : 3;
    } recon_plan;

    struct FlagsRefContigs {
        uint8_t sequential_ref_index : 1;  // v14: first ref_contig has ref_index 0, second has 1 etc. ref_index set to 0 on disk and should be set by piz.
        uint8_t compacted            : 1;  // v14.0.10: for stored reference (not reference file): contigs stored as CompactContig
        uint8_t unused               : 6;
    } ref_contigs;

} SectionFlags __attribute__((__transparent_union__));

typedef struct FlagsBgzf FlagsBgzf;

#define SECTION_FLAGS_NONE ((SectionFlags){ .flags = 0 })

typedef struct SectionHeader {         // 28 bytes
    union {
    uint32_t     magic; 
    uint32_t     uncomp_adler32;       // used in --verify-codec (not in file format)
    };
    union {
    uint32_t     v14_compressed_offset;// up to v14: number of bytes from the start of the header that is the start of compressed data (sizeof header + header encryption padding) (32b up to v14, but upper 16b always 0)
    uint32_t     z_digest;             // Adler32 of the compressed/encrypted data of this section (v15)
    };
    uint32_t     data_encrypted_len;   // = data_compressed_len + padding if encrypted, 0 if not
    uint32_t     data_compressed_len;
    uint32_t     data_uncompressed_len;
    uint32_t     vblock_i;             // VB with in file starting from 1 ; 0 for non-VB sections
    SectionType  section_type;         // 1 byte
    Codec        codec;                // 1 byte - primary codec in which this section is compressed
    Codec        sub_codec;            // 1 byte - sub codec, in case primary codec invokes another codec
    SectionFlags flags;                // 1 byte
} SectionHeader, *SectionHeaderP; 

typedef struct {
    SectionHeader;
    uint8_t  genozip_version;
    EncryptionType encryption_type;    // one of EncryptionType
    uint16_t data_type;                // one of DataType
    uint64_t recon_size_prim;          // data size of reconstructed file, if uncompressing as a single file in primary coordinates
    uint64_t num_lines_bound;          // number of lines in a bound file. "line" is data_type-dependent. For FASTQ, it is a read.
    uint32_t num_sections;             // number sections in this file (including this one)
    union {
        struct {                       // v14
            uint16_t vb_size;          // segconf.vb_size in MB (if not exact MB in zip - rounded up to nearest MB)
            char unused;                       
            CompIType num_txt_files;   // number of txt bound components in this file (1 if no binding). We don't count generated components (gencomp).
        };
        uint32_t v13_num_components;
    };
    union { // 16 bytes
        Digest genome_digest;          // DT_REF: Digest of genome as loaded to memory (algorithm determined by genozip_header.adler) (v15)
        Digest v14_REF_fasta_md5;      // DT_REF: MD5 of original FASTA file (v14, buggy)
        Digest FASTQ_v13_digest_bound; // DT_FASTQ: up to v13 "digest_bound": digest of concatenated pair of FQ (regarding other bound files in v13 - DVCF has digest 0, and other bound files are not reconstructable with v14+)    
    };
    #define PASSWORD_TEST "WhenIThinkBackOnAllTheCrapIlearntInHighschool"
    uint8_t  password_test[16];        // short encrypted block - used to test the validy of a password
#define FILE_METADATA_LEN 72
    char     created[FILE_METADATA_LEN];  // nul-terminated metadata
    Digest   license_hash;             // MD5(license_num)
#define REF_FILENAME_LEN 256
    union {
    char ref_filename[REF_FILENAME_LEN];   // external reference filename, nul-terimated. ref_filename[0]=0 if there is no external reference. DT_CHAIN: LUFT reference filename.
    char fasta_filename[REF_FILENAME_LEN]; // DT_REF: fasta file used to generated this reference
    }; 
    union {
        Digest ref_genome_digest;      // Uses REF_EXTERNAL: SectionHeaderGenozipHeader.genome_digest of the reference file
        Digest refhash_digest;         // DT_REF: digest of refhash (v15)
    };
    union { // 272 bytes - data-type specific
        struct {
            char prim_filename[REF_FILENAME_LEN]; // external primary coordinates reference file, nul-terimated. added v12.
            Digest prim_genome_digest; // SectionHeaderGenozipHeader.genome_digest of the primary reference file. added v12.
            char unused[0];
        } chain;
        struct {
            // copied from their respective values in segconf, and copied back to segconf in PIZ
            DictId segconf_seq_len_dict_id;   // SAM: dict_id of one of the Q?NAME contexts (the "length=" item), which is expected to hold the seq_len for this read. 0 if there is no such item. v14.
            uint32_t segconf_seq_len;         // SAM: "standard" seq_len (only applicable to some short-read files). v14.
            SagType segconf_sag_type;         // SAM: v14
            uint8_t segconf_seq_len_cm;       // SAM: v14   
            uint8_t segconf_ms_type      : 3; // SAM: v14 
            uint8_t segconf_has_MD_or_NM : 1; // SAM: PIZ should call sam_analyze_copied_SEQ for dependent reads unless explicitly told not to. v14.
            uint8_t segconf_bisulfite    : 1; // SAM: v14
            uint8_t segconf_is_paired    : 1; // SAM: v14
            uint8_t segconf_sag_has_AS   : 1; // SAM: v14
            uint8_t segconf_pysam_qual   : 1; // SAM: v14
            uint8_t segconf_cellranger   : 1; // SAM: v14
            uint8_t segconf_SA_HtoS      : 1; // SAM: v14
            uint8_t segconf_is_sorted    : 1; // SAM: v14
            uint8_t segconf_is_collated  : 1; // SAM: v14
            uint8_t segconf_MD_NM_by_un  : 1; // SAM: v14
            uint8_t segconf_predict_meth : 1; // SAM: v14
            uint8_t segconf_deep_qname1  : 1; // SAM: v15
            uint8_t segconf_deep_qname2  : 1; // SAM: v15
            uint8_t segconf_deep_no_qual : 1; // SAM: v15
            uint8_t unused_bits          : 7;
            char unused[255];
        } sam;

        struct { 
            DictId segconf_seq_len_dict_id;   // FASTQ: dict_id of one of the Q?NAME contexts, which is expected to hold the seq_len for this read. 0 if there is no such item. copied from segconf.seq_len_dict_id. added v14.
            char unused[264];
        } fastq;

        struct {
            uint8_t segconf_has_RGQ      : 1; // VCF: copied from segconf.has[FORMAT_RGQ]. added v14.
            uint8_t unused_bits          : 7;
            uint8_t unused[271];
        } vcf;
    };    
} SectionHeaderGenozipHeader, *SectionHeaderGenozipHeaderP;
typedef const SectionHeaderGenozipHeader *ConstSectionHeaderGenozipHeaderP;

// this footer appears AFTER the genozip header data, facilitating reading the genozip header in reverse from the end of the file
typedef struct {
    uint64_t genozip_header_offset;
    uint32_t magic;
} SectionFooterGenozipHeader, *SectionFooterGenozipHeaderP;

// NOTE: part of the file format - new flavors can be added, but their value cannot change (v15)
typedef enum __attribute__ ((__packed__)) { 
    QF_NO_ID=0, 
    QF_ILLUM_7=1, QF_ILLUM_7i=2, QF_ILLUM_7umi=3, QF_ILLUM_7_bc=4, QF_ILLUM_7gs=5, QF_ILLUM_5i=6, QF_ILLUM_5=7, QF_ILLUM_5rng=8,  
    QF_ILLUM_2bc=9, QF_ILLUM_1bc=10, QF_ILLUM_0bc=11, QF_ILLUM_X_0bc=12, QF_ILLUM_X_1bc=13, QF_ILLUM_X_2bc=14, QF_ILLUM_S_0bc=15, QF_ILLUM_S_1bc=16, QF_ILLUM_S_2bc=17, QF_ILLUM_7gsFQ=18, QF_ILLUM_7_2bc=19,
    QF_BGI_varlen=20, QF_BGI_r6=21, QF_BGI_r7=22, QF_BGI_r8=23, QF_BGI_ll7=24, QF_BGI_cl=25, QF_BGI_rgs8=26, QF_BGI_rgs8FQ=27, QF_BGI_Rgs_2bc=28,
    QF_PACBIO_3=30, QF_PACBIO_rng=31, QF_PACBIO_lbl=32, QF_PACBIO_pln=33,
    QF_NANOPORE=40, QF_NANOPORE_rng=41, QF_NANOPORE_ext=42,
    QF_ION_TORR_3=50, QF_ROCHE_454=51, QF_HELICOS=52, 
    QF_SRA_L=60, QF_SRA2=60, QF_SRA=62,
    QF_GENOZIP_OPT=70, QF_INTEGER=71, QF_HEX_CHR=72, QF_BAMSURGEON=73, QF_SEQAN=74, QF_CLC_GW=75, QF_STR_INT=76, QF_CONSENSUS=77,
    QF_ULTIMA_a=80, QF_ULTIMA_b6_bc=81, QF_ULTIMA_a_bc=83, QF_ULTIMA_c=84, QF_ULTIMA_c_bc=85, QF_ULTIMA_b6=86, QF_ULTIMA_d=87, QF_ULTIMA_d_bc=88, 
    QF_ULTIMA_b9=89, QF_ULTIMA_b9_bc=110, QF_ULTIMA_n=111,
    QF_SINGULAR=90, QF_SINGULR_1bc=92, 
    QF_ELEMENT=100, QF_ELEMENT_6=101, QF_ELEMENT_0bc=102, QF_ELEMENT_1bc=103, QF_ELEMENT_2bc=104,
    NUM_FLAVORS
} QnameFlavorId;

typedef enum __attribute__ ((__packed__)) { // these values and their order are part of the file format
                      CNN_NONE, CNN_SEMICOLON, CNN_COLON, CNN_UNDERLINE, CNN_HYPHEN, CNN_HASH, CNN_SPACE, CNN_PIPE, NUM_CNN
} QnameCNN /* 3 bits */;
#define CNN_TO_CHAR { 0,        ';',           ':',       '_',           '-',        '#',      ' ',       '|' }
#define CHAR_TO_CNN { [';']=CNN_SEMICOLON, [':']=CNN_COLON, ['_']=CNN_UNDERLINE, ['-']=CNN_HYPHEN, ['#']=CNN_HASH, [' ']=CNN_SPACE, ['|']=CNN_PIPE }

typedef struct __attribute__ ((__packed__)) { // 3 bytes
    QnameFlavorId id;
    uint8_t has_seq_len : 1;  // qname includes seq_len
    uint8_t unused_bit  : 1;  // 15.0.0 to 15.0.14 : is_mated: qname ends with /1 or /2 (never used in PIZ)
    uint8_t is_mated    : 1;  // qname's final two characters are separator + mate (1 or 2), eg "/1" or "|2" or ".1"
    uint8_t cnn         : 3;  // QnameCNN: terminate before the last character that is this, to canonoize 
    uint8_t unused_bits : 2;  
} QnameFlavorProp;

// The text file header section appears once in the file (or multiple times in case of bound file), and includes the txt file header 
typedef struct {
    SectionHeader;
    uint64_t txt_data_size;            // number of bytes in the original txt file
    uint64_t txt_num_lines;            // number of data (non-header) lines in the original txt file. Concat mode: entire file for first SectionHeaderTxtHeader, and only for that txt if not first
    uint32_t max_lines_per_vb;         // upper bound on how many data lines a VB can have in this file
    Codec    src_codec;                // codec of original txt file (none, bgzf, gz, bz2...)
    uint8_t  codec_info[3];            // codec specific info: for CODEC_BGZF, these are the LSB, 2nd-LSB, 3rd-LSB of the source BGZF-compressed file size
    Digest   digest;                   // digest of original single txt file (except modified or DVCF). v14: only if md5, not alder32 (adler32 digest, starting v14, is stored per VB) (bug in v14: this field is garbage instead of 0 for FASTQ_COMP_FQR2 if adler32)
    Digest   digest_header;            // MD5 or Adler32 of header
#define TXT_FILENAME_LEN 256
    char     txt_filename[TXT_FILENAME_LEN]; // filename of this single component. without path, 0-terminated. always in base form like .vcf or .sam, even if the original is compressed .vcf.gz or .bam
    uint64_t txt_header_size;          // size of header in original txt file (likely different than reconstructed size if dual-coordinates)  (v12)
    QnameFlavorProp flav_prop[NUM_QTYPES]; // SAM/BAM/FASTQ/KRAKEN. properties of QNAME flavor (v15)
} SectionHeaderTxtHeader, *SectionHeaderTxtHeaderP; 

typedef struct {
    SectionHeader;     
    union {
        uint32_t v11_first_line;       // up to v11 - first_line; if 0, this is the terminating section of the components
        uint32_t sam_prim_seq_len;     // SAM PRIM: total number of bases of SEQ in this VB (v14) 
    };    
    union {
        uint32_t v13_top_level_repeats;// v12/13: repeats of TOPLEVEL container in this VB. Up to v12 - called num_lines.
        uint32_t sam_prim_num_sag_alns;// SAM PRIM: number of alns (prim + depn) in SAGs this VB (v14)
    };
    uint32_t recon_size_prim;          // size of vblock as it appears in the default PRIMARY reconstruction
    uint32_t z_data_bytes;             // total bytes of this vblock in the genozip file including all sections and their headers 
    uint32_t longest_line_len;         // length of the longest line in this vblock 
    Digest   digest;                   // stand-alone Adler32 or commulative MD5 up to and including this VB. Up to v13: Adler32 was commulative too.
    union {
        uint32_t v13_num_lines_prim;   // v12/13: number of lines in default reconstruction (DVCF: in PRIMARY coords) (v12)
        uint32_t sam_prim_first_grp_i; // SAM PRIM: the index of first group of this PRIM VB, in z_file->sag_grps (v14)
    };
    
    union {
        uint32_t v13_dvcf_num_lines_luft; // v12/13: DVCF: number of lines in default reconstruction in LUFT coords (v12)
        uint32_t sam_prim_comp_qual_len;  // SAM PRIM: total size of SA Group's QUAL, as compressed in-memory in ZIP (v14)
    };

    union {
        uint32_t dvcf_recon_size_luft; // DVCF: size of vblock as it appears in the default LUFT reconstruction (v12)
        uint32_t sam_prim_qname_len;   // SAM PRIM: total length of QUAL in this VB (v14)
    };

    union {
        uint32_t sam_prim_comp_cigars_len; // SAM PRIM SAG_BY_SA: total size of sag's CIGARs, as compressed in-memory in ZIP (i.e. excluding CIGARs stored in OPTION_SA_CIGAR.dict) (v14)
        uint32_t sam_prim_solo_data_len;   // SAM PRIM SAG_BY_SOLO: size of solo_data
    };
    uint32_t longest_seq_len;          // SAM / FASTQ: largest seq_len in this VB (v15)
} SectionHeaderVbHeader, *SectionHeaderVbHeaderP; 
typedef const SectionHeaderVbHeader *ConstSectionHeaderVbHeaderP;

typedef struct {
    SectionHeader;           
    uint32_t num_snips;                // number of items in dictionary
    DictId   dict_id;           
} SectionHeaderDictionary, *SectionHeaderDictionaryP;

typedef struct {
    SectionHeader;           
    int64_t  nodes_param;              // an extra piece of data transferred to/from Context.counts_extra
    DictId   dict_id;           
} SectionHeaderCounts, *SectionHeaderCountsP;     

// used for SEC_LOCAL and SEC_B250
typedef struct {
    SectionHeader;
    LocalType ltype;           // SEC_LOCAL: goes into ctx.ltype - type of data for the ctx->local buffer
                               // SEC_250:   unused. up to 15.0.17 populated ctx->ltype and copied in PIZ to ctx->ltype if the corresponding local section is absent (not clear why)
    uint8_t param;             // Three options: 1. goes into ctx.local.param. (until v13: if flags.copy_local_param. since v14: always, except if ltype=LT_BITMAP) 
                               //                2. given to comp_uncompress as a codec parameter
                               //                3. starting 9.0.11 for ltype=LT_BITMAP: number of unused bits in top bits word
    B250Size b250_size : 2;    // b250 sections only: size of each b250 element (v14)
    uint8_t unused2    : 6;
    uint8_t unused;
    DictId dict_id;           
} SectionHeaderCtx, *SectionHeaderCtxP;         

// two ways of storing a range:
// uncompacted - we will have one section, SEC_REFERENCE, containing the data, and first/last_pos containing the coordinates of this range
// compacting we will have 2 sections:
// - SEC_REF_IS_SET - containing a bitmap (1 bit per base), and chrom,first/last_pos containing the coordinates of this range
// - SEC_REFERENCE - containing the compacted range (i.e. excluding bases that have a "0" in the bitmap), 
//   with chrom,first_pos same as SEC_REF_IS_SET and last_pos set so that (last_pos-first_pos+1) = number of '1' bits in SEC_REF_IS_SET
// SEC_REFERENCE (in both cases) contains 2 bits per base, and SEC_REF_IS_SET contains 1 bit per location.
typedef struct {
    SectionHeader;
    PosType64 pos;             // first pos within chrom (1-based) of this range         
    PosType64 gpos;            // first pos within genome (0-based) of this range
    uint32_t num_bases;        // number of bases (nucleotides) in this range
    uint32_t chrom_word_index; // index in contexts[CHROM].word_list of the chrom of this reference range    
} SectionHeaderReference, *SectionHeaderReferenceP;

typedef struct {
    SectionHeader;
    uint8_t num_layers;        // total number of layers
    uint8_t layer_i;           // layer number of this section (0=base layer, with the most bits)
    uint8_t layer_bits;        // number of bits in layer
    uint8_t ffu;
    uint32_t start_in_layer;   // start index within layer
} SectionHeaderRefHash, *SectionHeaderRefHashP;

// SEC_RECON_PLAN, contains ar array of ReconPlanItem
typedef struct SectionHeaderReconPlan {
    SectionHeader;
    VBIType conc_writing_vbs;  // max number of concurrent VBs in possesion of the writer thread needed to execute this plan    
    uint32_t vblock_mb;        // size of vblock in MB
} SectionHeaderReconPlan, *SectionHeaderReconPlanP;

// plan flavors (3 bit)
typedef enum { PLAN_RANGE=0, PLAN_FULL_VB=2, PLAN_INTERLEAVE=3/*PIZ-only*/, PLAN_TXTHEADER=4/*PIZ-only*/,
               PLAN_REMOVE_ME=5/*PIZ-only*/, PLAN_DOWNSAMPLE=6/*PIZ-only*/, PLAN_END_OF_VB=7 } PlanFlavor;
#define PLAN_FLAVOR_NAMES { "RANGE", "invalid", "FULL_VB", "INTERLEAVE", "TXTHEADER", "REMOVE_ME", "DOWNSAMPLE", "END_OF_VB" }
typedef struct {
    VBIType vb_i;               
    union {
        uint32_t start_line;   // used for RANGE
        uint32_t vb2_i;        // used for INTERLEAVE
        uint32_t comp_i;       // used for TXTHEADER
        uint32_t word2;        // generic access to the value
    }; 
    uint32_t num_lines : 29;   // used by RANGE, FULL_VB (writer only, not file), INTERLEAVE, DOWNSAMPLE. note: in v12/13 END_OF_VB, this field was all 1s.
    PlanFlavor flavor  : 3;  
} ReconPlanItem, *ReconPlanItemP;

// the data of SEC_SECTION_LIST is an array of the following type, as is the z_file->section_list_buf
typedef struct __attribute__ ((__packed__)) SectionEntFileFormat { // 19 bytes (since v15)
    uint32_t offset_delta;     // delta vs previous offset (always positive)
    int32_t vblock_i_delta;    // delta vs previous section's vblock_i (may be negative) 
    CompIType comp_i_plus_1;   // stores 0 if same as prev section, COMP_NONE if COMP_NONE, and otherwise comp_i+1.
    SectionType st;            // 1 byte
    union {                    // Section-Type-specific field
        DictId dict_id;        // DICT, LOCAL, B250 or COUNT sections (first occurance of this dict_id)
        struct {               // DICT, LOCAL, B250 or COUNT sections (2nd+ occurance of this dict_id)
            char is_dict_id;   // if zero, dict_id is copied from earlier section
            char unused3[3];
            uint32_t dict_sec_i; // section from which to copy the dict_id
        };
        struct { 
            uint32_t num_lines;      // VB_HEADER sections - number of lines in this VB. 
            uint32_t num_non_deep_lines; // VB_HEADER Deep SAM, component MAIN and PRIM: number of lines in this VB that are NOT deepable, i.e. SUPPLEMENTARY or SECONDARY or monochar or with SEQ.len=0  
        };
        uint64_t st_specific;  // generic access to the value
    };
    SectionFlags flags;        // same flags as in section header (since v12, previously "unused")
} SectionEntFileFormat;

typedef struct SectionEntFileFormatV14 { // 24 bytes (used up to v14)
    uint64_t offset;           // offset of this section in the file
    union {                    // Section-Type-specific field
        DictId dict_id;        // DICT, LOCAL, B250 or COUNT sections
        struct { 
            uint32_t num_lines;// VB_HEADER sections - number of lines in this VB. 
            uint32_t unused;
        };
        uint64_t st_specific;  // generic access to the value
    };
    VBIType vblock_i;          // 1-based
    SectionType st;            // 1 byte
    SectionFlags flags;        // same flags as in section header (since v12, previously "unused")
    uint8_t comp_i  : 2;       // used for component-related sections, 0-based (since v14, previously "unused"), (0 if COMP_NONE)
    uint8_t unused1 : 6;
    uint8_t unused2;         
} SectionEntFileFormatV14;     // up to v14 

// the data of SEC_RANDOM_ACCESS is an array of the following type, as is the z_file->ra_buf and vb->ra_buf
// we maintain one RA entry per vb per every chrom in the the VB
typedef struct RAEntry {
    VBIType vblock_i;          // the vb_i in which this range appears
    WordIndex chrom_index;     // before merge: node index into chrom context nodes, after merge - word index in CHROM dictionary
    PosType64 min_pos, max_pos;// POS field value of smallest and largest POS value of this chrom in this VB (regardless of whether the VB is sorted)
} RAEntry; 

// the data of SEC_REF_IUPACS (added v12)
typedef struct Iupac {
    PosType64 gpos;
    char iupac;
} Iupac;

typedef union {
    SectionHeader common;
    SectionHeaderGenozipHeader genozip_header;
    SectionHeaderTxtHeader txt_header;
    SectionHeaderVbHeader vb_header;
    SectionHeaderDictionary dict;
    SectionHeaderCounts counts;
    SectionHeaderCtx ctx;
    SectionHeaderReference reference;
    SectionHeaderRefHash ref_hash;
    SectionHeaderReconPlan recon_plan;
    char padding[sizeof(SectionHeaderGenozipHeader) + 15]; // SectionHeaderGenozipHeader is the largest, 15=crypt_max_padding_len()
} SectionHeaderUnion;

typedef union {
    SectionHeaderP common;
    SectionHeaderGenozipHeaderP genozip_header;
    SectionHeaderTxtHeaderP txt_header;
    SectionHeaderVbHeaderP vb_header;
    SectionHeaderDictionaryP dict;
    SectionHeaderCountsP counts;
    SectionHeaderCtxP ctx;
    SectionHeaderReferenceP reference;
    SectionHeaderRefHashP ref_hash;
    SectionHeaderReconPlanP recon_plan;
} SectionHeaderUnionP __attribute__((__transparent_union__));

#pragma pack()

// in-memory section 
typedef const struct SectionEnt {
    uint64_t offset;            // offset of this section in the file
    union {                     // Section-Type-specific field
        DictId dict_id;         // DICT, LOCAL, B250 or COUNT sections
        struct {
            uint32_t num_lines; // VB_HEADER sections - number of lines in this VB. 
            uint32_t num_deep_lines; // VB_HEADER SAM: number of lines in this VB that are not SUPPLEMENTARY or SECONDARY (v15) 
        };
        struct { 
            char is_dict_id;    // if zero, dict_id is copied from dict_sec_i (dict_id.id[0] cannot be zero = a DTYPE_FIELD dict_id cannot start with '@')
            char unused3[3];
            uint32_t dict_sec_i;// section from which to copy the dict_id
        };
        uint64_t st_specific;   // generic access to the value
    };
    VBIType vblock_i;           // 1-based
    uint32_t size;        
    CompIType comp_i;           // 0-based. 255 if COMP_NONE.
    SectionType st;             // 1 byte
    SectionFlags flags;         // same flags as in section header, since v12 (before was "unused")
} SectionEnt;

// alias types in SEC_DICT_ID_ALIASES (since v15)
typedef enum __attribute__ ((__packed__)) { ALIAS_NONE, ALIAS_CTX, ALIAS_DICT } AliasType;
#define ALIAS_TYPE_NAMES                  { "NONE",     "CTX",      "DICT"    }

// ---------
// ZIP stuff
// ---------

extern void sections_add_to_list (VBlockP vb, ConstSectionHeaderP header);
extern void sections_remove_from_list (VBlockP vb, uint64_t offset, uint64_t len);
extern void sections_list_concat (VBlockP vb);

// ---------
// PIZ stuff
// ---------

extern Section section_next (Section sec);

extern Section sections_first_sec (SectionType st, FailType soft_fail);
extern Section sections_last_sec (SectionType st, FailType soft_fail);
extern bool sections_next_sec3 (Section *sl_ent, SectionType st1, SectionType st2, SectionType st3);
#define sections_next_sec(sl_ent,st)       sections_next_sec3((sl_ent),(st),SEC_NONE,SEC_NONE)
#define sections_next_sec2(sl_ent,st1,st2) sections_next_sec3((sl_ent),(st1),(st2),SEC_NONE)

extern bool sections_prev_sec2 (Section *sl_ent, SectionType st1, SectionType st2);
#define sections_prev_sec(sl_ent,st) sections_prev_sec2((sl_ent),(st),SEC_NONE)

extern uint32_t sections_count_sections_until (SectionType st, Section first_sec, SectionType until_encountering);
static inline uint32_t sections_count_sections (SectionType st) { return sections_count_sections_until (st, 0, SEC_NONE); }

extern void sections_create_index (void);
extern Section sections_vb_header (VBIType vb_i);
extern VBIType sections_get_num_vbs (CompIType comp_i);
extern VBIType sections_get_num_vbs_(CompIType first_comp_i, CompIType last_comp_i);
extern VBIType sections_get_first_vb_i (CompIType comp_i);
extern Section sections_get_comp_txt_header_sec (CompIType comp_i);
extern Section sections_get_comp_recon_plan_sec (CompIType comp_i, bool is_luft_plan);
extern Section sections_get_comp_bgzf_sec (CompIType comp_i);
extern Section sections_get_next_vb_header_sec (CompIType comp_i, Section *vb_sec);
extern Section sections_one_before (Section sec);
extern bool is_there_any_section_with_dict_id (DictId dict_id);

// API for writer to create modified section list
extern void sections_new_list_add_vb (BufferP new_list, VBIType vb_i);
extern void sections_new_list_add_txt_header (BufferP new_list, CompIType comp_i);
extern void sections_new_list_add_bgzf (BufferP new_list);
extern void sections_new_list_add_global_sections (BufferP new_list);
extern void sections_commit_new_list (BufferP new_list);

extern void sections_list_memory_to_file_format (bool in_place);
extern void sections_list_file_to_memory_format (SectionHeaderGenozipHeaderP genozip_header);
extern bool sections_is_paired (void);

#define sections_has_dict_id(st) ((st) == SEC_B250 || (st) == SEC_LOCAL || (st) == SEC_DICT || (st) == SEC_COUNTS)
extern SectionType sections_st_by_name (char *name);
extern uint32_t st_header_size (SectionType sec_type);

// display functions
#define sections_read_prefix(P_prefix) ((P_prefix) ? 'P' : flag_loading_auxiliary ? 'L' : 'R')
extern void sections_show_header (ConstSectionHeaderP header, VBlockP vb /* optional if output to buffer */, uint64_t offset, char rw);
extern void noreturn genocat_show_headers (rom z_filename);
extern void sections_show_gheader (ConstSectionHeaderGenozipHeaderP header);
extern void sections_show_section_list (DataType dt);
extern rom st_name (SectionType sec_type);
extern rom lt_name (LocalType lt);
extern rom store_type_name (StoreType store);
extern DictId sections_get_dict_id (ConstSectionHeaderP header);

extern StrText vb_name (VBlockP vb);
#define VB_NAME vb_name(VB).s

extern StrText line_name (VBlockP vb);
#define LN_NAME line_name(VB).s

extern StrText comp_name_(CompIType comp_i);
#define comp_name(comp_i) comp_name_(comp_i).s

#define IS_DICTED_SEC(st) ((st)==SEC_DICT || (st)==SEC_B250 || (st)==SEC_LOCAL || (st)==SEC_COUNTS)
#define IS_VB_SEC(st)     ((st)==SEC_VB_HEADER || (st)==SEC_B250 || (st)==SEC_LOCAL)
#define IS_COMP_SEC(st)   (IS_VB_SEC(st) || (st)==SEC_TXT_HEADER || (st)==SEC_BGZF || (st)==SEC_RECON_PLAN)
#define IS_FRAG_SEC(st)   ((st)==SEC_DICT || (st)==SEC_TXT_HEADER || (st)==SEC_RECON_PLAN || (st)==SEC_REFERENCE || (st)==SEC_REF_IS_SET || (st)==SEC_REF_HASH) // global sections fragmented with a dispatcher, and hence use vb_i 

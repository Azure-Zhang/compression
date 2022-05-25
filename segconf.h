// ------------------------------------------------------------------
//   segconf.c
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"

// Documented range for users
#define MIN_VBLOCK_MEMORY  1    // in MB 
#define MAX_VBLOCK_MEMORY  2048 

// Range in developer options eg --vblock 10000B 
#define ABSOLUTE_MIN_VBLOCK_MEMORY ((uint64_t)1000) // in Bytes
#define ABSOLUTE_MAX_VBLOCK_MEMORY ((uint64_t)MAX_VBLOCK_MEMORY<<20)

typedef enum __attribute__ ((__packed__)) { TECH_UNKNOWN, TECH_ILLUM_7, TECH_ILLUM_5, TECH_PACBIO, TECH_ONP, TECH_454, TECH_BGI, TECH_IONTORR, TECH_HELICOS } SeqTech;

typedef enum __attribute__ ((__packed__)) { SQT_UNKNOWN, SQT_NUKE, SQT_AMINO, SQT_NUKE_OR_AMINO } SeqType;

typedef enum __attribute__ ((__packed__)) { PL_mux_by_DP_TEST, PL_mux_by_DP_NO, PL_mux_by_DP_YES } PLMuxByDP;

typedef enum __attribute__ ((__packed__)) { PS_NONE, PS_POS, PS_POS_REF_ALT, PS_UNKNOWN } PSType;

typedef enum __attribute__ ((__packed__)) { ms_NONE, ms_BIOBAMBAM, ms_MINIMAP2 } msType; // type of SAM ms:i field 

typedef enum __attribute__ ((__packed__)) { DP_DEFAULT, by_AD, by_SDP, by_INFO_DP } FormatDPMethod;

typedef enum __attribute__ ((__packed__)) { L3_UNKNOWN, L3_EMPTY, L3_COPY_DESC, L3_QF, NUM_L3s } FastqLine3Type;

typedef enum __attribute__ ((__packed__)) { XG_S_UNKNOWN, XG_WITHOUT_S, XG_WITH_S } XgIncSType;

// seg configuration set prior to starting to seg a file during segconfig_calculate or txtheader_zip_read_and_compress
typedef struct {

    // Seg parameters - general
    uint64_t vb_size;           // ZIP/PIZ: compression VBlock size in bytes (PIZ: passed in SectionHeaderGenozipHeader.vb_size)
    bool running;               // currently in segconf_calculate()
    bool has[MAX_DICTS];        // for select did_i's, states whether this field was encountered during segconf.running
    uint32_t line_len;          // approx line len
    float b250_per_line[MAX_DICTS]; // b250.len / num_lines
    #define AT_LEAST(did_i) ((uint64_t)(10.0 + (segconf.b250_per_line[did_i] * (float)(vb->lines.len32))))

    // read characteristics (SAM/BAM, KRAKEN and FASTQ)
    QnameFlavor qname_flavor, qname_flavor2;  
    SeqTech tech;

    // SAM/BAM and FASTQ
    uint32_t longest_seq_len;   // length of the longest seq_len in the segconf data 

    // SAM/BAM stuff
    bool sam_is_unmapped;       // all POS fields in the segconf block were 0
    bool sam_bowtie2;           
    bool has_bsseeker2;
    bool NM_is_integer;         // true if NM is integer, false if it binary
    bool has_TLEN_non_zero;
    bool has_DP_before_PL;
    enum { XA_NONE, XA_BWA, XA_IONTORRENT, XA_UNKNOWN } XA_type; // IonTorret and BWA have different XA:Z
    bool sam_is_collated;       // Every QNAME appears in two or more consecutive lines
    bool sam_is_sorted;         // every two consecutive lines that have the same RNAME, have non-decreasing POS
    bool sam_is_paired;         // file has a least one read that is marked as "last" in FLAG
    bool sam_multi_RG;          // evidence that file has more than one type of RG
    bool has_MD_or_NM;          // ZIP/ZIP: call sam_analyze_copied_SEQ for SEQ copied from prim unless no cigar, no seq or (PIZ) explicitly told not to.
    bool NM_after_MD;           // in all segconf lines that had both NM and MD, NM appeared after MD
    uint8_t MAPQ_value;         // used during segconf.running to calculate sam_mapq_has_single_value
    bool MAPQ_has_single_value; // all non-0 MAPQ have the same value
    msType sam_ms_type;         // ZIP/PIZ: type of ms:i 
    XgIncSType sam_XG_inc_S;    // Does XG include soft_clip[0]
    bool is_long_reads;
    uint32_t sam_cigar_len;     // approx average CIGAR len rounded up (during running==true - total len)
    uint32_t sam_seq_len;       // ZIP/PIZ: approx average (SEQ.len+hard-clips) rounded to the nearest (during running==true - total len)

    // SAM/BAM and FASTQ
    bool nontrivial_qual;       // true if we know that not all QUAL values are the same (as they are in newer PacBio files)

    // VCF stuff
    bool vcf_is_varscan;        // this VCF file was produced by VarScan
    uint64_t count_dosage[2];   // used to calculate pc_has_dosage
    float pc_has_dosage;        // % of the samples x lines that have a valid (0-2) dosage value [0.0,1.0]
    PSType ps_pid_type[2];      // [0]=PS [1]=PID
    bool use_null_DP_method;    // A method for predicting GT=./. by DP=.
    bool INFO_DP_by_FORMAT_DP;  // INFO_DP = Delta vs SUM(FORMAT_DP) 
    FormatDPMethod FORMAT_DP_method;
    PLMuxByDP PL_mux_by_DP;
    Mutex PL_mux_by_DP_mutex;
    uint64_t count_GQ_by_PL, count_GQ_by_GP; // used tp calculate GQ_by_PL, GQ_by_GP
    bool GQ_by_PL, GQ_by_GP;
    
    // FASTQ
    FastqLine3Type line3;       // format of line3
    QnameFlavor line3_flavor;   // in case of L3_QF 

    // FASTA stuff
    bool fasta_has_contigs;     // the sequences in this FASTA represent contigs (as opposed to reads) - in which case we have a FASTA_CONTIG dictionary and RANDOM_ACCESS
    SeqType seq_type;           // nucleotide or protein
    unsigned seq_type_counter;  // used for calculating seq_type 

    // Chain stuff
    bool chain_mismatches_ref;  // Some contigs mismatch the reference files, so this chain file cannot be used with --chain
} SegConf;

extern SegConf segconf; // ZIP: set based on segging a sample of a few first lines of the file
                        // PIZ: select fields are transferred through SectionHeaderGenozipHeader

extern void segconf_initialize (void);
extern void segconf_calculate (void);
extern void segconf_update_qual (STRp (qual));
extern bool segconf_is_long_reads(void);
extern void segconf_mark_as_used (VBlockP vb, unsigned num_ctxs, ...);

static inline void segconf_set_has (DidIType did_i)
{
    if (segconf.running) segconf.has[did_i] = true;
}

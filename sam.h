// ------------------------------------------------------------------
//   sam.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef SAM_INCLUDED
#define SAM_INCLUDED

#include "genozip.h"
#include "sections.h"

#pragma GENDICT_PREFIX SAM

// Fields
#pragma GENDICT SAM_RNAME=DTYPE_FIELD=RNAME // RNAME must be first
#pragma GENDICT SAM_QNAME=DTYPE_FIELD=QNAME // MAX_QNAME_ITEMS
#pragma GENDICT SAM_Q0NAME=DTYPE_1=Q0NAME // must have a did_i directly after container's 
#pragma GENDICT SAM_Q1NAME=DTYPE_1=Q1NAME 
#pragma GENDICT SAM_Q2NAME=DTYPE_1=Q2NAME
#pragma GENDICT SAM_Q3NAME=DTYPE_1=Q3NAME
#pragma GENDICT SAM_Q4NAME=DTYPE_1=Q4NAME
#pragma GENDICT SAM_Q5NAME=DTYPE_1=Q5NAME
#pragma GENDICT SAM_Q6NAME=DTYPE_1=Q6NAME 
#pragma GENDICT SAM_Q7NAME=DTYPE_1=Q7NAME 
#pragma GENDICT SAM_QmatNAME=DTYPE_1=QmatNAME // QmatNAME reserved for mate number (always the last dict_id in the container)
#pragma GENDICT SAM_FLAG=DTYPE_FIELD=FLAG
#pragma GENDICT SAM_POS=DTYPE_FIELD=POS
#pragma GENDICT SAM_MAPQ=DTYPE_FIELD=MAPQ
#pragma GENDICT SAM_CIGAR=DTYPE_FIELD=CIGAR
#pragma GENDICT SAM_RNEXT=DTYPE_FIELD=RNEXT
#pragma GENDICT SAM_PNEXT=DTYPE_FIELD=PNEXT
#pragma GENDICT SAM_TLEN=DTYPE_FIELD=TLEN
#pragma GENDICT SAM_OPTIONAL=DTYPE_FIELD=OPTIONAL
#pragma GENDICT SAM_SQBITMAP=DTYPE_FIELD=SQBITMAP
#pragma GENDICT SAM_NONREF=DTYPE_FIELD=NONREF    // these 4 fields must be in this order, right after SAM_SQBITMAP
#pragma GENDICT SAM_NONREF_X=DTYPE_FIELD=NONREF_X
#pragma GENDICT SAM_GPOS=DTYPE_FIELD=GPOS
#pragma GENDICT SAM_STRAND=DTYPE_FIELD=STRAND
#pragma GENDICT SAM_QUAL=DTYPE_FIELD=QUAL 
#pragma GENDICT SAM_DOMQRUNS=DTYPE_FIELD=DOMQRUNS // must be right after SAM_QUAL
#pragma GENDICT SAM_EOL=DTYPE_FIELD=EOL
#pragma GENDICT SAM_BAM_BIN=DTYPE_FIELD=BAM_BIN
#pragma GENDICT SAM_TOPLEVEL=DTYPE_FIELD=TOPLEVEL // must be called TOPLEVEL
#pragma GENDICT SAM_TOP2BAM=DTYPE_FIELD=TOP2BAM
#pragma GENDICT SAM_TOP2FQ=DTYPE_FIELD=TOP2FQ
#pragma GENDICT SAM_TOP2FQEX=DTYPE_FIELD=TOP2FQEX
#pragma GENDICT SAM_TAXID=DTYPE_FIELD=TAXID
#pragma GENDICT SAM_BUDDY=DTYPE_FIELD=BUDDY    // note: this MUST be the same dict_id ("BUDDY") for all data_types using buddy (they will have different did_i though), expected by reconstruct_from_buddy 
#pragma GENDICT SAM_MC_Z=DTYPE_FIELD=MC_Z  // used for consuming OPTION_MC_Z in case of translation to FASTQ

#pragma GENDICT OPTION_AM_i=DTYPE_2=AM:i // The smallest template-independent mapping quality in the template
#pragma GENDICT OPTION_AS_i=DTYPE_2=AS:i // SAM: Alignment score generated by aligner ; STAR: the local alignment score (paired for paired-end reads).
#pragma GENDICT OPTION_CC_Z=DTYPE_2=CC:Z // Reference name of the next hit
#pragma GENDICT OPTION_CM_i=DTYPE_2=CM:i // Edit distance between the color sequence and the color reference (see also NM)
#pragma GENDICT OPTION_E2_Z=DTYPE_2=E2:Z // The 2nd most likely base calls
#pragma GENDICT OPTION_2NONREF=DTYPE_2=N2ONREF // these 4 fields must be in this order, right after OPTION_E2_Z
#pragma GENDICT OPTION_N2ONREFX=DTYPE_2=n2ONREFX
#pragma GENDICT OPTION_2GPOS=DTYPE_FIELD=G2POS
#pragma GENDICT OPTION_S2TRAND=DTYPE_2=S2TRAND
#pragma GENDICT OPTION_FI_i=DTYPE_2=FI:i // The index of segment in the template
#pragma GENDICT OPTION_H0_i=DTYPE_2=H0:i // Number of perfect hits
#pragma GENDICT OPTION_H1_i=DTYPE_2=H1:i // Number of 1-difference hits (see also NM)
#pragma GENDICT OPTION_H2_i=DTYPE_2=H2:i // Number of 2-difference hits
#pragma GENDICT OPTION_LB_Z=DTYPE_2=LB:Z // Library
#pragma GENDICT OPTION_MC_Z=DTYPE_2=MC:Z // CIGAR string for mate/next segment
#pragma GENDICT OPTION_MD_Z=DTYPE_2=MD:Z // String encoding mismatched and deleted reference bases
#pragma GENDICT OPTION_MQ_i=DTYPE_2=MQ:i // Mapping quality of the mate/next segment
#pragma GENDICT OPTION_NH_i=DTYPE_2=NH:i // Number of reported alignments that contain the query in the current record
#pragma GENDICT OPTION_NM_i=DTYPE_2=NM:i // Edit distance to the reference

#pragma GENDICT OPTION_OA_Z=DTYPE_2=OA:Z // Original alignment
#pragma GENDICT OPTION_OA_RNAME=DTYPE_2=O0A_RNAME  
#pragma GENDICT OPTION_OA_STRAND=DTYPE_2=O1A_STRAND 
#pragma GENDICT OPTION_OA_POS=DTYPE_2=O2A_POS 
#pragma GENDICT OPTION_OA_CIGAR=DTYPE_2=O3A_CIGAR 
#pragma GENDICT OPTION_OA_NM=DTYPE_2=O4A_NM 
#pragma GENDICT OPTION_OA_MAPQ=DTYPE_2=O5A_MAPQ 

#pragma GENDICT OPTION_OC_Z=DTYPE_2=OC:Z // Original CIGAR
#pragma GENDICT OPTION_PG_Z=DTYPE_2=PG:Z // Program
#pragma GENDICT OPTION_PQ_i=DTYPE_2=PQ:i // Phred likelihood of the template
#pragma GENDICT OPTION_PU_Z=DTYPE_2=PU:Z // Platform unit
#pragma GENDICT OPTION_RG_Z=DTYPE_2=RG:Z // Read group

#pragma GENDICT OPTION_SA_Z=DTYPE_2=SA:Z // Other canonical alignments in a chimeric alignment
#pragma GENDICT OPTION_SA_RNAME=DTYPE_2=S0A_RNAME 
#pragma GENDICT OPTION_SA_STRAND=DTYPE_2=S1A_STRAND 
#pragma GENDICT OPTION_SA_POS=DTYPE_2=S2A_POS 
#pragma GENDICT OPTION_SA_CIGAR=DTYPE_2=S3A_CIGAR 
#pragma GENDICT OPTION_SA_NM=DTYPE_2=S4A_NM 
#pragma GENDICT OPTION_SA_MAPQ=DTYPE_2=S5A_MAPQ 

#pragma GENDICT OPTION_SM_i=DTYPE_2=SM:i   // Template-independent mapping quality
#pragma GENDICT OPTION_TC_i=DTYPE_2=TC:i   // The number of segments in the template
#pragma GENDICT OPTION_UQ_i=DTYPE_2=UQ:i   // Phred likelihood of the segment, conditional on the mapping being correct
#pragma GENDICT OPTION_U2_Z=DTYPE_2=U2:Z   // Phred probability of the 2nd call being wrong conditional on the best being wrong
#pragma GENDICT OPTION_D2OMQRUN=DTYPE_2=D2OMQRUN // must be right after OPTION_U2_Z
                
// bwa tags see here: http://bio-bwa.sourceforge.net/bwa.shtml : "SAM ALIGNMENT FORMAT"
#pragma GENDICT OPTION_X0_i=DTYPE_2=X0:i   // Number of best hits
#pragma GENDICT OPTION_X1_i=DTYPE_2=X1:i   // Number of suboptimal hits found by BWA
#pragma GENDICT OPTION_XC_i=DTYPE_2=XC:i   // ???
#pragma GENDICT OPTION_XN_i=DTYPE_2=XN:i   // Number of ambiguous bases in the referenece
#pragma GENDICT OPTION_XM_i=DTYPE_2=XM:i   // Number of mismatches in the alignment
#pragma GENDICT OPTION_XO_i=DTYPE_2=XO:i   // Number of gap opens
#pragma GENDICT OPTION_XG_i=DTYPE_2=XG:i   // Number of gap extentions
#pragma GENDICT OPTION_XT_A=DTYPE_2=XT:A   // Unique/Repeat/N/Mate-sw
#pragma GENDICT OPTION_XS_i=DTYPE_2=XS:i   // Suboptimal alignment score
#pragma GENDICT OPTION_XE_i=DTYPE_2=XE:i   // Number of supporting seeds
// type not known #pragma GENDICT OPTION_XF_?=DTYPE_2=XF:?  // Support from forward/reverse alignment  

#pragma GENDICT OPTION_XA_Z=DTYPE_2=XA:Z  // (OVERLAP WITH Ion Torrent XA:Z) Alternative hits; format: (chr,pos,CIGAR,NM;)*. 
#pragma GENDICT OPTION_XA_LOOKBACK=DTYPE_2=X^A_LOOKBACK
#pragma GENDICT OPTION_XA_RNAME=DTYPE_2=X0A_RNAME
#pragma GENDICT OPTION_XA_STRAND=DTYPE_2=X1A_STRAND 
#pragma GENDICT OPTION_XA_POS=DTYPE_2=X2A_POS 
#pragma GENDICT OPTION_XA_CIGAR=DTYPE_2=X3A_CIGAR 
#pragma GENDICT OPTION_XA_NM=DTYPE_2=X4A_NM 
#pragma GENDICT OPTION_XA_STRAND_POS=DTYPE_2=X5A_STRAND_POS 

// Ion Torrent Base Caller tags (source: "Torrent Suite Software 5.12 Help": http://192.167.218.6/ion-docs/GUID-965C5ED4-20C8-45D5-AF07-8B0008AF74AD.html)
#pragma GENDICT OPTION_ZA_i=DTYPE_2=ZA:i   // Number of library insert bases, where the library insert is defined as the sequence after the key and barcode adapter, and before the 3' adapter. (Only present if a 3' adapter was found.)
#pragma GENDICT OPTION_ZB_i=DTYPE_2=ZB:i   // Number of overlapping adapter bases. (Only present if a 3' adapter was found.)
#pragma GENDICT OPTION_ZC_B_i=DTYPE_2=ZC:B // B:i A vector of the following four values (only present if a 3' adapter was found): 1. The zero-based flow during which the first base of the adapter was incorporated (same as ZG). 2. The zero-based flow corresponding to the last insert base 3. Length of the last insert homopolymer 4. Zero-based index of adapter type found.
#pragma GENDICT OPTION_ZF_i=DTYPE_2=ZF:i   // The zero-indexed flow position corresponding to the first template base after 5' trimmed region
#pragma GENDICT OPTION_ZG_i=DTYPE_2=ZG:i   // The zero-based flow during which the first base of the adapter was incorporated. (Present only if a 3' adapter was found.)
#pragma GENDICT OPTION_ZM_B=DTYPE_2=ZM:B   // B:s Normalized signals, which include phasing effects. Stored as floor(256*value).
#pragma GENDICT OPTION_ZP_B=DTYPE_2=ZP:B   // B:f The estimated phase parameters for the read. The values are stored in the order CF (carry forward), IE (incomplete extension), and DR (droop).
#pragma GENDICT OPTION_ZT_Z=DTYPE_2=ZT:Z   // The trimmed 5’ unique molecular tag sequence. Written only if a tag was trimmed.
#pragma GENDICT OPTION_YT_Z=DTYPE_2=YT:Z   // The trimmed 3’ unique molecular tag sequence. Written only if a tag was trimmed.
#pragma GENDICT OPTION_ZE_Z=DTYPE_2=ZE:Z   // The 5’ trimmed sequence removed by the extra-trim-left command. Written only if a sequence was trimmed.
#pragma GENDICT OPTION_YE_Z=DTYPE_2=YE:Z   // The 3’ trimmed sequence removed by the extra-trim-right command. Written only if a sequence was trimmed.
#pragma GENDICT OPTION_ZK_Z=DTYPE_2=ZK:Z   // The trimmed 3' portion of read group specific identifiers that can vary within a read group. Written only if a tag was trimmed.
#pragma GENDICT OPTION_YK_Z=DTYPE_2=YK:Z   // The trimmed 3' portion of read group specific identifiers that can vary within a read group. Written only if a sequence was trimmed.

// Ion Torrent TMAP tags
// #pragma GENDICT OPTION_XA_Z=DTYPE_2=XA:Z // (OVERLAP WITH BWA XA:Z) The algorithm that produced this mapping and from what stage. The format is the algorithm name and the zero-based stage (separated by a dash).
// #pragma GENDICT OPTION_XM_i=DTYPE_2=XM:i // (OVERLAP WITH BWA XM:i) The target length, that is, the number of reference bases spanned by the alignment.
// #pragma GENDICT OPTION_XS_i=DTYPE_2=XS:i // (OVERLAP WITH BWA XS:i) The alignment score of the next-best suboptimal mapping.

// Bowtie tags. Source: http://bowtie-bio.sourceforge.net/manual.shtml#sam-bowtie-output
// #pragma GENDICT OPTION_XM_i=DTYPE_2=XM:i // (OVERLAP WITH BWA XM:i) Number of alignments (up to cap+1)
#pragma GENDICT OPTION_XA_i=DTYPE_2=XA:i    // Aligned read belongs to stratum <N>

// STAR aligner tags. Source: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
#pragma GENDICT OPTION_nM_i=DTYPE_2=nM:i    // the number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches in each mate.
#pragma GENDICT OPTION_jM_B=DTYPE_2=jM:B    // jM:B:c,M1,M2,... intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
#pragma GENDICT OPTION_jI_B=DTYPE_2=jI:B    // jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1-based).

// biobambam tags
#pragma GENDICT OPTION_mc_i=DTYPE_2=mc:i    // MateCoordinate (biobambam)
#pragma GENDICT OPTION_ms_i=DTYPE_2=ms:i    // MateBaseScore  (biobambam)

// Novocraft tags: http://www.novocraft.com/documentation/novosort-2/
#pragma GENDICT OPTION_Z5_i=DTYPE_2=Z5:i    // During the input phase for unsorted or name sorted alignments, Novosort calculates the signature of each read and adds this as a SAM tag (Z5:i:) to other segments of the template. Later, during the processing of the sorted alignments, we can determine the signature of a read and, from the Z5 tag, the signature of it’s pair. Reads are then grouped according to the two signatures, strand & library and duplicates detected within a group.
#pragma GENDICT OPTION_Zq_i=DTYPE_2=Zq:i

// added by GATK's BQSR (Base Quality Score Recalibration)
#pragma GENDICT OPTION_BD_Z=DTYPE_2=BD:Z    // Deletion base quality  (not used in newer versions of GATK)
#pragma GENDICT OPTION_BI_Z=DTYPE_2=BI:Z    // Insertion base quality (not used in newer versions of GATK)
#pragma GENDICT OPTION_BD_BI=DTYPE_2=BD_BI

#pragma GENDICT OPTION_TX_i=DTYPE_2=tx:i    // Genozip tag for taxonomy ID

// backward compatability for decompressing files compressed with older versions that had aliases to these destinations (used by ctx_initialize_predefined_ctxs)
#pragma GENDICT OPTION_CIGAR=DTYPE_2=@CIGAR // For files compressed with 12.0.37 or older which had aliases MC:Z, OC:Z -> @CIGAR
#pragma GENDICT SAM_E2_Z=DTYPE_FIELD=E2:Z   // This used to be the destination alias from OPTION_E2_Z 
#pragma GENDICT SAM_U2_Z=DTYPE_FIELD=U2:Z   // This used to be the destination alias from OPTION_U2_Z 

// ZIP Stuff
COMPRESSOR_CALLBACK(sam_zip_qual)
COMPRESSOR_CALLBACK(sam_zip_U2)
COMPRESSOR_CALLBACK(sam_zip_BD_BI)
extern void sam_zip_initialize (void);
extern bool sam_zip_is_unaligned_line (const char *line, int len);
extern bool sam_zip_dts_flag (void);

// HEADER stuff
extern bool sam_header_inspect (VBlockP txt_header_vb, BufferP txt_header, struct FlagsTxtHeader txt_header_flags);
extern void sam_header_finalize (void);
extern ContigPkgP sam_hdr_contigs;

// SEG Stuff
extern void sam_zip_initialize (void);
extern void sam_seg_initialize (VBlockP vb);
extern void sam_seg_finalize (VBlockP vb);
extern bool sam_seg_is_small (ConstVBlockP vb, DictId dict_id);
extern const char *sam_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);
extern uint32_t sam_seg_seq_len_from_cigar (STRp(cigar));
extern uint32_t sam_seg_get_seq_len_by_MD_field (STRp(md_str));

// PIZ Stuff
extern bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id);
extern bool sam_piz_read_one_vb (VBlockP vb, Section sl);
extern void sam_reconstruct_vb ();
extern CONTAINER_FILTER_FUNC (sam_piz_filter);
extern void sam_reconstruct_seq (VBlockP vb, ContextP ctx, const char *unused, unsigned unused2);
extern void sam_set_FLAG_filter (const char *optarg);
extern void sam_set_MAPQ_filter (const char *optarg);

// BAM Stuff
extern void bam_seg_initialize (VBlockP vb);
extern int32_t bam_is_header_done (bool is_eof);
extern int32_t bam_unconsumed (VBlockP vb, uint32_t first_i, int32_t *i);
extern void bam_read_vblock (VBlockP vb);
extern void bam_seg_initialize (VBlockP vb);
extern const char *bam_seg_txt_line (VBlockP vb_, const char *field_start_line, uint32_t remaining_txt_len, bool *has_special_eol);

// SAM-to-FASTQ stuff
CONTAINER_CALLBACK (sam_piz_container_cb);

// VB stuff
extern void sam_vb_release_vb();
extern void sam_vb_destroy_vb();
extern unsigned sam_vb_size (DataType dt);
extern unsigned sam_vb_zip_dl_size (void);

// Special - used for SAM & BAM
#define SAM_SPECIAL { sam_cigar_special_CIGAR, sam_piz_special_TLEN, sam_piz_special_BD_BI, sam_piz_special_SEQ_LEN, \
                      sam_piz_special_MD_old, bam_piz_special_FLOAT, bam_piz_special_BIN, bam_piz_special_NM,   \
                      sam_piz_special_MD, sam_piz_special_REF_CONSUMED, \
                      sam_piz_special_PNEXT_IS_PREV_POS, sam_piz_special_COPY_BUDDY_FLAG, sam_piz_special_COPY_BUDDY_TLEN, \
                      sam_piz_special_COPY_BUDDY_MC, sam_piz_special_CONSUME_MC_Z }
SPECIAL (SAM, 0,  CIGAR,             sam_cigar_special_CIGAR);
SPECIAL (SAM, 1,  TLEN,              sam_piz_special_TLEN);
SPECIAL (SAM, 2,  BDBI,              sam_piz_special_BD_BI);
SPECIAL (SAM, 3,  SEQ_LEN,           sam_piz_special_SEQ_LEN);           // Reconstructs seq_len. Note: called "AS" until 12.0.37 and used to reconstruct AS:i
SPECIAL (SAM, 4,  MD_old,            sam_piz_special_MD_old);            // used in files compressed with Genozip up to 12.0.36
SPECIAL (SAM, 5,  FLOAT,             bam_piz_special_FLOAT);             // used in BAM to represent float optional values
SPECIAL (SAM, 6,  BIN,               bam_piz_special_BIN);   
SPECIAL (SAM, 7,  NM,                bam_piz_special_NM);                // introduced 12.0.37
SPECIAL (SAM, 8,  MD,                sam_piz_special_MD);                // introduced 12.0.37
SPECIAL (SAM, 9,  REF_CONSUMED,      sam_piz_special_REF_CONSUMED);      // introduced 12.0.41: Reconstructs ref_consumed
SPECIAL (SAM, 10, PNEXT_IS_PREV_POS, sam_piz_special_PNEXT_IS_PREV_POS); // introduced 12.0.41
SPECIAL (SAM, 11, COPY_BUDDY_FLAG,   sam_piz_special_COPY_BUDDY_FLAG);   // introduced 12.0.41
SPECIAL (SAM, 12, COPY_BUDDY_TLEN,   sam_piz_special_COPY_BUDDY_TLEN);   // introduced 12.0.41
SPECIAL (SAM, 13, COPY_BUDDY_MC,     sam_piz_special_COPY_BUDDY_MC);     // introduced 12.0.41
SPECIAL (SAM, 14, CONSUME_MC_Z,      sam_piz_special_CONSUME_MC_Z);      // introduced 12.0.41
#define NUM_SAM_SPECIAL 15

#define SAM_LOCAL_GET_LINE_CALLBACKS                      \
    { DT_SAM,  _OPTION_BD_BI,       sam_zip_BD_BI      }, \
    { DT_SAM,  _SAM_QUAL,           sam_zip_qual       }, \
    { DT_SAM,  _OPTION_U2_Z,        sam_zip_U2         }, 

#define BAM_LOCAL_GET_LINE_CALLBACKS                      \
    { DT_BAM,  _OPTION_BD_BI,       sam_zip_BD_BI      }, \
    { DT_BAM,  _SAM_QUAL,           sam_zip_qual       }, \
    { DT_BAM,  _OPTION_U2_Z,        sam_zip_U2         }, 

// Important: Numbers (and order) of translators cannot be changed, as they are part of the file format
// (they are included in the TOP2BAM container)
// translator numbers must start from 1 - 0 is reserved for "none"
TRANSLATOR (SAM, BAM,   1,  I8,         container_translate_I8)   // reconstruct binary little endian functions
TRANSLATOR (SAM, BAM,   2,  U8,         container_translate_U8)   // 
TRANSLATOR (SAM, BAM,   3,  LTEN_I16,   container_translate_LTEN_I16) 
TRANSLATOR (SAM, BAM,   4,  LTEN_U16,   container_translate_LTEN_U16) 
TRANSLATOR (SAM, BAM,   5,  LTEN_I32,   container_translate_LTEN_I32) 
TRANSLATOR (SAM, BAM,   6,  LTEN_U32,   container_translate_LTEN_U32) 
TRANSLATOR (SAM, BAM,   7,  FLOAT,      sam_piz_sam2bam_FLOAT)      // reconstructs SAM-stored textual floating point as little endian 32bit float
TRANSLATOR (SAM, BAM,   8,  ARRAY_SELF, sam_piz_sam2bam_ARRAY_SELF) // remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR (SAM, BAM,   9,  RNAME,      sam_piz_sam2bam_RNAME)      // reconstructs the b250 index or -1 if "*"
TRANSLATOR (SAM, BAM,   10, POS,        sam_piz_sam2bam_POS)        // reconstructs Little Endian U32 0-based POS. 
TRANSLATOR (SAM, BAM,   11, SEQ,        sam_piz_sam2bam_SEQ)        // textual SEQ to BAM-format SEQ 
TRANSLATOR (SAM, BAM,   12, QUAL,       sam_piz_sam2bam_QUAL)       // textual QUAL to BAM-format QUAL 
TRANSLATOR (SAM, BAM,   13, TLEN,       sam_piz_sam2bam_TLEN)       // place TLEN last_value in BAM alignment 
TRANSLATOR (SAM, BAM,   14, OPTIONAL,   sam_piz_sam2bam_OPTIONAL)   // used up to v11, kept for for backward compatability as old files expect it
TRANSLATOR (SAM, BAM,   15, OPTIONAL_SELF, sam_piz_sam2bam_OPTIONAL_SELF) // transform prefixes in Optional Container from SAM to BAM format 
TRANSLATOR (SAM, FASTQ, 16, SEQ,        sam_piz_sam2fastq_SEQ)      // reverse-complement the sequence if needed, and drop if "*"
TRANSLATOR (SAM, FASTQ, 17, QUAL,       sam_piz_sam2fastq_QUAL)     // reverse the QUAL if reverse-complemented and drop fastq records with QUAL="*"
TRANSLATOR (SAM, FASTQ, 18, FLAG,       sam_piz_sam2fastq_FLAG)     // emit 1 if (FLAGS & 0x40) or 2 of (FLAGS & 0x80)

#define NUM_SAM_TRANS   19 // including "none"
#define SAM_TRANSLATORS { NULL /* none */, container_translate_I8, container_translate_U8, container_translate_LTEN_I16, \
                          container_translate_LTEN_U16, container_translate_LTEN_I32, container_translate_LTEN_U32, \
                          sam_piz_sam2bam_FLOAT, sam_piz_sam2bam_ARRAY_SELF, sam_piz_sam2bam_RNAME, sam_piz_sam2bam_POS, sam_piz_sam2bam_SEQ, \
                          sam_piz_sam2bam_QUAL, sam_piz_sam2bam_TLEN, sam_piz_sam2bam_OPTIONAL, sam_piz_sam2bam_OPTIONAL_SELF, \
                          sam_piz_sam2fastq_SEQ, sam_piz_sam2fastq_QUAL, sam_piz_sam2fastq_FLAG }

TXTHEADER_TRANSLATOR (sam_header_bam2sam);
TXTHEADER_TRANSLATOR (sam_header_sam2bam);
TXTHEADER_TRANSLATOR (txtheader_sam2fq);

#define SAM_CONTIG_FMT "@SQ	SN:%.*s	LN:%"PRId64

#endif

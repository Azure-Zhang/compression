// ------------------------------------------------------------------
//   sam_piz.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <math.h>
#include "sam_private.h"
#include "seg.h"
#include "context.h"
#include "piz.h"
#include "reconstruct.h"
#include "strings.h"
#include "dict_id.h"
#include "codec.h" // must be included before reference.h
#include "reference.h"
#include "regions.h"
#include "aligner.h"
#include "file.h"
#include "container.h"
#include "coverage.h"
#include "bases_filter.h"
#include "endianness.h"
#include "lookback.h"
#include "qname.h"

bool sam_piz_read_one_vb (VBlockP vb, Section sl)
{
    if (CTX(OPTION_XA_Z)->dict.len) { // this file has XA
        lookback_init (vb, CTX (OPTION_XA_RNAME),  STORE_INDEX);
        lookback_init (vb, CTX (OPTION_XA_STRAND), STORE_INDEX);
        lookback_init (vb, CTX (OPTION_XA_POS),    STORE_INT);
    }

    return true; // all good
}

// returns true if section is to be skipped reading / uncompressing
bool sam_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    // in collect_coverage, we need only these contexts (+ others, if combined with other filters)
    if (flag.collect_coverage && sections_has_dict_id (st) &&
        (     dict_id.num != _SAM_TOPLEVEL && 
              dict_id.num != _SAM_QNAME    && // only main QNAME context, not its Q?NAME - just for sam_piz_filter to set buddy_line_i
              dict_id.num != _SAM_RNAME    && 
              dict_id.num != _SAM_FLAG     && 
              dict_id.num != _SAM_CIGAR    &&
              dict_id.num != _SAM_BUDDY    && // BUDDY, OPTIONAL and SAM_MC_Z, OPTION_MC_Z are needed for CIGAR
              dict_id.num != _SAM_OPTIONAL && // sam_piz_filter drops it and reconstructs SAM_MC_Z instead
              dict_id.num != _SAM_MC_Z     &&
              dict_id.num != _OPTION_MC_Z  &&
            !(dict_id.num == _SAM_MAPQ     && flag.sam_mapq_filter) && 
            !(dict_id.num == _SAM_SQBITMAP && flag.bases) && 
            !(dict_id.num == _SAM_NONREF   && flag.bases) && 
            !(dict_id.num == _SAM_NONREF_X && flag.bases) && 
            !(dict_id.num == _SAM_GPOS     && flag.bases) && 
            !(dict_id.num == _SAM_STRAND   && flag.bases) && 
            !(dict_id.num == _SAM_TAXID    && flag.kraken_taxid != TAXID_NONE) && 
            !(dict_id.num == _SAM_QNAME    && kraken_is_loaded) && 
            !(dict_id.num == _SAM_BUDDY    && kraken_is_loaded) && // needed by QNAME
            !(dict_id_is_qname_sf(dict_id) && kraken_is_loaded) && 
            !(dict_id.num == _SAM_POS      && flag.regions) && // BUDDY and PNEXT are needed for POS
            !(dict_id.num == _SAM_PNEXT    && flag.regions))) return true;

    // if --count, we only need TOPLEVEL and the fields needed for the available filters (--regions, --FLAG, --MAPQ, --kraken, --taxid)
    if (flag.count && sections_has_dict_id (st)   && 
        (     dict_id.num != _SAM_TOPLEVEL && 
              dict_id.num != _SAM_TOP2BAM  && 
              dict_id.num != _SAM_TOP2FQ   && 
              dict_id.num != _SAM_TOP2FQEX && 
              dict_id.num != _SAM_RNAME    &&  // easier to always have RNAME
            !(flag.out_dt == DT_BAM        && flag.bases)&&  // if output is BAM we need the entire BAM record to correctly analyze the SEQ for IUPAC, as it is a structure.
            !(dict_id.num == _SAM_FLAG     && flag.sam_flag_filter) &&
            !(dict_id.num == _SAM_MAPQ     && flag.sam_mapq_filter) && 
            !(dict_id.num == _SAM_SQBITMAP && flag.bases) && 
            !(dict_id.num == _SAM_NONREF   && flag.bases) && 
            !(dict_id.num == _SAM_NONREF_X && flag.bases) && 
            !(dict_id.num == _SAM_GPOS     && flag.bases) && 
            !(dict_id.num == _SAM_STRAND   && flag.bases) && 
            !(dict_id.num == _SAM_TAXID    && flag.kraken_taxid != TAXID_NONE) && 
            !(dict_id.num == _SAM_QNAME    && kraken_is_loaded) && 
            !(dict_id_is_qname_sf(dict_id) && kraken_is_loaded) && 
            !(dict_id.num == _SAM_POS      && flag.regions) &&
            !(dict_id.num == _SAM_PNEXT    && flag.regions) &&
            !(dict_id.num == _SAM_BUDDY    && flag.regions) // BUDDY and PNEXT are needed for POS
            )) return true;

    return false;
}

// set --FLAG filtering from command line argument
void sam_set_FLAG_filter (const char *optarg)
{
    #define FLAG_ERR "Bad argument of --FLAG: \"%s\". It should be a decimal or hexadecimal integer prefixed with one of + - ^ (+:INCLUDE_IF_ALL ; -:INCLUDE_IF_NONE ; ^:EXCLUDE_IF_ALL"
    switch (optarg[0]) {
        case '+': flag.sam_flag_filter = SAM_FLAG_INCLUDE_IF_ALL  ; break;
        case '-': flag.sam_flag_filter = SAM_FLAG_INCLUDE_IF_NONE ; break;
        case '^': flag.sam_flag_filter = SAM_FLAG_EXCLUDE_IF_ALL  ; break;
        default : ASSERT (false, FLAG_ERR, optarg);
    }

    ASSERT (str_get_int_range_allow_hex16 (&optarg[1], strlen (optarg)-1, 1, 65535, &flag.FLAG), FLAG_ERR, optarg);
}

// set --MAPQ filtering from command line argument
void sam_set_MAPQ_filter (const char *optarg)
{
    #define MAPQ_ERR "Bad argument of --MAPQ: \"%s\". It should be a number 0-255 (INCLUDE lines with MAPQ of at least this) or ^ (eg ^1) (EXCLUDE lines with MAPQ of at least this)"
    
    if (optarg[0] == '^') {
        flag.sam_mapq_filter = SAM_MAPQ_EXCLUDE_IF_AT_LEAST;
        optarg++;
    }
    else 
        flag.sam_mapq_filter = SAM_MAPQ_INCLUDE_IF_AT_LEAST;

    ASSERT (str_get_int_range8 (optarg, 0, 0, 255, &flag.MAPQ), MAPQ_ERR, optarg);
}

static inline void sam_piz_update_coverage (VBlockP vb, const uint16_t sam_flag, uint32_t soft_clip)
{
    ARRAY (uint64_t, read_count, vb->read_count);
    ARRAY (uint64_t, coverage, vb->coverage);
    uint64_t *coverage_special   = AFTERENT (uint64_t, vb->coverage)   - NUM_COVER_TYPES;
    uint64_t *read_count_special = AFTERENT (uint64_t, vb->read_count) - NUM_COVER_TYPES;
    WordIndex chrom_index = vb->last_index(SAM_RNAME);

    if (chrom_index == WORD_INDEX_NONE ||
             sam_flag & SAM_FLAG_UNMAPPED)       { coverage_special[CVR_UNMAPPED]      += vb->seq_len; read_count_special[CVR_UNMAPPED]     ++; }
    else if (sam_flag & SAM_FLAG_FAILED_FILTERS) { coverage_special[CVR_FAILED]        += vb->seq_len; read_count_special[CVR_FAILED]       ++; }
    else if (sam_flag & SAM_FLAG_DUPLICATE)      { coverage_special[CVR_DUPLICATE]     += vb->seq_len; read_count_special[CVR_DUPLICATE]    ++; }
    else if (sam_flag & SAM_FLAG_SECONDARY)      { coverage_special[CVR_SECONDARY]     += vb->seq_len; read_count_special[CVR_SECONDARY]    ++; }
    else if (sam_flag & SAM_FLAG_SUPPLEMENTARY)  { coverage_special[CVR_SUPPLEMENTARY] += vb->seq_len; read_count_special[CVR_SUPPLEMENTARY]++; }
    else {
        coverage_special[CVR_SOFT_CLIP] += soft_clip;
        coverage[chrom_index] += vb->seq_len - soft_clip;
        read_count[chrom_index]++;
    }
}

// For collated files, in lines where PNEXT equal the previous line's POS.
SPECIAL_RECONSTRUCTOR (sam_piz_special_PNEXT_IS_PREV_POS)
{
    new_value->i = CTX(SAM_POS)->last_value.i - CTX(SAM_POS)->last_delta;
    if (reconstruct) RECONSTRUCT_INT (new_value->i);
    
    return true; // new value
}

// Case 1: BIN is set to SPECIAL, we will set new_value here to -1 and wait for CIGAR to calculate it, 
//         as we need vb->ref_consumed - sam_cigar_special_CIGAR will update the reconstruced value
// Case 2: BIN is an textual integer snip - its BIN.last_value will be set as normal and transltor will reconstruct it
SPECIAL_RECONSTRUCTOR (bam_piz_special_BIN)
{
    ctx->semaphore = true; // signal to sam_cigar_special_CIGAR to calculate
    return false; // no new value
}

// note of float reconstruction:
// When compressing SAM, floats are stored as a textual string, reconstruced natively for SAM and via sam_piz_sam2bam_FLOAT for BAM.
//    Done this way so when reconstructing SAM, the correct number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries, encoded as uint32, and stringified to a snip. They are reconstructed,
//    either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. Done this way so BAM binary float is reconstructd precisely.
SPECIAL_RECONSTRUCTOR (bam_piz_special_FLOAT)
{
    int64_t n;
    ASSERT (str_get_int (snip, snip_len, &n), "failed to read integer in %s", ctx->tag_name);
    
    union {
        uint32_t i;
        float f;
    } machine_en = { .i = (uint32_t)n };

    if (!reconstruct) goto finish;

    // binary reconstruction in little endian - BAM format
    if (flag.out_dt == DT_BAM) {
        uint32_t n32_lten = LTEN32 (machine_en.i); // little endian (BAM format)
        RECONSTRUCT (&n32_lten, sizeof (uint32_t)); // in binary - the float and uint32 are the same
    }

    // textual reconstruction - SAM format 
    else { 
        #define NUM_SIGNIFICANT_DIGITS 6 // 6 significant digits, as samtools does
        
        // calculate digits before and after the decimal point
        double log_f = log10 (machine_en.f >= 0 ? machine_en.f : -machine_en.f);
        unsigned int_digits = (log_f >= 0) + (unsigned)log_f;
        unsigned dec_digits = MAX_(0, NUM_SIGNIFICANT_DIGITS - int_digits);
        
        // reconstruct number with exactly NUM_SIGNIFICANT_DIGITS digits
        sprintf (AFTERENT (char, vb->txt_data), "%.*f", dec_digits, machine_en.f); 
        unsigned len = strlen (AFTERENT (char, vb->txt_data)); 
        vb->txt_data.len += len;

        // remove trailing decimal zeros:  "5.500"->"5.5" ; "5.0000"->"5" ; "50"->"50"
        if (dec_digits) {
            unsigned trailing_zeros=0;
            for (int i=vb->txt_data.len-1; i >= vb->txt_data.len-dec_digits; i--)
                if (*ENT (char, vb->txt_data, i) == '0') 
                    trailing_zeros++;
                else
                    break;
            
            vb->txt_data.len -= (dec_digits==trailing_zeros) ? dec_digits+1 : trailing_zeros;
        }
    }

finish:
    new_value->f = (double)machine_en.f;
    return true; // have new value
}

//-----------------------------------------------------------------
// Translator functions for reconstructing SAM data into BAM format
//-----------------------------------------------------------------

// output the word_index of RNAME, which is verified in ref_contigs_get_ref_chrom during seg
// to be the same as the reference id 
TRANSLATOR_FUNC (sam_piz_sam2bam_RNAME)
{
    STR0(snip);
    ctx_get_snip_by_word_index (ctx, ctx->last_value.i, &snip, &snip_len);

    // if it is '*', reconstruct -1
    if (snip_len == 1 && *snip == '*') 
        RECONSTRUCT_BIN32 (-1);

    // if its RNEXT and =, emit the last index of RNAME
    else if (ctx->dict_id.num == _SAM_RNEXT && snip_len == 1 && *snip == '=') 
        RECONSTRUCT_BIN32 (CTX(SAM_RNAME)->last_value.i);

    // otherwise - output the word_index which was stored here because of flags.store=STORE_INDEX set in seg 
    else     
        RECONSTRUCT_BIN32 (ctx->last_value.i); 
    
    return 0;
}

// output, in binary form, POS-1 as BAM uses 0-based POS
TRANSLATOR_FUNC (sam_piz_sam2bam_POS)
{
    RECONSTRUCT_BIN32 (ctx->last_value.i - 1);
    return 0;
}

// translate OPTIONAL SAM->BAM - called as translator-only item on within the Optional reconstruction
// fix prefix eg MX:i: -> MXs
TRANSLATOR_FUNC (sam_piz_sam2bam_OPTIONAL_SELF)
{
    ContainerP con = (ContainerP)recon;

    // if this translator is called due to SAM->FASTQEXT, cancel translation for the OPTIONAL container and its items
    if (flag.out_dt == DT_FASTQ) {
        con->no_translation = true; // turn off translation for this container and its items
        return 0; 
    }

    if (recon_len == -1) return 0; // no Optional data in this alignment

    char *prefixes_before = &recon[con_sizeof (*con)] + 2; // +2 to skip the empty prefixes of container wide, and item[0]
    char *prefixes_after = prefixes_before;

    uint32_t num_items = con_nitems (*con);

    for (unsigned i=1; i < num_items; i++, prefixes_before+=6, prefixes_after+=4) {
        prefixes_after[0] = prefixes_before[0]; // tag[0] 
        prefixes_after[1] = prefixes_before[1]; // tag[1]
        prefixes_after[2] = prefixes_before[3]; // type
        prefixes_after[3] = SNIP_CONTAINER; // end of prefix

        // a SAM 'i' translate to one of several BAM types using the translator code
        // that may be 0->6 (NONE to SAM2BAM_LTEN_U32)
        if (prefixes_after[2] == 'i') 
            prefixes_after[2] = "\0cCsSiI"[con->items[i].translator];
    }

    return -2 * (num_items-1); // change in prefixes_len
}

// translate OPTIONAL SAM->BAM - called after Optional reconstruction is done
TRANSLATOR_FUNC (sam_piz_sam2bam_OPTIONAL)
{
    /* up to v11 we used this translator to set alignment.block_size. due to container logic change, it
       has now moved to sam_piz_container_cb. we keep this translator function for backward
       compatability with v11 bam.genozip files */
    return 0;
}

// note of float reconstruction:
// When compressing SAM, floats are stored as a textual string, reconstruced natively for SAM and via sam_piz_sam2bam_FLOAT for BAM.
//    Done this way so when reconstructing SAM, the correct number of textual digits is reconstructed.
// When compressing BAM, floats are stored as 32-bit binaries, encoded as uint32, and stringified to a snip. They are reconstructed,
//    either as textual for SAM or binary for BAM via bam_piz_special_FLOAT. Done this way so BAM binary float is reconstructd precisely.
TRANSLATOR_FUNC (sam_piz_sam2bam_FLOAT)
{
    union {
        float f; // 32 bit float
        uint32_t i;
    } value;
    
    ASSERT0 (sizeof (value)==4, "expecting value to be 32 bits"); // should never happen

    value.f = (float)ctx->last_value.f;
    RECONSTRUCT_BIN32 (value.i);

    return 0;
}

// remove the comma from the prefix that contains the type, eg "i,"->"i"
TRANSLATOR_FUNC (sam_piz_sam2bam_ARRAY_SELF)
{
    ContainerP con = (ContainerP)recon;
    char *prefixes = &recon[con_sizeof (*con)];

    // remove the ',' from the prefix, and terminate with CON_PX_SEP_SHOW_REPEATS - this will cause
    // the number of repeats (in LTEN32) to be outputed after the prefix
    prefixes[1] = CON_PX_SEP_SHOW_REPEATS; // prefixes is now { type, CON_PX_SEP_SHOW_REPEATS }
    
    return -1; // change in prefixes length
}

//------------------------------------------------------------------------------------
// Translator and filter functions for reconstructing SAM / BAM data into FASTQ format
//------------------------------------------------------------------------------------

TXTHEADER_TRANSLATOR (txtheader_sam2fq)
{
    txtheader_buf->len = 0; // fastq has no header
}

// filtering during reconstruction: called by container_reconstruct_do for each sam alignment (repeat)
CONTAINER_CALLBACK (sam_piz_container_cb)
{
    if (is_top_level) {
        // case SAM to BAM translation: set alignment.block_size (was in sam_piz_sam2bam_OPTIONAL until v11)
        if (dict_id.num == _SAM_TOP2BAM) { 
            BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)ENT (char, vb->txt_data, vb->line_start);
            alignment->block_size = vb->txt_data.len - vb->line_start - sizeof (uint32_t); // block_size doesn't include the block_size field itself
            alignment->block_size = LTEN32 (alignment->block_size);
        }
        
        // case SAM to FASTQ translation: drop line if this is not a primary alignment (don't show its secondary or supplamentary alignments)
        else if ((dict_id.num == _SAM_TOP2FQ || dict_id.num == _SAM_TOP2FQEX)
        && ((uint16_t)vb->last_int(SAM_FLAG) & (SAM_FLAG_SECONDARY | SAM_FLAG_SUPPLEMENTARY)))
            vb->drop_curr_line = "not_primary";

        // --taxid: filter out by Kraken taxid (SAM, BAM, FASTQ)
        if (flag.kraken_taxid && !vb->drop_curr_line 
        && (   (kraken_is_loaded  && !kraken_is_included_loaded (vb, last_txt(vb, SAM_QNAME), vb->last_txt_len (SAM_QNAME)))// +1 in case of FASTQ to skip "@"
            || (!kraken_is_loaded && !kraken_is_included_stored (vb, SAM_TAXID, !flag.collect_coverage && !flag.count)))) 
            vb->drop_curr_line = "taxid";

        // --FLAG
        if (flag.sam_flag_filter && !vb->drop_curr_line ) {

            uint16_t this_sam_flag = (uint16_t)vb->last_int (SAM_FLAG);
        
            bool all_flags_set = (this_sam_flag & flag.FLAG) == flag.FLAG;
            bool no_flags_set  = (this_sam_flag & flag.FLAG) == 0;
            if ((flag.sam_flag_filter == SAM_FLAG_INCLUDE_IF_ALL  && !all_flags_set)
            ||  (flag.sam_flag_filter == SAM_FLAG_INCLUDE_IF_NONE && !no_flags_set)
            ||  (flag.sam_flag_filter == SAM_FLAG_EXCLUDE_IF_ALL  &&  all_flags_set))
                vb->drop_curr_line = "FLAG";
        }

        // --MAPQ
        if (flag.sam_mapq_filter && !vb->drop_curr_line) {
            
            if (dict_id.num == _SAM_TOP2FQ || dict_id.num == _SAM_TOP2FQEX)
                reconstruct_from_ctx (vb, SAM_MAPQ, 0, false); // when translating to FASTQ, MAPQ is normally not reconstructed

            uint8_t this_mapq = (uint8_t)vb->last_int (SAM_MAPQ);
        
            if ((flag.sam_mapq_filter == SAM_MAPQ_INCLUDE_IF_AT_LEAST && this_mapq < flag.MAPQ) ||
                (flag.sam_mapq_filter == SAM_MAPQ_EXCLUDE_IF_AT_LEAST && this_mapq >= flag.MAPQ))
                vb->drop_curr_line = "MAPQ";
        }

        // --bases
        if (flag.bases && !vb->drop_curr_line && 
            !(TXT_DT(DT_BAM) ? iupac_is_included_bam   (last_txt (vb, SAM_SQBITMAP), ((BAMAlignmentFixed *)recon)->l_seq)
                             : iupac_is_included_ascii (last_txt (vb, SAM_SQBITMAP), vb->last_txt_len (SAM_SQBITMAP))))
            vb->drop_curr_line = "bases";
        
        // count coverage, if needed    
        if ((flag.show_sex || flag.show_coverage) && is_top_level && !vb->drop_curr_line)
            sam_piz_update_coverage (vb, vb->last_int(SAM_FLAG), VB_SAM->soft_clip);

        if (flag.idxstats && !vb->drop_curr_line) {
            if (vb->last_int(SAM_FLAG) & SAM_FLAG_UNMAPPED)   
                (*ENT (uint64_t, vb->unmapped_read_count, vb->last_index(SAM_RNAME)))++;
            else
                (*ENT (uint64_t, vb->read_count, vb->last_index(SAM_RNAME)))++;
        }
    }
}

// filter is called before reconstruction of a repeat or an item, and returns false if item should 
// not be reconstructed. contexts are not consumed.
CONTAINER_FILTER_FUNC (sam_piz_filter)
{
    // BAM: set buddy at the beginning of each line, as QNAME is reconstructed much later
    if (dict_id.num == _SAM_TOP2BAM && item == 0 && 
        CTX(SAM_QNAME)->b250.len) { // might be 0 in special cases, like flag.count
        STR(snip);
        PEEK_SNIP(SAM_QNAME);
        if (snip_len && *snip == SNIP_COPY_BUDDY)
            reconstruct_set_buddy(vb);
    }

    // collect_coverage: set buddy_line_i here, since we don't reconstruct QNAME
    else if (dict_id.num == _SAM_QNAME && flag.collect_coverage) {
        STR(snip);
        LOAD_SNIP(SAM_QNAME);
        if (snip_len && *snip == SNIP_COPY_BUDDY)
            reconstruct_set_buddy(vb);
        return false; // don't reconstruct QNAME
    }

    // XA: insert RNAME, STRAND and POS lookbacks
    else if (dict_id.num == _OPTION_XA_Z && item == 4)  // importantly, after RNAME and STRAND_POS
        sam_piz_XA_field_insert_lookback (vb);
    
    // collect_coverage: rather than reconstructing optional, reconstruct SAM_MC_Z that just consumes MC:Z if it exists
    else if (dict_id.num == _SAM_OPTIONAL && flag.collect_coverage) {
        reconstruct_from_ctx (vb, SAM_MC_Z, 0, false);
        return false; // don't reconstruct OPTIONAL
    }

    return true; // go ahead and reconstruct
}

// emit 1 if (FLAGS & 0x40) or 2 of (FLAGS & 0x80) except if --FLAG is specified too (--FLAG can be
// used to get only R1 or R2)
TRANSLATOR_FUNC (sam_piz_sam2fastq_FLAG)
{
    if (flag.sam_flag_filter) return 0;
    
    uint16_t sam_flag = (uint16_t)vb->last_int(SAM_FLAG);

    if (sam_flag & (SAM_FLAG_IS_FIRST | SAM_FLAG_IS_LAST)) {
        
        recon -= item_prefix_len + 1; // move to before prefix and previous field's separator (\n for regular fq or \t for extended)
        memmove (recon+2, recon, recon_len + item_prefix_len + 1); // make room for /1 or /2 
        recon[0] = '/';
        recon[1] = (sam_flag & SAM_FLAG_IS_FIRST) ? '1' : '2';
        vb->txt_data.len += 2; 
    }

    return 0;
}


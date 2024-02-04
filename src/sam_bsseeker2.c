// ------------------------------------------------------------------
//   sam_bsseeker2.c
//   Copyright (C) 2022-2024 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// Compresses auxilliary fields generated by BS-Seeker2

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"

// ----------------------------------------------------------------------------------------------
// XO:Z: BSSeeker2: Orientation, from forward/reverted
// ----------------------------------------------------------------------------------------------

void sam_seg_bsseeker2_XO_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XO), unsigned add_bytes)
{
    // predicting XO to be one of 4 values: +FR -FR +FW -FW, based on a combination of FLAG.rev_comp and multi_segs
    
    // case: value of XO fails the prediction - seg as normal snip     
    ASSSEG (XO_len == 3 && (XO[0]=='+' || XO[0]=='-') && XO[1]=='F' && (XO[2]=='R' || XO[2]=='W'),
            "XO:Z=%.*s but expecting one of four values: +FR -FR +FW -FW", STRf(XO));

    seg_by_did (VB, ((char[]){ SNIP_SPECIAL, SAM_SPECIAL_BSSEEKER2_XO, 
                               vb->bisulfite_strand?'*' : ((XO[0] == '-') == dl->FLAG.rev_comp)?'^' : XO[0], // if we have vb->bisulfite_strand - take for it
                               (XO[2] == 'R') == dl->FLAG.multi_segs ? '*' : XO[2] }), // prediction: 'R' iff FLAG.multi_segs                  
                4, OPTION_XO_Z, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_BSSEEKER2_XO)
{
    if (reconstruct) {
        RECONSTRUCT1 (snip[0] == '*' ? "+-"[VB_SAM->bisulfite_strand=='G'] 
                    : snip[0] == '^' ? "+-"[last_flags.rev_comp] 
                    :                  snip[0]);
        RECONSTRUCT1 ('F');
        RECONSTRUCT1 (snip[1] == '*' ? "WR"[last_flags.multi_segs]         : snip[1]);
    }

    return NO_NEW_VALUE;
}

// ----------------------------------------------------------------------------------------------
// XG:Z: BSSeeker2: genome sequences, with 2bp extended on both ends, from 5' to 3'
// ----------------------------------------------------------------------------------------------

// some files appear to include the left S clip in the calculations, while others don't. We support both options,
// and try to keep the snip consistent in the VB, for better compressability (if we're lucky, it will be determined
// by segconf and remain constant)
static rom sam_seg_XG_Z_analyze_test_lens (VBlockSAMP vb, STRp(XG))
{
    if (XG_len < 6 || XG[2] != '_' || XG[XG_len-3] != '_') return "malformed";

    thool *XG_inc_S = &CTX(OPTION_XG_Z)->XG_inc_S;

    if (vb->soft_clip[0]) {
        if (vb->ref_consumed + vb->soft_clip[0] + 6 == XG_len) 
            *XG_inc_S = yes;
        else if (vb->ref_consumed + 6 == XG_len)
            *XG_inc_S = no;
        else
            return "XG_wrong_length_with_S"; // +6 for the 4 flanking bases and 2 underscore characters.

        if (segconf.running) 
            segconf.sam_XG_inc_S = *XG_inc_S;
    }
    
    // no soft_clip in this line - any value of XG_inc_S will work, we try to pick a consistent value for better compression
    else {
        if (vb->ref_consumed + 6 != XG_len) return "XG_wrong_length_without_S"; 
    
        if (segconf.sam_XG_inc_S != unknown)
            *XG_inc_S = segconf.sam_XG_inc_S;    
        else if (*XG_inc_S == unknown)
            *XG_inc_S = yes; // arbitrary
        // else keep value of previous line
    }

    return NULL; // all good
}

// called for analyzing XG before segging SEQ and later XG
// - Verifies the bases in XG are identical to the reference, or if not in the reference yet (in REF_INTERNAL), adds them
// Example: XG:Z:AC_CCCGTCCCTACTAAAATACAAAAATTAGCCCAGCTTGGTGGTGGGCACCTGTAATCTTAGCTACTGCAGAGACTGAGGCAGGAGAATCGCTTGAACCCAGGAGGTGGAGGTT_GC
void sam_seg_bsseeker2_XG_Z_analyze (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XG), PosType32 line_pos)
{
    START_TIMER;
    decl_ctx (OPTION_XG_Z);

    if (!has(XG_Z)) return;
    
    RefLock lock = REFLOCK_NONE;
    rom result = NULL; // initialize to "success"
    ctx->XG.len = 0;    // initialize value, to remain in case of failure before generating ctx->XG
    int32_t inc_soft_clip = 0;

    #define FAIL(reason) ({ result=reason ; goto done; })

    rom res = sam_seg_XG_Z_analyze_test_lens (vb, STRa(XG));
    if (res) FAIL (res);

    if (segconf.running) return;

    buf_alloc_exact (VB, ctx->XG, XG_len-2/*two underscores removed*/, char, "XG");

    // copy (straight or revcomp), but remove underscores
    if (dl->FLAG.rev_comp) {
        str_revcomp_actg (Bc(ctx->XG, 0), &XG[XG_len-2], 2); 
        str_revcomp_actg (Bc(ctx->XG, 2), &XG[3], XG_len-6); 
        str_revcomp_actg (Bc(ctx->XG, XG_len-4), XG, 2); 
    }        
    else {
        memcpy (Bc(ctx->XG, 0), XG, 2); 
        memcpy (Bc(ctx->XG, 2), &XG[3], XG_len-6); 
        memcpy (Bc(ctx->XG, XG_len-4), &XG[XG_len-2], 2); 
    }

    inc_soft_clip = (ctx->XG_inc_S == yes) ? vb->soft_clip[0] : 0; // soft clip length, if we need to include it

    PosType32 start_pos = line_pos - inc_soft_clip - 2;
    if (start_pos < 1) FAIL("start_pos");

    rom xg = B1STc(ctx->XG);
    PosType32 after_pos = start_pos + vb->ref_consumed + inc_soft_clip + 4/*flanking*/;

    // get range
    RangeP range = ref_seg_get_range (VB, gref, vb->chrom_node_index, STRa(vb->chrom_name), start_pos, after_pos - start_pos, 
                                      WORD_INDEX_NONE,
                                      IS_REF_EXTERNAL ? NULL : &lock); // lock if we might modify reference or is_set

    if (!range) FAIL("no_range"); // either hash contention in REF_INTERNAL or this chromosome is missing in the reference file 
    if (range->last_pos < after_pos-1) FAIL("multi_range"); // sequence spans two ranges - can only happen in REF_INTERNAL

    decl_acgt_decode;
    for (PosType32 pos = start_pos ; pos < after_pos; pos++, xg++) {
        uint32_t pos_index = pos - range->first_pos; // index within range

        bool is_populated = (flag.reference & REF_ZIP_LOADED) || ref_is_nucleotide_set (range, pos_index);

        // case: reference already contains a base - but unfortunately it is not "base" - so MD is not reconstractable from reference
        if (is_populated) { if (*xg != REF (pos_index))    FAIL ("ref_mismatch"); }
        else              { if (*xg!='A' && *xg!='C' && *xg!='G' && *xg!='T') FAIL ("not_ACGT");     }
    }

    // if successful (all bases are the same, or not populated): populate missing bases in REF_INTERNAL
    if (IS_REF_INTERNAL && IS_MAIN(vb)) {
        xg = B1STc(ctx->XG);
        
        for (PosType32 pos = start_pos ; pos < after_pos; pos++, xg++) {
            uint32_t pos_index = pos - range->first_pos; // index within range
            bool is_populated = ref_is_nucleotide_set (range, pos_index);

            // case: reference is not set yet - set it now
            if (!is_populated)
                ref_set_nucleotide (range, pos_index, *xg);
        }
    }

    // set is_set - we will need these bases in the reference to reconstruct XG
    if (flag.reference & REF_STORED && IS_MAIN(vb)) 
        bits_set_region (&range->is_set, start_pos - range->first_pos, after_pos - start_pos); 

    ctx_set_encountered (VB, CTX(OPTION_XG_Z)); // = verified

done:
    if (result && flag.show_wrong_xg)
        iprintf ("%s: RNAME=%.*s POS=%d FLAG=%u CIGAR=\"%s\" ref_consumed%s=%u XG_len-6=%u Special XG not suitable (reason: \"%s\") (no harm)\n", 
                 LN_NAME, STRf(vb->chrom_name), line_pos, dl->FLAG.value, vb->last_cigar, (ctx->XG_inc_S == yes ? "+soft_clip[0]" : ""), vb->ref_consumed + inc_soft_clip, XG_len-6, result);

    ref_unlock (gref, &lock);
    
    COPY_TIMER (sam_seg_bsseeker2_XG_Z_analyze);
    #undef FAIL
}

void sam_seg_bsseeker2_XG_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XG), unsigned add_bytes)
{
    decl_ctx (OPTION_XG_Z);

    if (ctx_encountered_in_line (VB, OPTION_XG_Z)) // encountered in this line & verified
        seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_BSSEEKER2_XG, '1'+ctx->XG_inc_S }, 3, ctx, add_bytes);
    else
        seg_add_to_local_string (VB, ctx, STRa(XG), LOOKUP_SIMPLE, add_bytes);
}

SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_BSSEEKER2_XG)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    if (!reconstruct) goto done;
    decl_acgt_decode;

    // case: XG already reconstructed in sam_piz_special_BSSEEKER2_XM and held in vb->XM
    if (ctx->XG.len) { // note: XG.len is initialized for every line in sam_reset_line
        RECONSTRUCT (Bc(ctx->XG, 0), 2);
        RECONSTRUCT1('_');
        RECONSTRUCT (Bc(ctx->XG, 2), ctx->XG.len-4);
        RECONSTRUCT1('_');
        RECONSTRUCT (Bc(ctx->XG, ctx->XG.len-2), 2);
        goto done;
    }

    bool XG_inc_S = snip[0] - '1';
    int32_t inc_soft_clip = XG_inc_S ? vb->soft_clip[0] : 0; // soft clip length, if we need to include it

    ConstRangeP range = ref_piz_get_range (VB, gref, HARD_FAIL);

    uint32_t idx = CTX(SAM_POS)->last_value.i - range->first_pos - inc_soft_clip;

    // 2 left-flanking bases
    RECONSTRUCT1 (REF (idx-2));
    RECONSTRUCT1 (REF (idx-1));
    RECONSTRUCT1 ('_');

    char *recon = BAFTtxt;
    uint32_t recon_len = vb->ref_consumed + inc_soft_clip;
    for (uint32_t i=0; i < recon_len; i++)
        recon[i] = REF (idx + i);

    Ltxt += recon_len;

    // 2 right-flanking bases
    RECONSTRUCT1 ('_');
    RECONSTRUCT1 (REF (idx + recon_len));
    RECONSTRUCT1 (REF (idx + recon_len + 1));

    if (last_flags.rev_comp)
        str_revcomp_actg (recon-3, recon-3, recon_len + 6);

    done: return NO_NEW_VALUE;
}
 
// ----------------------------------------------------------------------------------------------
// XM:Z: BSSeeker2: number of sites for mismatch
// ----------------------------------------------------------------------------------------------

// X=methylated CG x=un-methylated CG Y=methylated CHG y=un-methylated CHG Z=methylated CHH z=un-methylated CHH
// Example: XM:Z:-------------------------z--z--X----------------------z------------zy--------z--z--z--------yx---------z--z-------z-z---------z------------X---------z--z--X---------

static inline char XM_predict_fwd (BamCigarOpType op, char xg0, char xg1, char xg2, char seq)
{
    if      (op == BC_I || op == BC_D)         return '-'; // insertion - no reference for this base
    else if (xg0=='C' && xg1=='G')             return seq=='C'?'X' : seq=='T'?'x' : '-'; // CG
    else if (xg0=='C' && xg1!='G' && xg2=='G') return seq=='C'?'Y' : seq=='T'?'y' : '-'; // CHG
    else if (xg0=='C' && xg1!='G' && xg2!='G') return seq=='C'?'Z' : seq=='T'?'z' : '-'; // CHH
    else                                       return '-';
}

static inline char XM_predict_rev (BamCigarOpType op, char xg0, char xg1, char xg2, char seq)
{
    if      (op == BC_I || op == BC_D)         return '-'; // insertion - no reference for this base
    else if (xg0=='G' && xg1=='C')             return seq=='G'?'X' : seq=='A'?'x' : '-'; // CG
    else if (xg0=='G' && xg1!='C' && xg2=='C') return seq=='G'?'Y' : seq=='A'?'y' : '-'; // CHG
    else if (xg0=='G' && xg1!='C' && xg2!='C') return seq=='G'?'Z' : seq=='A'?'z' : '-'; // CHH
    else                                       return '-';
}

#define XM_UPDATE_OP \
    if (!op.n) { /* sam_cigar_analyze already verifed that CIGAR is good, so we are going to just assume it here */ \
        op = cigar[op_i++]; \
        if ((op.op == BC_S && op_i!=1) || op.op == BC_H || op.op == BC_P || op.op == BC_N) { \
            op.n = 1; i--; \
            continue; /* skip this op, get next one */ \
        } \
    }

void sam_seg_bsseeker2_XM_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(XM), unsigned add_bytes)
{
    if (segconf.running || vb->cigar_missing || vb->seq_missing) goto fallback;

    ContextP ctx_xg = CTX(OPTION_XG_Z); 

    thool XG_inc_S = ctx_xg->XG_inc_S;

    int32_t inc_soft_clip = (XG_inc_S == yes) ? vb->soft_clip[0] : 0; // soft clip length, if we need to include it
        
    uint32_t expected_xm_len = vb->ref_consumed + dl->SEQ.len - vb->ref_and_seq_consumed;
    if (XG_inc_S == no) expected_xm_len -= vb->soft_clip[0] + vb->soft_clip[1];

    if (expected_xm_len != XM_len)  { // XM includes entries for I,D,M, maybe left-S
        if (flag.show_wrong_xm)
            iprintf ("%s: XM not special bc bad length (no harm): cigar=\"%s\". Expecting: ref_consumed + seq_len - ref_and_seq_consumed=%d == XM_len=%d\n",
                     LN_NAME, vb->last_cigar, vb->ref_consumed + dl->SEQ.len - vb->ref_and_seq_consumed, XM_len);
        goto fallback;
    }

    // note: we don't require XG to be verified, but we require it to be the expected length
    if (vb->ref_consumed + inc_soft_clip != ctx_xg->XG.len - 4) {
        if (flag.show_wrong_xm)
            iprintf ("%s: XM not special bc bad length (no harm): cigar=\"%s\". Expecting: (vb->ref_consumed%s)=%d != (XG.len-4)=%d\n",
                     LN_NAME, vb->last_cigar, (XG_inc_S==yes ? " + vb->soft_clip[0]":""), vb->ref_consumed + inc_soft_clip, (int)ctx_xg->XG.len - 4);
        goto fallback;
    }

    rom xg  = Bc (ctx_xg->XG, 2); // skip 2 flanking bases
    rom seq = vb->textual_seq_str;
    bool rev_comp = dl->FLAG.rev_comp;

    ARRAY (BamCigarOp, cigar, vb->binary_cigar);
    int op_i = (XG_inc_S == no && cigar[0].op == BC_S);
    BamCigarOp op = {};

    for (int32_t i=0; i < XM_len; i++, op.n--) {

        XM_UPDATE_OP;

        char predicted_xm = rev_comp ? XM_predict_rev (op.op, xg[0], xg[-1], xg[-2], *seq)
                                     : XM_predict_fwd (op.op, xg[0], xg[+1], xg[+2], *seq);

        uint32_t xm_i = rev_comp ? XM_len-i-1 : i;
        char this_xm = XM[xm_i];

        if (predicted_xm != this_xm) {
            if (flag.show_wrong_xm)
                iprintf ("%s: XM mis-predicted (no harm): xm_i=%u revcomp=%u multi=%u: cigar=%s op=%c xg=%c seq=%c predicted_xm=%c xm=%c XM=%.*s\n",
                         LN_NAME, xm_i, rev_comp, dl->FLAG.multi_segs, vb->last_cigar, cigar_op_to_char[op.op], *xg, *seq, predicted_xm, this_xm, STRf(XM));
            goto fallback;
        }

        if (op.op != BC_I) xg++;
        if (op.op != BC_D) seq++;
    }

    // we include the argument XG_inc_S because we might have this special even if XG failed verification and has no special
    seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_BSSEEKER2_XM, '1' + XG_inc_S }, 3, OPTION_XM_Z, add_bytes);
    return;

fallback:
    seg_add_to_local_string (VB, CTX(OPTION_XM_Z), STRa(XM), LOOKUP_SIMPLE, add_bytes);
}

SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_BSSEEKER2_XM)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    if (!reconstruct) goto done;

    bool XG_inc_S = snip[0] - '1';

    STR(xg);
    reconstruct_peek (VB, CTX (OPTION_XG_Z), pSTRa(xg));

    bool rev_comp = last_flags.rev_comp;

    // case: reconstructed to txt_data: copy data peeked 1. because our recon will overwrite it 2. so we don't need to reconstruct it again
    BufferP xg_buf = &CTX(OPTION_XG_Z)->XG;
    buf_alloc_exact (vb, *xg_buf, xg_len-2, char, "xg"); // remove the two underscores
    memcpy (Bc (*xg_buf, 0), xg, 2);
    memcpy (Bc (*xg_buf, 2), xg+3, xg_len - 6);
    memcpy (Bc (*xg_buf, xg_buf->len-2), xg+xg_len-2, 2);
    xg = B1STc(*xg_buf);
    xg_len = xg_buf->len;
    int32_t xg_i=0;

    rom seq = sam_piz_get_textual_seq (VB);
    char *recon = BAFTtxt;
    ARRAY (BamCigarOp, cigar, vb->binary_cigar);
    int op_i = (!XG_inc_S && cigar[0].op == BC_S);
    BamCigarOp op = {};

    uint32_t xm_len = vb->ref_consumed + vb->seq_len - vb->ref_and_seq_consumed;
    if (!XG_inc_S) xm_len -= vb->soft_clip[0] + vb->soft_clip[1];

    for (int32_t i=0; i < xm_len; i++, op.n--) {

        XM_UPDATE_OP

        #define XG_FWD(x) xg[(x)+2] 
        #define XG_REV(x) complem((int)xg[1+(xg_len-4)-(x)])
                
        char xm = rev_comp ? XM_predict_rev (op.op, XG_REV(xg_i), XG_REV(xg_i-1), XG_REV(xg_i-2), *seq)
                           : XM_predict_fwd (op.op, XG_FWD(xg_i), XG_FWD(xg_i+1), XG_FWD(xg_i+2), *seq);

        recon[rev_comp ? xm_len-i-1 : i] = xm;

        if (op.op != BC_I) xg_i++;
        if (op.op != BC_D) seq++;
    }

    Ltxt += xm_len;

    done: return NO_NEW_VALUE;
}
  
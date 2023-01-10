// ------------------------------------------------------------------
//   sam_bismark.c
//   Copyright (C) 2022-2023 Genozip Limited. Patent pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited,
//   under penalties specified in the license.

// Compresses auxilliary fields generated by Bismark

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"

// Bismark: The converted reference used to align this read: 2 possible values: "CT" or "GA"
void sam_seg_bismark_XG_Z (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(xg), unsigned add_bytes)
{
    ASSSEG (xg_len==2 && ((xg[0]=='C' && xg[1]=='T') || (xg[0]=='G' && xg[1]=='A')), xg, "Invalid XG:Z=%.*s, expecting CT or GA", STRf(xg));

    if (vb->bisulfite_strand)
        seg_by_did (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_BISMARK_XG }, 2, OPTION_XG_Z, add_bytes);

    else 
        seg_by_did (VB, STRa(xg), OPTION_XG_Z, add_bytes);
}

SPECIAL_RECONSTRUCTOR (sam_piz_special_BISMARK_XG)
{
    ASSPIZ0 (VB_SAM->bisulfite_strand, "XG:Z cannot be reconstructed because bisulfite_strand is not known - likely SEQ was not reconstructed");
    
    if (reconstruct) 
        RECONSTRUCT (VB_SAM->bisulfite_strand == 'C' ? "CT" : "GA", 2);

    return NO_NEW_VALUE;
}

// enter methylatble bases into the INTERNAL reference in their unconverted form 
// (not currently used as bisulfite features are disabled for REF_INTERNAL (bug 648), and not thoroughly tested)
void sam_seg_bismark_XM_Z_analyze (VBlockSAMP vb, ZipDataLineSAM *dl)
{
    if (!IS_REF_INTERNAL || // analyzing sets bases in an internal reference - not needed if not internal
        has_MD ||           // analyzing MD sets the same bases
        !has_XM || !vb->bisulfite_strand || vb->comp_i != SAM_COMP_MAIN) return;

    STR(xm);
    sam_seg_get_aux_Z (vb, vb->idx_XM_Z, pSTRa(xm), IS_BAM_ZIP);
    uint32_t xm_i = 0;

    RangeP range = NULL;
    RefLock lock = REFLOCK_NONE;
    uint32_t ref_consumed = vb->ref_consumed; // M/=/X and D
    SamPosType pos = dl->POS;

    for_buf (BamCigarOp, op, vb->binary_cigar) 
        if (op->op==BC_M || op->op==BC_E || op->op==BC_X) 
            for (uint32_t i=0; i < op->n; i++) {
                if (xm[xm_i++] != '.') {
                    rom error = sam_seg_analyze_set_one_ref_base (vb, false, pos, vb->bisulfite_strand, ref_consumed, &range, &lock); 
                    if (error == ERR_ANALYZE_RANGE_NOT_AVAILABLE) return; // possibly pos/ref_consumed go beyond end of range
                }
                pos++;
                ref_consumed--;
            }
        

        else if (op->op == BC_D || op->op == BC_N) {
            pos += op->n;
            ref_consumed -= op->n;
        }
        
        else if (op->op == BC_I || op->op == BC_S)
            xm_i += op->n;

    ASSSEG (xm_i == xm_len, xm, "Mismatch between XM:Z=\"%.*s\" and CIGAR=\"%s\"", STRf(xm), 
            dis_binary_cigar (vb, B1ST(BamCigarOp, vb->binary_cigar), vb->binary_cigar.len32, &vb->scratch).s);

    if (range) ref_unlock (gref, &lock);
}

typedef enum { XM_AS_PREDICTED, XM_DIFF, XM_IN_LOCAL } XmSnip; // v14 - part of the file format

// Z/z=methylated/unmethylated CpG ; X/x=CHG ; H/h=CHH ; U/u=undetermined methylation type
void sam_seg_bismark_XM_Z (VBlockSAMP vb, ZipDataLineSAM *dl, Did did_i, int special_code, STRp(xm), unsigned add_bytes)
{
    START_TIMER

    XmSnip xm_type = XM_AS_PREDICTED; // optimistic
    ContextP ctx = CTX(did_i);

    if (segconf.running) goto no_diff;

    if (!str_issame_(STRa(xm), STRb(vb->meth_call))) {
        ASSSEG (xm_len == dl->SEQ.len, xm, "Expecting XM:Z.len=%u == SEQ.len=%u", xm_len, dl->SEQ.len);

        if (flag.show_wrong_xm) {
            iprintf ("%s: QNAME=\"%.*s\" bisulfite_strand=%c %s\n", 
                     LN_NAME, STRfw(dl->QNAME), vb->bisulfite_strand, vb->meth_call.len ? "" : "(no meth_call)");
            iprintf ("XM: %.*s\n", xm_len, xm);
            if (vb->meth_call.len) iprintf ("SQ: %.*s\n", STRfb(vb->meth_call));
        }

        if (vb->meth_call.len) {
            ASSSEG (xm_len == vb->meth_call.len, xm, "Expecting XM:Z.len=%u == meth_call.len=%u", xm_len, vb->meth_call.len32);
            
            // note: our diff is compact, according to actual difference observed in Bismark and BSBolt:
            // z,x,h,u - same case, different type
            // ^ - same type, different case
            // We diff only bases that in our prediction are non-. 
            ARRAY (char, meth_call, vb->meth_call);
            char *next = meth_call;
            for (uint32_t i=0; i < xm_len; i++) {
                char c = meth_call[i];  // prediction
                char x = xm[i];         // actual
            
                if (c == '.' && x != '.') goto no_diff; //  we don't diff if a site is non-methylatable according to our prediction, but yet contains a non-.

                else if (c == '.') continue;

                else if (x == c) *next++ = 0;
                    
                else if (x == '.') *next++ = '.';

                else if (LOWER_CASE(x) == LOWER_CASE(c)) *next++ = '^'; // c and x differ in case, but same type

                else if (IS_SLETTER(x) == IS_SLETTER(c)) *next++ = LOWER_CASE(x); // c and x are same case, but differ in type

                else goto no_diff; // differ both in case and in type - not supported for diff
            }
                                                                  
            xm_type = XM_DIFF;
            seg_add_to_local_fixed (VB, ctx, meth_call, next - meth_call, LOOKUP_NONE, 0);
        }
        
        else no_diff: {
            xm_type = XM_IN_LOCAL;
            seg_add_to_local_fixed (VB, ctx, STRa(xm), LOOKUP_NONE, 0);
        }
    }

    seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, special_code, '0' + xm_type }, 3, ctx, add_bytes);

    COPY_TIMER(sam_seg_bismark_XM_Z);
}

// Z/z=methylated/unmethylated CpG ; X/x=CHG ; H/h=CHH ; U/u=undetermined methylation type
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_BISMARK_XM)
{
    XmSnip xm_type = snip[0] - '0';
    VBlockSAMP vb = (VBlockSAMP)vb_;

    switch (xm_type) {
        case XM_AS_PREDICTED:
            RECONSTRUCT_BUF (vb->meth_call);
            break;

        case XM_IN_LOCAL:
            RECONSTRUCT_NEXT (ctx, vb->seq_len);
            break;

        case XM_DIFF: {
            char *diff = Bc(ctx->local, ctx->next_local);
            ARRAY (char, meth_call, vb->meth_call);
            char *recon = BAFTtxt;

            // note: in our diff, we only consider the potential methlation sites, '.'s are never diffed
            for (uint32_t i=0; i < meth_call_len; i++) {
                char c = meth_call[i], d=0;
                
                if (c == '.') 
                    *recon++ = '.';
                
                else if (!(d = *diff++)) // prediction is correct
                    *recon++ = c;

                else if (d == '.')
                    *recon++ = '.';

                else if (IS_SLETTER(d))  // differ in type, but same case 
                    *recon++ = IS_SLETTER (c) ? d : UPPER_CASE(d); 

                else if (d == '^')       // differ in case, but sames type
                    *recon++ = IS_SLETTER (c) ? UPPER_CASE(c) : LOWER_CASE(c); 

                else
                    ASSPIZ0 (false, "Corrupt XM:Z diff");
            }

            ctx->next_local = BNUM (ctx->local, diff);
            vb->txt_data.len32 += meth_call_len;
    
            break;
        }

        default:
            ASSPIZ (false, "Invalid XM snip=%u", xm_type);
    }

    return NO_NEW_VALUE; 
}

void sam_bismark_zip_update_meth_call (VBlockSAMP vb, RangeP range, uint32_t range_len, int32_t idx, bool methylated,
                                       uint32_t M_i, uint32_t Mseg_len, // the position in the M segment and its length
                                       uint32_t i) // our position within the M segment
{
    if (vb->bisulfite_strand == 'C') { // C->T
        bool has_next = vb->binary_cigar.next < vb->binary_cigar.len;
        BamCigarOp next = has_next ? *B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next) : (BamCigarOp){};
        if (next.op == BC_S || next.op == BC_H) has_next = false;

        bool nxt_is_I = has_next && next.op == BC_I;
        bool nxt_is_D = has_next && next.op == BC_D;
        
        int32_t idx1 = (i < Mseg_len-1 || !has_next) ? RR_IDX (idx + 1) 
                     : nxt_is_D                      ? RR_IDX (idx + 1 + next.n) 
                     :                                 0;

        int32_t idx2 = (i < Mseg_len-2 || !has_next) ? RR_IDX (idx + 2) 
                     : nxt_is_D                      ? RR_IDX (idx + 2 + next.n) 
                     :                                 0;

        // for simplicity, ignoring edge cases:
        //    1D and idx2: xM1D1M1D - idx2 needs to skip 2 Ds
        //                 xM1D1M1I - check for idx2 results in u/U due to the I

        *Bc (vb->meth_call, M_i + i) = 
            (i == Mseg_len-1 && nxt_is_I) ? (methylated ? 'U' : 'u') :
            (REF(idx1) == 'G'           ) ? (methylated ? 'Z' : 'z') :
            (i == Mseg_len-2 && nxt_is_I) ? (methylated ? 'U' : 'u') :
            (REF(idx2) == 'G'           ) ? (methylated ? 'X' : 'x') :
                                            (methylated ? 'H' : 'h') ;
    }

    else { // G->A
        bool has_prev = vb->binary_cigar.next >= 2;
        BamCigarOp prev = has_prev ? *B(BamCigarOp, vb->binary_cigar, vb->binary_cigar.next-2) : (BamCigarOp){};
        if (prev.op == BC_S || prev.op == BC_H) has_prev = false;

        bool prv_is_I = has_prev && prev.op == BC_I;
        bool prv_is_D = has_prev && prev.op == BC_D;

        int32_t idx1 = (i >= 1 || !has_prev) ? RR_IDX (idx - 1) 
                     : prv_is_D              ? RR_IDX (idx - 1 - prev.n) 
                     :                         0;

        int32_t idx2 = (i >= 2 || !has_prev) ? RR_IDX (idx - 2) 
                     : prv_is_D              ? RR_IDX (idx - 2 - prev.n) 
                     :                         0;

        // for simplicity, ignoring edge cases:
        //   xD1M1DxM - idx2 needs to skip 2 Ds
        //   xI1M1DxM - check for idx2 results in u/U due to the I
        
        *Bc (vb->meth_call, M_i + i) = 
            (i == 0 && prv_is_I) ? (methylated ? 'U' : 'u') :
            (REF(idx1) == 'C' )  ? (methylated ? 'Z' : 'z') :
            (i == 1 && prv_is_I) ? (methylated ? 'U' : 'u') :
            (REF(idx2) == 'C' )  ? (methylated ? 'X' : 'x') :
                                   (methylated ? 'H' : 'h') ;
    }
}

void sam_bismark_piz_update_meth_call (VBlockSAMP vb, bytes ref, int32_t idx, uint32_t seq_i, 
                                              uint32_t M_i, uint32_t Mseg_len, // our position in the M segment and its length
                                              const BamCigarOp *cigar, char bisulfite, bool methylated)
{
    START_TIMER;
    
    // prediction logic, precisely mirroring sam_seg_bisulfite_M 
    
    if (bisulfite == 'C') { // C->T
        bool has_next = cigar < BAFT(BamCigarOp, vb->binary_cigar) && cigar->op != BC_S && cigar->op != BC_H; // "cigar" is the next cigar
        bool nxt_is_I = has_next && cigar->op == BC_I;
        bool nxt_is_D = has_next && cigar->op == BC_D;

        // ref 2 bases before and after (see sam_reconstruct_SEQ_get_ref_bytemap), so no need to round-robin
        int32_t idx1 = (M_i < Mseg_len-1 || !has_next) ? idx + 1
                     : nxt_is_D                        ? idx + 1 + cigar->n 
                     :                                   0;

        int32_t idx2 = (M_i < Mseg_len-2 || !has_next) ? idx + 2 
                     : nxt_is_D                        ? idx + 2 + cigar->n 
                     :                                   0;

        *Bc (vb->meth_call, seq_i) = 
            (M_i == Mseg_len-1 && nxt_is_I) ? (methylated ? 'U' : 'u') :
            (ref[idx1] == 2/*G*/          ) ? (methylated ? 'Z' : 'z') :
            (M_i == Mseg_len-2 && nxt_is_I) ? (methylated ? 'U' : 'u') :
            (ref[idx2] == 2/*G*/          ) ? (methylated ? 'X' : 'x') :
                                              (methylated ? 'H' : 'h') ;
    }

    else { // G->A
        cigar -= 2; // cigar is one before current cigar
        bool has_prev = cigar >= B1ST (BamCigarOp, vb->binary_cigar) && cigar->op != BC_S && cigar->op != BC_H;
        bool prv_is_I = has_prev && cigar->op == BC_I;
        bool prv_is_D = has_prev && cigar->op == BC_D;

        int32_t idx1 = (M_i >= 1 || !has_prev) ? idx - 1 
                     : prv_is_D                ? idx - 1 - cigar->n 
                     :                           0;

        int32_t idx2 = (M_i >= 2 || !has_prev) ? idx - 2 
                     : prv_is_D                ? idx - 2 - cigar->n 
                     :                           0;

        *Bc (vb->meth_call, seq_i) = 
            (M_i == 0 && prv_is_I) ? (methylated ? 'U' : 'u') :
            (ref[idx1] == 1/*C*/ ) ? (methylated ? 'Z' : 'z') :
            (M_i == 1 && prv_is_I) ? (methylated ? 'U' : 'u') :
            (ref[idx2] == 1/*C*/ ) ? (methylated ? 'X' : 'x') :
                                     (methylated ? 'H' : 'h') ;
    }

    COPY_TIMER (sam_bismark_piz_update_meth_call);
}

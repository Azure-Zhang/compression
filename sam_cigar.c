// ------------------------------------------------------------------
//   sam_cigar.c
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

// a module for handling CIGAR and MC:Z

#include "genozip.h"
#include "sam_private.h"
#include "strings.h"
#include "reference.h"
#include "segconf.h"
#include "seg.h"
#include "piz.h"
#include "reconstruct.h"
#include "codec.h"
#include "profiler.h"
#include "md5.h"
#include "htscodecs/rANS_static4x16.h"

static const bool cigar_valid_op[256] = { ['M']=true, ['I']=true, ['D']=true, ['N']=true, ['S']=true, ['H']=true, ['P']=true, ['=']=true, ['X']=true }; 

const char cigar_op_to_char[16] = "MIDNSHP=Xabcdefg"; // BAM to SAM (a-g are invalid value)

static const uint8_t cigar_char_to_op[96] = { ['M']=BC_M, ['I']=BC_I, ['D']=BC_D, ['N']=BC_N, ['S']=BC_S, 
                                              ['H']=BC_H, ['P']=BC_P, ['=']=BC_E, ['X']=BC_X, ['*']=BC_NONE }; 

#undef S 
#define S (c == 'S')
#define H (c == 'H')
#define I (c == 'I')
#define D (c == 'D')
#define M (c == 'M')
#define N (c == 'N')
#define P (c == 'P')
#define E (c == '=')
#define X (c == 'X')

//---------
// Shared
//---------

void sam_cigar_S_to_H (STRc(cigar))
{
    // replace left clipping - find first op
    char *c=cigar; while (IS_DIGIT(*c)) c++;
    if (*c == 'S') *c = 'H';

    // replace right clipping
    if (cigar[cigar_len-1] == 'S') cigar[cigar_len-1] = 'H';
}

void sam_cigar_H_to_S (STRc(cigar))
{
    // replace left clipping - find first op
    char *c=cigar; while (IS_DIGIT(*c)) c++;
    if (*c == 'H') *c = 'S';

    // replace right clipping
    if (cigar[cigar_len-1] == 'H') cigar[cigar_len-1] = 'S';
}

// gets seq_len implied by cigar or squanked cigar segments: "M24S" "M14S" "S" "". 
uint32_t sam_cigar_get_seq_len (STRp(cigar))
{
    uint32_t n=0, seq_len=0;

    for (uint32_t i=0; i < cigar_len; i++) {
        char c = cigar[i];
        if (IS_DIGIT(c)) 
            n = n*10 + c-'0';

        else { // op
            if (M || I || S || E || X) seq_len += n;
            n = 0;
        }
    }
    return seq_len;
}

void sam_cigar_binary_to_textual (VBlockSAMP vb, uint16_t n_cigar_op, const uint32_t *cigar, Buffer *textual_cigar /* out */)
{
    if (!n_cigar_op) {
        buf_alloc (vb, textual_cigar, 2, 0, char, 100, textual_cigar->name ? NULL : "textual_cigar");
        BNXTc (*textual_cigar) = '*';
        goto finish;
    }

    // calculate length
    uint32_t len=0;
    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t op_len = LTEN32 (cigar[i]) >> 4; // maximum is 268,435,455
        if      (op_len < 10)       len += 2; // 1 for the op, 1 for the number
        else if (op_len < 100)      len += 3; 
        else if (op_len < 1000)     len += 4; 
        else if (op_len < 10000)    len += 5; 
        else if (op_len < 100000)   len += 6; 
        else if (op_len < 1000000)  len += 7;
        else if (op_len < 10000000) len += 8;
        else                        len += 9;
    }

    buf_alloc (vb, textual_cigar, len + 1 /* for \0 */, 100, char, 0, textual_cigar->name ? NULL : "textual_cigar");

    for (uint16_t i=0; i < n_cigar_op; i++) {
        uint32_t subcigar = LTEN32 (cigar[i]);
        uint32_t op_len = subcigar >> 4;

        textual_cigar->len += str_int (op_len, BAFTc (*textual_cigar));
        BNXTc (*textual_cigar) = cigar_op_to_char[subcigar & 0xf];
    }

finish:
    *BAFTc (*textual_cigar) = 0; // nul terminate
    
    if (command == ZIP)
        vb->last_cigar = B1STc (*textual_cigar);
}
/*
void sam_cigar_textual_to_binary (VBlockSAMP vb, STRp(cigar))
{
   

    if (cigar_len==1 && *cigar == '*') return; // empty binary cigar

    uint32_t n_ops = 0;
    for (int i=0; i < cigar_len; i++)
        if (!IS_DIGIT (cigar[i])) n_ops++;

    ARRAY_alloc (BamCigarOp, ops, n_ops, false, vb->binary_cigar, vb, "binary_cigar");

    for (int op_i=0; op_i < n_ops; op_i++) {
        uint32_t n=0; // do arith on proper integer, not bit field
        while (IS_DIGIT(*cigar)) { n = 10*n + *cigar - '0' ; cigar++; }
        ops[op_i] = (BamCigarOp) { .op = cigar_char_to_op[*cigar++], .n = n } ;
    }
}
*/
// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
void sam_cigar_analyze (VBlockSAMP vb, STRp(cigar)/* textual */, bool cigar_is_in_textual_cigar, uint32_t *seq_consumed)
{
    *seq_consumed = 0; // everything else is initialized in sam_reset_line
    
    ASSERT (cigar[0] != '*' || cigar_len == 1, "Invalid CIGAR: %.*s", cigar_len, cigar); // a CIGAR start with '*' must have 1 character

    // ZIP case: if the CIGAR is "*", later sam_cigar_seg_textual uses the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*", so this condition is always false
    if (cigar[0] == '*') {
        vb->cigar_missing = true;
        return;
    }

    // PIZ case: CIGAR string starts with '-' (indicating missing SEQ) - just skip the '-' for now
    else if (*cigar == '-') {
        vb->seq_missing = true;
        cigar++;
        cigar_len--;
    }

    // store original textual CIGAR for use of sam_piz_special_MD, as in BAM it will be translated ; also cigar might point to buddy data in ctx->per_line - ctx->per_line might be realloced as we store this line's CIGAR in it 
    if (command == PIZ && !cigar_is_in_textual_cigar) {
        buf_add_more (VB, &vb->textual_cigar, STRa(cigar), "textual_cigar");
        *BAFTc (vb->textual_cigar) = 0; // nul-terminate (buf_add_more allocated space for it)
    }

    // create the BAM-style cigar data in binary_cigar. 
    buf_alloc (vb, &vb->binary_cigar, 0, cigar_len/2 /* max possible n_cigar_op */, BamCigarOp, 2, "binary_cigar");
    ARRAY (BamCigarOp, bam_ops, vb->binary_cigar);

    uint32_t n=0;
    uint32_t op_i=0;
    for (uint32_t i=0; i < cigar_len; i++) {

        char c = cigar[i];

        if (IS_DIGIT(c)) 
            n = n*10 + (c - '0');

        else {
            ASSINP (n, "Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=\"%.*s\"", 
                    txt_name, c, cigar_len, cigar);    

            #define SEQ_CONSUMED     *seq_consumed            += n
            #define REF_CONSUMED     vb->ref_consumed         += n
            #define SEQ_REF_CONSUMED vb->ref_and_seq_consumed += n
            #define COUNT(x)         vb->x                    += n
            // Bug 546: Per SAM spec, a CIGAR may have H clips enclosing S clips eg 2H3S5M4S2H. 
            #define LAST_OR_1ST ({ ASSERT(!op_i || i==cigar_len-1, "'%c' can only appear as the first or last op in the CIGAR string. cigar=\"%.*s\"", c, cigar_len, cigar); })

            switch (c) { 
                case 'M' : 
                case '=' : 
                case 'X' : SEQ_CONSUMED ; REF_CONSUMED ; SEQ_REF_CONSUMED          ; break ;
                case 'I' : SEQ_CONSUMED ; COUNT(mismatch_bases)                    ; break ;
                case 'D' : REF_CONSUMED ; COUNT(mismatch_bases)                    ; break ;
                case 'N' : REF_CONSUMED                                            ; break ;
                case 'S' : LAST_OR_1ST ; SEQ_CONSUMED ; COUNT(soft_clip[op_i > 0]) ; break ; // Note: a "121S" (just one op S or H) is considered a left-clip (eg as expected by sam_seg_XG_Z_analyze)
                case 'H' : LAST_OR_1ST ; COUNT (hard_clip[op_i > 0])               ; break ;
                case 'P' :                                                           break ;
                case '*' : SEQ_CONSUMED                                            ; break ; // Note: '*' is when CIGAR is "151*" (PIZ only) - alignment with no CIGAR but a SEQ
                default  : ASSINP (false, "Invalid CIGAR in %s: invalid operation '%c'(ASCII %u). CIGAR=\"%.*s\"", txt_name, c, c, cigar_len, cigar);
            }

            // convert character CIGAR op to BAM cigar field op  "MIDNSHP=X" -> 012345678 ; * is our private value of BC_NONE
            bam_ops[op_i++] = (BamCigarOp){ .op = cigar_char_to_op[(uint8_t)c], .n = n };

            n = 0;
        }
    }          

    // note: piz: in case of eg "151*" - *seq_consumed will be updated to the length, but binary_cigar will be empty
    vb->binary_cigar.len = op_i - (op_i==1 && bam_ops[0].op == BC_NONE);  // 0 if CIGAR is "*" 

    if (!vb->binary_cigar.len) 
        vb->cigar_missing = true;

    if (command == ZIP) 
        DATA_LINE (vb->line_i)->ref_consumed = vb->ref_consumed; // consumed by sam_seg_predict_TLEN 

    // PIZ reconstructing: we store ref_consumed in ctx->piz_ctx_specific_buf because ctx->history is already taken for storing the CIGAR string
    else if (!vb->preprocessing) {
        uint64_t line_i = vb->line_i - vb->first_line;
        if (!line_i)
            buf_alloc_zero (vb, &CTX(SAM_CIGAR)->piz_ctx_specific_buf, 0, vb->lines.len, uint32_t, 0, "piz_ctx_specific_buf"); // initialize to exactly one per line.

        *B32 (CTX(SAM_CIGAR)->piz_ctx_specific_buf, line_i) = vb->ref_consumed;
    }
    
    ASSINP (!n, "Invalid CIGAR in %s: expecting it to end with an operation character. CIGAR=\"%.*s\"", 
            txt_name, cigar_len, cigar);

    ASSINP (!seq_consumed || *seq_consumed, "Invalid CIGAR in %s: CIGAR implies 0-length SEQ. CIGAR=\"%.*s\"", txt_name, cigar_len, cigar);

    // note: if there's ever any need to support both S and H, we need to fix sam_sa_seg_depn_find_sagroup too
    ASSINP (!((vb->hard_clip[0] || vb->hard_clip[1]) && (vb->soft_clip[0] || vb->soft_clip[1])),
            "Invalid CIGAR in %s: has both S and H. CIGAR=\"%.*s\"", txt_name, cigar_len, cigar);
}

// analyze the binary cigar when segging BAM - faster
void bam_seg_cigar_analyze (VBlockSAMP vb, uint32_t *seq_consumed)
{
    *seq_consumed = 0; // everything else is initialized in sam_reset_line
    ARRAY (BamCigarOp, cigar, vb->binary_cigar);

    // ZIP case: if the CIGAR is "*", later sam_cigar_seg_textual uses the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*", so this condition is always false
    if (!cigar_len) {
        vb->cigar_missing = true;
        return;
    }

    for (uint32_t op_i=0; op_i < cigar_len; op_i++) {

        #define SEQ_CONSUMED_     *seq_consumed            += cigar[op_i].n
        #define REF_CONSUMED_     vb->ref_consumed         += cigar[op_i].n
        #define SEQ_REF_CONSUMED_ vb->ref_and_seq_consumed += cigar[op_i].n
        #define COUNT_(x)         vb->x                    += cigar[op_i].n
        // Bug 546: Per SAM spec, a CIGAR may have H clips enclosing S clips eg 2H3S5M4S2H. 
        #undef LAST_OR_1ST
        #define LAST_OR_1ST ({ ASSERT(!op_i || op_i==cigar_len-1, "'%c' can only appear as the first or last op in the CIGAR string. cigar=\"%.*s\"", \
                                      cigar_op_to_char[cigar[op_i].op], (int)vb->textual_cigar.len, vb->textual_cigar.data); })

        switch (cigar[op_i].op) { 
            case BC_M : 
            case BC_E : 
            case BC_X : SEQ_CONSUMED_ ; REF_CONSUMED_ ; SEQ_REF_CONSUMED_         ; break ;
            case BC_I : SEQ_CONSUMED_ ; COUNT_(mismatch_bases)                    ; break ;
            case BC_D : REF_CONSUMED_ ; COUNT_(mismatch_bases)                    ; break ;
            case BC_N : REF_CONSUMED_                                             ; break ;
            case BC_S : LAST_OR_1ST ; SEQ_CONSUMED_ ; COUNT_(soft_clip[op_i > 0]) ; break ; // Note: a "121S" (just one op S or H) is considered a left-clip (eg as expected by sam_seg_XG_Z_analyze)
            case BC_H : LAST_OR_1ST ; COUNT_ (hard_clip[op_i > 0])                ; break ;
            case BC_P :                                                             break ;
            default   : ASSINP (false, "Invalid CIGAR in %s vb=%u line=%u: invalid operation %u", txt_name, vb->vblock_i, (uint32_t)vb->line_i, cigar[op_i].op);
        }
    }          

    DATA_LINE (vb->line_i)->ref_consumed = vb->ref_consumed; // consumed by sam_seg_predict_TLEN 

    ASSINP (!seq_consumed || *seq_consumed, "Invalid CIGAR in %s vb=%u line=%u: CIGAR implies 0-length SEQ CIGAR=\"%.*s\"", 
            txt_name, vb->vblock_i, (uint32_t)vb->line_i, (int)vb->textual_cigar.len, vb->textual_cigar.data);

    // note: if there's ever any need to support both S and H, we need to fix sam_sa_seg_depn_find_sagroup too
    ASSINP (!((vb->hard_clip[0] || vb->hard_clip[1]) && (vb->soft_clip[0] || vb->soft_clip[1])),
            "Invalid CIGAR in %s: has both S and H. CIGAR=\"%.*s\"", txt_name, (int)vb->textual_cigar.len, vb->textual_cigar.data);
}

bool sam_cigar_is_valid (STRp(cigar))
{
    uint32_t i=0;
    while (i < cigar_len) {

        uint32_t num_digits=0;
        for (; i < cigar_len && IS_DIGIT(cigar[i]) ; i++) num_digits++;

        if (!num_digits) return false;

        if (i == cigar_len || !cigar_valid_op[(int)cigar[i++]])
            return false;
    }
    return true;
}

// reverses a CIGAR, eg "40S111M"->"111M40S". returns false if not a valid CIGAR string.
bool sam_cigar_reverse (char *dst, STRp(cigar))
{
    if (cigar_len==1 && *cigar=='*') {
        *dst = '*';
        return true;
    }

    rom c = &cigar[cigar_len-1];
    while (c >= cigar) {

        char cigar_op = *(c--);
        if (!cigar_valid_op[(int)cigar_op]) return false;

        uint32_t num_digits=0;
        for (; c >= cigar && IS_DIGIT(*c) ; c--) num_digits++;
        
        if (!num_digits) return false;

        memcpy (dst, c+1, num_digits);
        dst[num_digits] = cigar_op;
        
        dst += num_digits + 1;
    }

    return true;
}

//---------------------------------------------------------------==------------------------------------
// Squanking - removing the longest number from the CIGAR string if it can be recovered from elsewhere:
// - for SA/XA/OA CIGARs - from the seq_len implied by the primary CIGARs
// - for short read data - from the segconf.seq_len
//-----------------------------------------------------------------------------------------------------
typedef enum { SEQ_LEN_FROM_MAIN='0', SEQ_LEN_FROM_SEGCONF='1'} SeqLenSource; // these values are part of the file format

static bool squank_seg (VBlockSAMP vb, ContextP ctx, STRp(cigar), uint32_t only_if_seq_len/*0=always*/,
                        SeqLenSource seq_len_source, uint32_t add_bytes)
{
    int32_t n=-1, max_n=-1; // -1 to be careful: n=0 is not expected, but IS a valid number by the SAM/BAM spec
    uint32_t start_n=0, segment1_len=0, start_segment2=0, seq_len=0;
    
    for (uint32_t i=0; i < cigar_len; i++) {
        char c = cigar[i];
        if (IS_DIGIT(c)) {
            if (n==-1) { // new number
                n = 0;
                start_n = i;
            }
            n = n*10 + c-'0';
        }

        else { // op
            if (M || I || S || E || X) {
                seq_len += n;

                if (n > max_n) {
                    max_n          = n;
                    segment1_len   = start_n;
                    start_segment2 = i; 
                }
            }
            n = -1;
        }
    }

    // case: we can squank
    if (!only_if_seq_len || only_if_seq_len == seq_len) {
        buf_alloc (vb, &ctx->local, cigar_len+1, 0, char, CTX_GROWTH, "local");
        buf_add (&ctx->local, cigar, segment1_len);
        BNXTc (ctx->local) = 0;
 
        // if squanking cuts out an S value, we can also remove the 'S' op as we can deduce it
        if (cigar[start_segment2] == 'S') start_segment2++; 

        if (start_segment2 < cigar_len)
            buf_add (&ctx->local, &cigar[start_segment2], cigar_len - start_segment2);
        
        BNXTc (ctx->local) = 0;

        ctx->local_num_words++;

        if (seq_len_source == SEQ_LEN_FROM_MAIN) // SA/XA/OA CIGAR (seq_len is compared MAIN's)
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_SQUANK, SEQ_LEN_FROM_MAIN }, 3, ctx, add_bytes);
        
        else // MAIN CIGAR - go through CIGAR special first
            seg_by_ctx (VB, (char[]){ SNIP_SPECIAL, SAM_SPECIAL_CIGAR, SNIP_LOOKUP }, 3, ctx, add_bytes);

        return true; // success
    }

    // case: add to local as is, with a LOOKUP
    else 
        return false; // cannot squank
}

// Lookup squanted CIGAR from local
SPECIAL_RECONSTRUCTOR_DT (sam_piz_special_SQUANK)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    bool is_main_cigar = (snip[0] == SEQ_LEN_FROM_SEGCONF);

    int32_t seq_len = is_main_cigar ? segconf.sam_seq_len // MAIN vs the "standard" seq_len (useful for short reads in which most reads are the same length)
                                    : vb->seq_len + VB_SAM->hard_clip[0] + VB_SAM->hard_clip[1]; // SA/OA/XA vs MAIN: hard-clips in the MAIN CIGAR are counted as well

    LOAD_SNIP_FROM_LOCAL (ctx); // segment1 of squank
    STR(segment1);
    STRset (segment1, snip);    // copy since LOAD_SNIP_FROM_LOCAL is hard-coded to "snip"

    LOAD_SNIP_FROM_LOCAL (ctx); // segment2 of squank

    if (!is_main_cigar && !reconstruct) goto done; // nothing more to do 

    int32_t segment1_seq_len = sam_cigar_get_seq_len (STRa(segment1));
    int32_t segment2_seq_len = sam_cigar_get_seq_len (STRa(snip)); 
    int32_t missing_len = seq_len - segment1_seq_len - segment2_seq_len;
    ASSPIZ (missing_len >= 0, "Expecting missing_len=%d > 0. seq_len=%d segment1_seq_len=%d segment2_seq_len=%d",
            missing_len, seq_len, segment1_seq_len, segment2_seq_len);
            
    // reconstruct always if coming from MAIN - it is needed for sam_cigar_analyze even if reconstruct = false
    Buffer *buf = &vb->txt_data;

    if (is_main_cigar) {
        ASSERTNOTINUSE (vb->scratch);
        buf_alloc (vb, &vb->scratch, 0, segment1_seq_len + segment2_seq_len + str_int_len (missing_len), char, 2, "scratch");
        buf = &vb->scratch;
    }
                
    if (segment1_len) 
        buf_add (buf, segment1, segment1_len);
    
    buf_add_int_as_text (buf, missing_len);

    if (!snip_len || IS_DIGIT(snip[0])) 
        BNXT (char, *buf) = 'S'; // reconstruct removed S - see squank_seg

    if (snip_len) 
        buf_add (buf, snip, snip_len);

done:
    return false; // no new value
}

//---------
// SEG
//---------

// used for XA, OA, SA, and also CIGAR field in PRIM VBs
bool sam_seg_0A_cigar_cb (VBlockP vb, ContextP ctx, STRp (cigar), uint32_t repeat)
{
    // complicated CIGARs are better off in local - anything more than eg 112M39S 
    // note: we set no_stons=true in sam_seg_initialize_0X so we can use local for this rather than singletons
    if (repeat != (uint32_t)-1 && // note: -1==prim cigar - cannot squank as it is used to determine seq_len
        cigar_len > MAX_CIGAR_LEN_IN_DICT && 
        squank_seg (VB_SAM, ctx, STRa(cigar), DATA_LINE(vb->line_i)->SEQ.len, SEQ_LEN_FROM_MAIN, cigar_len))
        {} // squank succeeded - nothing to do

    else if (cigar_len > MAX_CIGAR_LEN_IN_DICT)
        seg_add_to_local_text (vb, ctx, STRa(cigar), true, cigar_len);
 
    // short CIGAR
    else 
        seg_by_ctx (vb, STRa(cigar), ctx, cigar_len);

    return true;
}

static void sam_cigar_seg_prim_cigar (VBlockSAMP vb, STRp(textual_cigar))
{
    ContextP sa_cigar_ctx = CTX(OPTION_SA_CIGAR);

    sam_seg_0A_cigar_cb (VB, sa_cigar_ctx, STRa(textual_cigar), (uint32_t)-1 /*-1=prim cigar*/);
    sa_cigar_ctx->txt_len      -= textual_cigar_len; // remove "add_bytes" - already account for in SAM_CIGAR
    sa_cigar_ctx->counts.count += textual_cigar_len; // count CIGAR field contribution to OPTION_SA_CIGAR, so sam_stats_reallocate can allocate the z_data between CIGAR and SA:Z
}

void sam_cigar_seg_textual (VBlockSAMP vb, ZipDataLineSAM *dl, uint32_t last_cigar_len, STRp(seq_data), STRp(qual_data))
{
    START_TIMER
    
    ContextP ctx = CTX(SAM_CIGAR);
    bool qual_is_available = (qual_data_len != 1 || *qual_data != '*');
    bool seq_is_available  = (seq_data_len  != 1 || *seq_data  != '*');

    ASSSEG (!(seq_is_available && *seq_data=='*'), seq_data, "seq_data=%.*s (seq_len=%u), but expecting a missing seq to be \"*\" only (1 character)", 
            seq_data_len, seq_data, seq_data_len);

    char cigar_snip[last_cigar_len + 50];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    uint32_t cigar_snip_len=2;

    // case: SEQ is "*" - we add a '-' to the CIGAR
    if (!seq_is_available) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->SEQ.len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->SEQ.len) { // CIGAR is not available
        ASSSEG (!seq_data_len || !qual_is_available || seq_data_len==dl->QUAL.len, seq_data,
                "Bad line: SEQ length is %u, QUAL length is %u, unexpectedly differ. SEQ=%.*s QUAL=%.*s", 
                seq_data_len, dl->QUAL.len, seq_data_len, seq_data, dl->QUAL.len, qual_data);    

        dl->SEQ.len = MAX_(seq_data_len, dl->QUAL.len); // one or both might be not available and hence =1

        cigar_snip_len += str_int (dl->SEQ.len, &cigar_snip[cigar_snip_len]);
    } 
    else { // CIGAR is available - just check the seq and qual lengths
        ASSSEG (!seq_is_available || seq_data_len == dl->SEQ.len, seq_data,
                "Bad line: according to CIGAR, expecting SEQ length to be %u but it is %u. SEQ=%.*s", 
                dl->SEQ.len, seq_data_len, seq_data_len, seq_data);

        ASSSEG (!qual_is_available || qual_data_len == dl->SEQ.len, qual_data,
                "Bad line: according to CIGAR, expecting QUAL length to be %u but it is %u. QUAL=%.*s", 
                dl->SEQ.len, dl->QUAL.len, dl->QUAL.len, qual_data);    
    }

    // store the CIGAR in DataLine for use by a buddy MC:Z and SA:Z
    dl->CIGAR = (TxtWord){ .index = BNUMtxt (vb->last_cigar), .len = last_cigar_len }; // in SAM (but not BAM) vb->last_cigar points into txt_data

    // case: we buddy non-trival CIGARs with MC:Z. We don't buddy eg "151M" bc this will add rather than reduce entropy.
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1

    // case: DEPN or PRIM line.
    // Note: in DEPN, cigar already verified in sam_sa_seg_depn_find_sagroup to be the same as in SA alignment
    if (sam_seg_has_SA_Group(vb)) {

        sam_seg_against_sa_group_bool (vb, ctx, vb->soft_clip[0] || vb->soft_clip[1], last_cigar_len+1); // +1 for \t

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (sam_is_prim_vb) 
            sam_cigar_seg_prim_cigar (vb, vb->last_cigar, last_cigar_len);
    }

    // case: copy from buddy
    else if (last_cigar_len > 4 && vb->buddy_line_i != -1 && segconf.has[OPTION_MC_Z] && !segconf.running && 
        cigar_snip_len == 2 && // we don't buddy if CIGAR or SEQ are "*"
        buddy_dl->MC.len == last_cigar_len && 
        !memcmp (vb->last_cigar, Bc (vb->txt_data, buddy_dl->MC.index), last_cigar_len)) {

        cigar_snip[cigar_snip_len++] = COPY_BUDDY; // always at cigar_snip[2]
        seg_by_did_i (VB, STRa(cigar_snip), SAM_CIGAR, last_cigar_len+1); // +1 for \t
    }

    // case: long CIGAR and SEQ and CIGAR are not missing, with the "standard" sam_seq_len (normally only works for short reads)
    else if (last_cigar_len > MAX_CIGAR_LEN_IN_DICT && cigar_snip_len == 2 && dl->SEQ.len == segconf.sam_seq_len)
        squank_seg (vb, ctx, vb->last_cigar, last_cigar_len, 0/*always*/, SEQ_LEN_FROM_SEGCONF, last_cigar_len+1); 

    // case: long CIGAR and SEQ or CIGAR are missing or short CIGAR
    else { 
        memcpy (&cigar_snip[cigar_snip_len], vb->last_cigar, last_cigar_len);
        
        cigar_snip_len += last_cigar_len;

        if (last_cigar_len > MAX_CIGAR_LEN_IN_DICT) {
            seg_add_to_local_text (VB, ctx, STRa(cigar_snip), true, last_cigar_len+1);
        }
        else 
            seg_by_ctx (VB, STRa(cigar_snip), ctx, last_cigar_len+1); // +1 for \t
    }
    
    if (segconf.running) {
        segconf.sam_cigar_len += last_cigar_len;
        segconf.sam_seq_len   += seq_data_len;
    }

    COPY_TIMER(sam_cigar_seg);
}


void sam_cigar_seg_binary (VBlockSAMP vb, ZipDataLineSAM *dl, uint32_t l_seq, rom cigar, uint32_t n_cigar_op)
{
    START_TIMER

    ContextP ctx = CTX(SAM_CIGAR);
    char cigar_snip[vb->textual_cigar.len + 20];
    cigar_snip[0] = SNIP_SPECIAL;
    cigar_snip[1] = SAM_SPECIAL_CIGAR;
    uint32_t cigar_snip_len=2;

    // case: we have no sequence - we add a '-' to the CIGAR
    if (!l_seq) cigar_snip[cigar_snip_len++] = '-';

    // case: CIGAR is "*" - we get the dl->SEQ.len directly from SEQ or QUAL, and add the length to CIGAR eg "151*"
    if (!dl->SEQ.len) { // CIGAR is not available
        dl->SEQ.len = l_seq;
        cigar_snip_len += str_int (MAX_(dl->SEQ.len, 1), &cigar_snip[cigar_snip_len]); // if seq_len=0, then we add "1" because we have the * in seq and qual (and consistency with sam_cigar_seg_textual)
    }
    
    // case: we buddy non-trival CIGARs with MC:Z. We don't buddy eg "151M" bc this will add rather than reduce entropy.
    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    uint32_t add_bytes = n_cigar_op * sizeof (uint32_t) /* cigar */ + sizeof (uint16_t) /* n_cigar_op */;

    // case: DEPN or PRIM line.
    // Note: in DEPN, cigar already verified in sam_sa_seg_depn_find_sagroup to be as in SA alignment
    if (sam_seg_has_SA_Group(vb)) {

        sam_seg_against_sa_group_bool (vb, ctx, (vb->soft_clip[0] || vb->soft_clip[1]), add_bytes);

        // in PRIM, we also seg it as the first SA alignment (used for PIZ to load alignments to memory, not used for reconstructing SA)
        if (sam_is_prim_vb) 
            sam_cigar_seg_prim_cigar (vb, STRb(vb->textual_cigar));
    }

    // case: copy from buddy
    else if (vb->textual_cigar.len > 4 && vb->buddy_line_i != -1 && segconf.has[OPTION_MC_Z] && !segconf.running && 
        buddy_dl->MC.len == vb->textual_cigar.len && 
        cigar_snip_len == 2 && // we don't buddy if CIGAR or SEQ are "*"
        !memcmp (vb->textual_cigar.data, Bc (vb->txt_data, buddy_dl->MC.index), vb->textual_cigar.len)) {

        cigar_snip[cigar_snip_len++] = COPY_BUDDY; // always at cigar_snip[2]
        seg_by_ctx (VB, STRa(cigar_snip), ctx, add_bytes);
    }

    // case: long CIGAR and SEQ and CIGAR are not missing, with the "standard" sam_seq_len (normally only works for short reads)
    else if (vb->textual_cigar.len > MAX_CIGAR_LEN_IN_DICT && cigar_snip_len == 2 && dl->SEQ.len == segconf.sam_seq_len)
        squank_seg (vb, ctx, STRb(vb->textual_cigar), 0/*always*/, SEQ_LEN_FROM_SEGCONF, add_bytes); 

    // case: long CIGAR and SEQ or CIGAR are missing or short CIGAR
    else {
        memcpy (&cigar_snip[cigar_snip_len], vb->textual_cigar.data, vb->textual_cigar.len);
        cigar_snip_len += vb->textual_cigar.len;

        if (vb->textual_cigar.len > MAX_CIGAR_LEN_IN_DICT) 
            seg_add_to_local_text (VB, ctx, STRa(cigar_snip), true, add_bytes);
        else 
            seg_by_ctx (VB, STRa(cigar_snip), ctx, add_bytes);
    }

    // store a copy of the CIGAR in buddy_textual_cigars for use by a buddy MC:Z
    if (segconf.has[OPTION_MC_Z] && !segconf.running) {
        dl->CIGAR =(TxtWord){ .index = vb->buddy_textual_cigars.len, .len = vb->textual_cigar.len }; // in BAM dl->CIGAR points into buddy_textual_cigars
        buf_add_buf (VB, &vb->buddy_textual_cigars, &vb->textual_cigar, char, "buddy_textual_cigars");
    }

    if (segconf.running) {
        segconf.sam_cigar_len += cigar_snip_len;
        segconf.sam_seq_len   += dl->SEQ.len;
    }

    COPY_TIMER(sam_cigar_seg);
}

uint32_t sam_cigar_get_MC_ref_consumed (STRp(mc))
{
    // get ref_and_seq_consumed
    uint32_t n=0;
    uint32_t ref_and_seq_consumed=0;
    for (uint32_t i=0; i < mc_len; i++) {

        char c = mc[i];
        char lookup = cigar_lookup_sam[(uint8_t)c];
        if (!lookup) return 0; // invalid CIGAR - unrecognized character

        lookup &= 0x0f; // remove validity bit

        if (lookup == CIGAR_DIGIT) 
            n = n*10 + (c - '0');
        
        else {
            if (!n) return 0; // invalid CIGAR - no number before op

            if ((lookup & CIGAR_CONSUMES_REFERENCE)) 
                ref_and_seq_consumed += n;

            n=0;
        }
    }
    return ref_and_seq_consumed;
}

// MC:Z "CIGAR string for mate/next segment" (https://samtools.github.io/hts-specs/SAMtags.pdf)
void sam_cigar_seg_MC (VBlockSAMP vb, ZipDataLineSAM *dl, STRp(mc), uint32_t add_bytes)
{
    segconf_set_has (OPTION_MC_Z);

    ZipDataLineSAM *buddy_dl = DATA_LINE (vb->buddy_line_i); // an invalid pointer if buddy_line_i is -1
    Buffer *buddy_cigar_buf = IS_BAM_ZIP ? &vb->buddy_textual_cigars : &vb->txt_data; // buddy_dl->MC points into this buffer

    // we buddy non-trival CIGARs. We don't buddy eg "151M" bc this will add rather than reduce entropy.
    if (mc_len > 4 && !segconf.running && vb->buddy_line_i != -1 && 
        mc[0] != '*' &&
        buddy_dl->CIGAR.len == mc_len && !memcmp (mc, Bc (*buddy_cigar_buf, buddy_dl->CIGAR.index), mc_len)) 

        seg_by_did_i (VB, STRa(MC_buddy_snip), OPTION_MC_Z, add_bytes); // copy MC from earlier-line buddy CIGAR
    
    else if (mc_len > 7) 
        seg_add_to_local_text (VB, CTX(OPTION_MC_Z), STRa(mc), true, add_bytes);
    else
        seg_by_did_i (VB, STRa(mc), OPTION_MC_Z, add_bytes);    

    dl->MC = TXTWORD(mc); 

    ctx_set_last_value (VB, CTX(OPTION_MC_Z), (ValueType){ .i = sam_cigar_get_MC_ref_consumed (STRa(mc)) } );
}

//---------
// PIZ
//---------

// CIGAR - calculate vb->seq_len from the CIGAR string, and if original CIGAR was "*" - recover it
SPECIAL_RECONSTRUCTOR_DT (sam_cigar_special_CIGAR)
{
    VBlockSAMP vb = (VBlockSAMP)vb_;
    
    const uint16_t sam_flag = vb->last_int(SAM_FLAG);

    // case: we copy the snip from buddy MC:Z
    if (snip[0] == COPY_BUDDY) 
        reconstruct_from_buddy_get_textual_snip (VB, CTX (OPTION_MC_Z), pSTRa(snip));

    if (snip[0] == SNIP_LOOKUP) { // squank into vb->scratch
        sam_piz_special_SQUANK (VB, ctx, (char[]){ SEQ_LEN_FROM_SEGCONF }, 1, new_value, reconstruct); 
        snip     = vb->scratch.data;
        snip_len = vb->scratch.len;
    }

    // calculate seq_len (= l_seq, unless l_seq=0), ref_consumed and (if bam) vb->textual_cigar and vb->binary_cigar
    sam_cigar_analyze (vb, STRa(snip), false, &vb->seq_len); 

    if ((flag.out_dt == DT_SAM || (flag.out_dt == DT_FASTQ && flag.extended_translation)) 
    &&  reconstruct) {
        if (snip[snip_len-1] == '*') // eg "151*" - zip added the "151" to indicate seq_len - we don't reconstruct it, just the '*'
            RECONSTRUCT1 ('*');
        
        else if (snip[0] == '-') // eg "-151M" or "-151*" - zip added the "-" to indicate a '*' SEQ field - we don't reconstruct it
            RECONSTRUCT (snip + 1, snip_len - 1);

        else 
            RECONSTRUCT (snip, snip_len);    
    }

    // BAM - output vb->binary_cigar generated in sam_cigar_analyze
    else if (flag.out_dt == DT_BAM) {
        // now we have the info needed to reconstruct bin, l_read_name, n_cigar_op and l_seq
        BAMAlignmentFixed *alignment = (BAMAlignmentFixed *)Bc (vb->txt_data, vb->line_start);
        alignment->l_read_name = BAFTc (vb->txt_data) - alignment->read_name;
        alignment->n_cigar_op  = LTEN16 (vb->binary_cigar.len);
        alignment->l_seq       = (snip[0] == '-') ? 0 : LTEN32 (vb->seq_len);

        LTEN_u32_buf (&vb->binary_cigar, NULL);
        RECONSTRUCT (vb->binary_cigar.data, vb->binary_cigar.len * sizeof (BamCigarOp));
        LTEN_u32_buf (&vb->binary_cigar, NULL); // restore

        // if BIN is SAM_SPECIAL_BIN, inst.semaphone is set by bam_piz_special_BIN - a signal to us to calculate
        ContextP sam_bam_bin_ctx = CTX(SAM_BAM_BIN);
        if (sam_bam_bin_ctx->semaphore) {
            sam_bam_bin_ctx->semaphore = false;

            PosType pos = CTX(SAM_POS)->last_value.i;
            bool segment_unmapped = (sam_flag & SAM_FLAG_UNMAPPED);
            PosType last_pos = segment_unmapped ? pos : (pos + vb->ref_consumed - 1);
            
            uint16_t bin = bam_reg2bin (pos, last_pos); // zero-based, half-closed half-open [start,end)
            alignment->bin = LTEN16 (bin); // override the -1 previously set by the translator
        }
    }
    
    else if (flag.out_dt == DT_FASTQ) {
        // only analyze, but don't reconstruct CIGAR in FASTQ
    }

    buf_free (vb->scratch);
    return false; // no new value
}   

// copy from buddy CIGAR to MC. If reconstructing to BAM, we convert the binary CIGAR to textual MC:Z.*
SPECIAL_RECONSTRUCTOR (sam_piz_special_COPY_BUDDY_MC)
{
    if (!reconstruct) return false;

    VBlockSAMP sam_vb = (VBlockSAMP)vb;

    // fall back to normal COPY_BUDDY in case of SAM
    if (flag.out_dt != DT_BAM) return reconstruct_from_buddy (vb, ctx, STRa(snip), reconstruct, new_value); // SAM or FASTQ output

    // get CIGAR field value previously reconstructed in BAM binary format
    STR(bam_cigar);
    reconstruct_from_buddy_get_textual_snip (vb, CTX (SAM_CIGAR), pSTRa(bam_cigar));
    
    // convert binary CIGAR to textual MC:Z
    uint32_t n_cigar_op = bam_cigar_len / sizeof (uint32_t);
    sam_cigar_binary_to_textual (sam_vb, n_cigar_op, (uint32_t*)bam_cigar, &vb->txt_data);

    return false; // no new value 
}

// invoked from TOP2FQ (but not TOP2FQEX, bc it reconstructs AUX) to consume MC if it exists in this line, in case this line is
// a buddy line of a future line in which case this MC will be copied to the future line's CIGAR 
SPECIAL_RECONSTRUCTOR (sam_piz_special_CONSUME_MC_Z)
{
    ContextP mc_ctx = CTX(OPTION_MC_Z);
    if (!mc_ctx->flags.store_per_line) goto done; // MC is not buddied
    
    ContextP opt_ctx = CTX(SAM_AUX);
    WordIndex opt_word_index = WORD_INDEX_NONE;

    // get AUX container
    snip_len=0;
    if (opt_ctx->b250.len ||
        (!opt_ctx->b250.len && !opt_ctx->local.len && opt_ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        opt_word_index = LOAD_SNIP(opt_ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (snip_len==1 && *snip == SNIP_LOOKUP)
            snip_len=0;
    }
    
    // case: a singleton (SNIP_LOOKUP) or all data is local
    if (!snip_len) 
        LOAD_SNIP_FROM_LOCAL (opt_ctx);

    if (!snip_len || snip[0] != SNIP_CONTAINER) goto done; // not a container

    ContainerP opt_con = container_retrieve (vb, opt_ctx, opt_word_index, snip+1, snip_len-1, 0, 0);

    // check if this line has an optional tag MC:Z
    bool found = false;
    for (uint32_t item_i=0; item_i < con_nitems (*opt_con); item_i++)
        if (opt_con->items[item_i].dict_id.num == _OPTION_MC_Z) {
            found = true;
            break;
        }
    if (!found) goto done; // AUX has no MC:Z tag

    // store MC:Z in history to be used by future line CIGAR. Note: since MC:Z is not reconstructed, ctx->history will point
    // to either dict or local (depending on where this snip originates) rather than the normal txt_data.
    snip_len = 0;

    // case: MC:Z is in dict (refered to from a b250)
    if (mc_ctx->b250.len ||
        (!mc_ctx->b250.len && !mc_ctx->local.len && mc_ctx->dict.len)) {  // all_the_same case - no b250 or local, but have dict      
        LOAD_SNIP(mc_ctx->did_i); // note: if we have no b250, local but have dict, this will be word_index=0 (see ctx_get_next_snip)

        if (snip_len==1 && *snip == SNIP_LOOKUP)
            snip_len=0;

        else 
            *B(HistoryWord, mc_ctx->history, vb->line_i - vb->first_line) = 
                (HistoryWord){ .char_index = BNUM (mc_ctx->dict, snip), .snip_len = snip_len, .lookup = LookupDict };
    }
    
    // case: MC:Z is in local
    if (!snip_len) {
        uint32_t char_index = LOAD_SNIP_FROM_LOCAL (mc_ctx);
        *B(HistoryWord, mc_ctx->history, vb->line_i - vb->first_line) = 
            (HistoryWord){ .char_index = char_index, .snip_len = snip_len, .lookup = LookupLocal }; 
    }

done:
    return false; // no new value 
}

// called from sam_piz_special_pull_from_SAGROUP for reconstructing the main CIGAR field of a PRIM / DEPN line
void sam_reconstruct_main_cigar_from_SA_Group (VBlockSAMP vb, bool substitute_S, bool reconstruct)
{
    // we generate the CIGAR in vb->scratch. sam_cigar_special_CIGAR will reconstruct it (possibly binary) in txt_data. 
    ASSERTNOTINUSE (vb->scratch);
    const SAAlnType *a = vb->sa_aln;
    rom cigar_snip;
    uint32_t cigar_len;

    // case: cigar is stored in dict 
    if (a->cigar.piz.is_word) {
        ctx_get_snip_by_word_index_do (CTX(OPTION_SA_CIGAR), a->cigar.piz.index, &cigar_snip, &cigar_len, __FUNCLINE);
        buf_add_more (VB, &vb->scratch, cigar_snip, cigar_len, "scratch");
    }

    // case: cigar is stored in local  
    else {
        cigar_len = a->cigar.piz.len_lo | (a->cigar.piz.len_hi << ALN_CIGAR_LEN_BITS_LO);
        buf_alloc (vb, &vb->scratch, 0, cigar_len,  char, 0, "scratch");

        // case: compressed
        if (a->cigar.piz.comp_len) {      
            uint8_t *comp = B8(z_file->sa_cigars, a->cigar.piz.index);
            uint32_t uncomp_len = cigar_len;
            void *success = rans_uncompress_to_4x16 (VB, comp, a->cigar.piz.comp_len,
                                                     B1ST(uint8_t, vb->scratch), &uncomp_len); 
            ASSPIZ (success && uncomp_len == cigar_len, "rans_uncompress_to_4x16 failed to decompress an SA Aln CIGAR data: vb_i=%u line_i=%"PRIu64" grp_i=%u aln_i=%"PRIu64" success=%u comp_len=%u uncomp_len=%u expected_uncomp_len=%u cigar_index=%"PRIu64" comp[10]=%s",
                    vb->vblock_i, vb->line_i, ZGRP_I(vb->sa_grp), ZALN_I(a), !!success, (uint32_t)a->cigar.piz.comp_len, uncomp_len, cigar_len, (uint64_t)a->cigar.piz.index, str_hex10 (comp, a->cigar.piz.comp_len).s);
        }

        // case: not compressed
        else 
            memcpy (B1ST(uint8_t, vb->scratch), Bc(z_file->sa_cigars, a->cigar.piz.index), cigar_len);
    }

    char *cigar = B1STc (vb->scratch);

    // case: we need to replace soft-clipping (S) with hard-clipping (H)
    if (substitute_S) sam_cigar_S_to_H (STRa(cigar));

    sam_cigar_special_CIGAR (VB, CTX(SAM_CIGAR), STRa(cigar), NULL, reconstruct);
// printf ("xxx cigar=%.*s\n", cigar_len, cigar);
    buf_free (vb->scratch);
}

// called from sam_sa_reconstruct_SA_from_SA_Group for reconstructing a CIGAR in an SA:Z field of a PRIM/DEPN line
void sam_reconstruct_SA_cigar_from_SA_Group (VBlockSAMP vb, SAAlnType *a)
{
    if (a->cigar.piz.is_word) {
        STR(cigarS);
        ctx_get_snip_by_word_index (CTX(OPTION_SA_CIGAR), a->cigar.piz.index, cigarS);
        RECONSTRUCT_SEP (cigarS, cigarS_len, ',');
    }

    else {
        uint32_t cigar_len = a->cigar.piz.len_lo | (a->cigar.piz.len_hi << ALN_CIGAR_LEN_BITS_LO);

        if (a->cigar.piz.comp_len) { // compressed
            uint32_t uncomp_len = cigar_len;
    
            void *success = rans_uncompress_to_4x16 (VB, B8(z_file->sa_cigars, a->cigar.piz.index), a->cigar.piz.comp_len,
                                                     BAFT(uint8_t, vb->txt_data), &uncomp_len); 
            ASSPIZ (success && uncomp_len == cigar_len, "rans_uncompress_to_4x16 failed to decompress an SA Aln CIGAR data: vb_i=%u line_i=%"PRIu64" grp_i=%u aln_i=%"PRIu64" success=%u comp_len=%u uncomp_len=%u expected_uncomp_len=%u cigar_index=%"PRIu64,
                    vb->vblock_i, vb->line_i, ZGRP_I(vb->sa_grp), ZALN_I(a), !!success, (uint32_t)a->cigar.piz.comp_len, uncomp_len, cigar_len, (uint64_t)a->cigar.piz.index);
// printf("xxx SAcigar comp=%.*s\n", cigar_len, BAFTtxt);
            vb->txt_data.len += cigar_len;
        }

        else { // not compressed
// printf("xxx SAcigar noncomp=%.*s\n", cigar_len, B8(z_file->sa_cigars, a->cigar));
            RECONSTRUCT (B8(z_file->sa_cigars, a->cigar.piz.index), cigar_len);
}
        RECONSTRUCT1 (',');
    }
}

// PIZ: main thread (not thread-safe): called from sam_show_sa_one_grp for getting first few characters of alignment cigar
rom sam_piz_display_aln_cigar (const SAAlnType *a)
{
    static char cigar[SA_CIGAR_DISPLAY_LEN+1];
    memset (cigar, 0, sizeof(cigar));

    if (a->cigar.piz.is_word) {
        ContextP ctx = ZCTX(OPTION_SA_CIGAR);

        if (a->cigar.piz.index < ctx->word_list.len) {
            STR(cigarS);
            ctx_get_snip_by_word_index (ctx, a->cigar.piz.index, cigarS);
            memcpy (cigar, cigarS, MIN_(cigarS_len, SA_CIGAR_DISPLAY_LEN));
        }
        else
            strcpy (cigar, "BAD_WORD");
    }

    else {
        uint32_t cigar_len = ALN_CIGAR_LEN(a);
        uint32_t uncomp_len = MIN_(SA_CIGAR_DISPLAY_LEN, cigar_len); // possibly shorter than original cigar

        if (a->cigar.piz.comp_len) { // compressed
            void *success = rans_uncompress_to_4x16 (evb, B8(z_file->sa_cigars, a->cigar.piz.index), a->cigar.piz.comp_len,
                                                    (uint8_t *)cigar, &uncomp_len); 
            if (success && uncomp_len) cigar[uncomp_len] = '\0';
        }

        else // not compressed
            memcpy (cigar, B8(z_file->sa_cigars, a->cigar.piz.index), uncomp_len);
    }

    return cigar;
}

//---------------------------------------------------------------------------------------------------
// CIGAR signature
// Note: the signature is in-memory and is not written to the genozip file, so can be changed at will
//---------------------------------------------------------------------------------------------------

CigarSignature cigar_sign (STRp(cigar))
{
    START_TIMER;
    CigarSignature sig;

    // case: cigar is not longer than the signature - the cigar IS the signature
    if (cigar_len <= CIGAR_SIG_LEN) {
        memcpy (sig.bytes, cigar, cigar_len);
        memset (sig.bytes + cigar_len, 0, CIGAR_SIG_LEN - cigar_len);
    }

    // case: long cigar - use MD5 (note: I tried using Adler32 and got contention in real data)
    else {
        Digest digest = md5_do (STRa(cigar));
        memcpy (sig.bytes, digest.bytes, CIGAR_SIG_LEN);
    }
COPY_TIMER_VB(evb,tmp1);
    return sig;
}

bool cigar_is_same_signature (CigarSignature sig1, CigarSignature sig2) 
{
    return !memcmp (sig1.bytes, sig2.bytes, CIGAR_SIG_LEN);
}

DisCigarSig cigar_display_signature (CigarSignature sig)
{
    DisCigarSig dis;
    
    str_to_hex (sig.bytes, CIGAR_SIG_LEN, dis.s);

    return dis;
}


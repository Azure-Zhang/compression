// ------------------------------------------------------------------
//   seg_sam.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "move_to_front.h"
#include "header.h"
#include "file.h"
#include "random_access.h"
#include "endianness.h"
#include "strings.h"

#define DATA_LINE(vb,i) (&((ZipDataLineSAM *)((vb)->data_lines))[(i)])

// called from seg_all_data_lines
void seg_sam_initialize (VBlock *vb_)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;

    buf_alloc (vb, &vb->md_data, 12 * vb->num_lines, 1, "md_data", vb->vblock_i);
    buf_alloc (vb, &vb->random_pos_data, vb->num_lines * sizeof (uint32_t), 1, "random_pos_data", vb->vblock_i);    
}             

// if this row's tlen is the negative of the previous row - store "*" instead of the tlen - thereby cutting
// the number of dictionary words by half and improving compression
static inline void seg_sam_tlen_field (VBlockSAM *vb, const char *tlen, unsigned tlen_len, unsigned vb_line_i)
{
    ASSERT (tlen_len, "%s: empty TLEN field found in vb_i=%u vb_line_i=%u", global_cmd, vb->vblock_i, vb_line_i);

    bool tlen_is_positive = (tlen[0] != '-');

    // case: tlen is a negative number, with absolute value the same as positive vb->tlen_snip
    // e.g. tlen="-1000" ; previous line="1000"
    bool this_neg_prev_pos = vb->last_tlen_abs && 
                             (vb->last_tlen_abs_len == tlen_len-1) && 
                             !tlen_is_positive && vb->last_tlen_is_positive &&
                             !memcmp (vb->last_tlen_abs, &tlen[1], tlen_len - 1);

    bool this_pos_prev_neg = vb->last_tlen_abs && 
                             (vb->last_tlen_abs_len == tlen_len) && 
                             tlen_is_positive && !vb->last_tlen_is_positive &&
                             !memcmp (vb->last_tlen_abs, tlen, tlen_len);

    // case: tlen is the negative of previous line - replace it with "*"
    if (this_neg_prev_pos || this_pos_prev_neg) {
        seg_one_field (vb, "*", 1, vb_line_i, SAM_TLEN);
        vb->txt_section_bytes[SEC_SAM_QNAME_B250] += tlen_len - 1; // seg_one_field only added 2 (1 for \t), we should add the rest
    }

    // case: tlen is not negative of previous line - add as is
    else 
        seg_one_field ((VBlockP)vb, tlen, tlen_len, vb_line_i, SAM_TLEN);

    vb->last_tlen_abs         = &tlen[!tlen_is_positive]; // store absolute value
    vb->last_tlen_is_positive = tlen_is_positive;
    vb->last_tlen_abs_len     = tlen_len - !tlen_is_positive;
}

// returns length of string ending with separator, or -1 if separator was not found
static inline int seg_sam_get_next_subitem (const char *str, int str_len, char separator)
{
    for (int i=0; i < str_len; i++) {
        if (str[i] == separator) return i;
        if (str[i] == ',' || str[i] == ';') return -1; // wrong separator encountered
    }
    return -1;
}

#define DO_SSF(ssf,sep) \
        ssf = &field[i]; \
        ssf##_len = seg_sam_get_next_subitem (&field[i], field_len-i, sep); \
        if (ssf##_len == -1) goto error; /* bad format */ \
        i += ssf##_len + 1; /* skip snip and separator */        

#define DEC_SSF(ssf) const char *ssf; \
                     int ssf##_len;

static inline uint32_t seg_sam_add_to_random_pos_data (VBlockSAM *vb, SectionType sec, 
                                                       const char *snip, unsigned snip_len, unsigned add_bytes,
                                                       const char *field_name)
{
    bool standard_random_pos_encoding = true;

    // it is eligable for standard encoding only if the length is within the range
    if (!snip_len || snip_len > 10) standard_random_pos_encoding = false; // more than 10 digits is for sure bigger than 4GB=32bits

    // a multi-digit number cannot have a leading zero
    if (snip_len > 1 && snip[0] == '0') standard_random_pos_encoding = false;

    // it is eligable for standard encoding only if all the characters are digits
    uint64_t n64=0; // 64 bit just in case we go above 32 bit
    if (standard_random_pos_encoding)
        for (unsigned i=0; i < snip_len; i++) {
            if (!IS_DIGIT (snip[i])) {
                standard_random_pos_encoding = false;
                break;
            }
            n64 = n64 * 10 + (snip[i] - '0');
        }

    // it is eligable for standard encoding only if it fits in 32 bit (0xffffffff is reserved for a future escape, if needed)
    if (n64 > 0xfffffffe) standard_random_pos_encoding = false;

    ASSERT (standard_random_pos_encoding, "%s: Error in %s, vb_i=%u: Bad position data in field %s - expecting an integer between 0 and %u without leading zeros, but instead seeing \"%.*s\"",
            global_cmd, txt_name, vb->vblock_i, field_name, 0xfffffffe, snip_len, snip);

    buf_alloc_more (vb, &vb->random_pos_data, 1, 0, uint32_t, 2);

    // 32 bit number - this compresses better than textual numbers with LZMA
    uint32_t n32 = (uint32_t)n64;
    NEXTENT (uint32_t, vb->random_pos_data) = BGEN32 (n32);

    vb->txt_section_bytes[sec] += add_bytes;

    return n32;
}

static unsigned seg_sam_SA_or_OA_field (VBlockSAM *vb, const char *field, unsigned field_len, unsigned vb_line_i,
                                        const char *field_name)
{
    // OA and SA format is: (rname ,pos ,strand ,CIGAR ,mapQ ,NM ;)+ . in OA - NM is optional (but its , is not)
    // Example SA:Z:chr13,52863337,-,56S25M70S,0,0;chr6,145915118,+,97S24M30S,0,0;chr18,64524943,-,13S22M116S,0,0;chr7,56198174,-,20M131S,0,0;chr7,87594501,+,34S20M97S,0,0;chr4,12193416,+,58S19M74S,0,0;
    // See: https://samtools.github.io/hts-specs/SAMtags.pdf
    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(strand); DEC_SSF(cigar); DEC_SSF(mapq); DEC_SSF(nm); 
    int repeats = 0;

    for (unsigned i=0; i < field_len; repeats++) {
        DO_SSF (rname,  ','); 
        DO_SSF (pos,    ','); 
        DO_SSF (strand, ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (mapq,   ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (strand_len != 1 || (strand[0] != '+' && strand[0] != '-')) goto error; // invalid format

        // pos - add to the random pos data together with all other random pos data (originating from POS, PNEXT etc).
        seg_sam_add_to_random_pos_data (vb, SEC_SAM_RAND_POS_DATA, pos, pos_len, pos_len+1, field_name);

        // nm : we store in the same dictionary as the Optional subfield NM
        seg_one_subfield ((VBlockP)vb, nm, nm_len, vb_line_i, (DictIdType)dict_id_OPTION_NM, 
                          SEC_SAM_OPTNL_SF_B250, nm_len+1);

        // strand : we store in in the STRAND dictionary (used by SA, OA, XA)
        seg_one_subfield ((VBlockP)vb, strand, 1, vb_line_i, (DictIdType)dict_id_OPTION_STRAND, 
                          SEC_SAM_OPTNL_SF_B250, 2);

        // rname, cigar and mapq: we store in the same dictionary as the primary fields
        seg_one_field ((VBlockP)vb, rname, rname_len, vb_line_i, SAM_RNAME);
        seg_one_field ((VBlockP)vb, cigar, cigar_len, vb_line_i, SAM_CIGAR);
        seg_one_field (vb, mapq,  mapq_len, vb_line_i, SAM_MAPQ);
    }

    return repeats;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we return 0 and the subfield will be stored as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSERT (!repeats, "Invalid format in repeat #%u of field SA or OA, vb_i=%u vb_line_i=%u. data: %.*s",
            repeats+1, vb->vblock_i, vb_line_i, field_len, field);
    return 0;
}

static unsigned seg_sam_XA_field (VBlockSAM *vb, const char *field, unsigned field_len, unsigned vb_line_i)
{
    // XA format is: (chr,pos,CIGAR,NM;)*  pos starts with +- which is strand
    // Example XA:Z:chr9,-60942781,150M,0;chr9,-42212061,150M,0;chr9,-61218415,150M,0;chr9,+66963977,150M,1;
    // See: http://bio-bwa.sourceforge.net/bwa.shtml
    DEC_SSF(rname); DEC_SSF(pos); DEC_SSF(cigar); DEC_SSF(nm); 
    int repeats = 0;

    for (unsigned i=0; i < field_len; repeats++) {
        DO_SSF (rname,  ','); 
        DO_SSF (pos,    ','); 
        DO_SSF (cigar,  ','); 
        DO_SSF (nm,     ';'); 

        // sanity checks before adding to any dictionary
        if (pos_len < 2 || (pos[0] != '+' && pos[0] != '-')) goto error; // invalid format - expecting pos to begin with the strand

        // pos - add to the pos data together with all other pos data (POS, PNEXT etc),
        // strand - add to separate STRAND dictionary, to not adversely affect compression of POS.
        // there is no advantage to storing the strand together with pos as they are not correlated.
        seg_sam_add_to_random_pos_data (vb, SEC_SAM_RAND_POS_DATA, pos+1, pos_len-1, pos_len, "XA"); // one less for the strand, one more for the ,

        // strand
        seg_one_subfield ((VBlockP)vb, pos, 1, vb_line_i, (DictIdType)dict_id_OPTION_STRAND, 
                          SEC_SAM_OPTNL_SF_B250, 1);

        // nm : we store in the same dictionary as the Optional subfield NM
        seg_one_subfield ((VBlockP)vb, nm, nm_len, vb_line_i, (DictIdType)dict_id_OPTION_NM, 
                          SEC_SAM_OPTNL_SF_B250, nm_len+1);

        // rname and cigar: we store in the same dictionary as the primary fields
        seg_one_field (vb, rname, rname_len, vb_line_i, SAM_RNAME);
        seg_one_field (vb, cigar, cigar_len, vb_line_i, SAM_CIGAR);
    }

    return repeats;

error:
    // if the error occurred on on the first repeat - this file probably has a different
    // format - we return 0 and the subfield will be stored as a normal subfield
    // if it occurred on the 2nd+ subfield, after the 1st one was fine - we reject the file
    ASSERT (!repeats, "Invalid format in repeat #%u of field XA (0-based), vb_i=%u vb_line_i=%u. data: %.*s",
            repeats, vb->vblock_i, vb_line_i, field_len, field);
    return 0;
}

uint32_t seg_sam_get_seq_len_by_MD_field (const char *md_str, unsigned md_str_len, bool *is_numeric)
{
    uint32_t result=0, curr_num=0;

    for (unsigned i=0; i < md_str_len; i++) {   
        if (IS_DIGIT (md_str[i])) 
            curr_num = curr_num * 10 + (md_str[i] - '0');

        else {
            result += curr_num + 1; // number terminates here + one character
            curr_num = 0;
        }
    }

    result += curr_num; // in case the string ends with a number

    if (is_numeric) *is_numeric = (result == curr_num); // md is only digits

    return result;
}

// if the the case where sequence length as calculated from the MD is the same as that calculated
// from the CIGAR/SEQ/QUAL (note: this is required by the SAM spec but nevertheless genozip doesn't require it):
// MD is shortened to replace the last number with a *, since it can be calculated knowing the length. The result is that
// multiple MD values collapse to one, e.g. "MD:Z:119C30" and "MD:Z:119C31" both become "MD:Z:119C*" hence improving compression.
// In the case where the MD is simply a number "151" and drop it altogther and keep just an empty string.
static inline bool seg_sam_get_shortened_MD (const char *md_str, unsigned md_str_len, uint32_t seq_len
,                                            char *new_md_str, unsigned *new_md_str_len)
{
    bool is_numeric;
    uint32_t seq_len_by_md = seg_sam_get_seq_len_by_MD_field (md_str, md_str_len, &is_numeric);

    if (seq_len_by_md != seq_len) return false;  // MD string doesn't comply with SAM spec and is therefore not changed

    // case - MD is just a number eg "151" that is equal seq_len - we make it an empty string
    if (is_numeric) {
        *new_md_str_len = 0;
        return true;
    }
    
    // case - MD ends with a number eg "119C31" - we replace it with "119C*"
    else if (IS_DIGIT (md_str[md_str_len-1])) {
        memcpy (new_md_str, md_str, md_str_len);
        int i=md_str_len-1; for (; i>=0; i--)
            if (!IS_DIGIT (md_str[i])) break;
        
        new_md_str[i+1] = '*';
        *new_md_str_len = i+2;
        return true;
    }

    return false; // MD doesn't end with a number and is hence unchanged (this normally doesn't occur as the MD would finish with 0)
}

// AS and XS are values (at least as set by BWA) at most the seq_len, and AS is often equal to it. we modify
// it to be new_value=(value-seq_len) 
static inline void seg_sam_AS_field (VBlockSAM *vb, ZipDataLineSAM *dl, DictIdType dict_id, 
                                     const char *snip, unsigned snip_len, uint32_t vb_line_i)
{
    bool positive_delta = true;

    // verify that its a unsigned number
    for (unsigned i=0; i < snip_len; i++)
        if (!IS_DIGIT (snip[i])) positive_delta = false;

    int32_t as;
    if (positive_delta) {
        as = atoi (snip); // type i is signed 32 bit by SAM specification
        if (dl->seq_len < as) positive_delta=false;
    }

    // if possible, store the positive delta - preceded by a "*"
    if (positive_delta) {
        char new_snip[20];    
        unsigned delta_len;
        new_snip[0] = '*'; 
        str_uint (dl->seq_len-as, &new_snip[1], &delta_len);

        seg_one_subfield ((VBlock *)vb, new_snip, delta_len+1, vb_line_i, dict_id, SEC_SAM_OPTNL_SF_B250, snip_len + 1); 
    }

    // not possible - just store unmodified
    else
        seg_one_subfield ((VBlock *)vb, snip, snip_len, vb_line_i, dict_id, SEC_SAM_OPTNL_SF_B250, snip_len + 1); // +1 for \t
}

// process an optional subfield, that looks something like MX:Z:abcdefg. We use "MX" for the field name, and
// the data is abcdefg. The full name "MX:Z:" is stored as part of the OPTIONAL dictionary entry
static void seg_sam_optional_field (VBlockSAM *vb, ZipDataLineSAM *dl, const char *field, unsigned field_len, 
                                    char separator, unsigned vb_line_i)
{
    ASSERT (field_len, "Error in %s: line invalidly ends with a tab. vb_i=%u vb_line_i=%u", 
            txt_name, vb->vblock_i, vb_line_i);

    ASSERT (field_len >= 6 && field[2] == ':' && field[4] == ':', "Error in %s: invalid optional field format: %.*s",
            txt_name, field_len, field);

    DictIdType dict_id = dict_id_sam_optnl_sf (dict_id_make (field, 4));
    const char *value = &field[5]; // the "abcdefg" part of "MX:Z:abcdefg"
    unsigned value_len = field_len - 5;

    // some fields have "repeats" - multiple instances of the same format of data separated by a ;
    unsigned repeats = 0; 
    if (dict_id.num == dict_id_OPTION_SA || dict_id.num == dict_id_OPTION_OA)
        repeats = seg_sam_SA_or_OA_field (vb, value, value_len, vb_line_i, dict_id.num == dict_id_OPTION_SA ? "SA" : "OA");

    else if (dict_id.num == dict_id_OPTION_XA) 
        repeats = seg_sam_XA_field (vb, value, value_len, vb_line_i);

    // special subfield - in the subfield, store just the (char)-1 + the number of repeats.
    // the data itself was already stored in the subsubfields in seg_sam_*_field
    if (repeats) {
        char repeats_str[20];
        sprintf (repeats_str, "%c%u", -1, repeats); // (char)-1 to indicate repeats (invalid char per SAM specification, so it cannot appear organically)
        seg_one_subfield ((VBlock *)vb, repeats_str, strlen (repeats_str), vb_line_i, dict_id, SEC_SAM_OPTNL_SF_B250, 1 /* \t */); 
    }

    else { // non-repeating fields

        // fields containing CIGAR format data
        if (dict_id.num == dict_id_OPTION_MC || dict_id.num == dict_id_OPTION_OC) 
            seg_one_field (vb, value, value_len, vb_line_i, SAM_CIGAR);

        // MD's logical length is normally the same as seq_len, we use this to optimize it.
        // In the common case that it is just a number equal the seq_len, we replace it with an empty string.
        else if (dict_id.num == dict_id_OPTION_MD) {
            // if MD value can be derived from the seq_len, we don't need to store - store just an empty string
            char new_md[MAX_SAM_MD_LEN];
            unsigned new_md_len;
            bool md_is_changeable = (value_len <= MAX_SAM_MD_LEN);

            if (md_is_changeable) 
                md_is_changeable = seg_sam_get_shortened_MD (value, value_len, dl->seq_len, new_md, &new_md_len);

            seg_add_to_data_buf ((VBlockP)vb, &vb->md_data, SEC_SAM_MD_DATA, 
                                 md_is_changeable ? new_md : value, 
                                 md_is_changeable ? new_md_len : value_len,
                                 '\t', value_len+1);
        }

        // AS is a value (at least as set by BWA) at most the seq_len, and often equal to it. we modify
        // it to be new_AS=(AS-seq_len) 
        else if (dict_id.num == dict_id_OPTION_AS) 
            seg_sam_AS_field (vb, dl, dict_id, value, value_len, vb_line_i);
        
        // mc:i: (output of bamsormadup? - mc in small letters) appears to a pos value usually close to POS.
        // we encode as a delta.
        else if (dict_id.num == dict_id_OPTION_mc) {
            uint8_t mc_did_i = mtf_get_ctx_by_dict_id (vb->mtf_ctx, vb->dict_id_to_did_i_map, &vb->num_dict_ids, NULL, dict_id, SEC_SAM_OPTNL_SF_DICT)->did_i;

            seg_pos_field ((VBlockP)vb, vb->last_pos, mc_did_i, SEC_SAM_OPTNL_SF_B250, value, value_len, vb_line_i, "mc:i");
        }

        // E2 - SEQ data (note: E2 doesn't have a dictionary)
        else if (dict_id.num == dict_id_OPTION_E2) { 
            ASSERT (value_len == dl->seq_len, 
                    "Error in %s: Expecting E2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                    txt_name, dl->seq_len, value_len, value_len, value);

            dl->e2_data_start = value - vb->txt_data.data;
            dl->e2_data_len   = value_len;
            vb->txt_section_bytes[SEC_SEQ_DATA] += dl->seq_len + 1; // +1 for \t
        }
        // U2 - QUAL data (note: U2 doesn't have a dictionary)
        else if (dict_id.num == dict_id_OPTION_U2) {
            ASSERT (value_len == dl->seq_len, 
                    "Error in %s: Expecting U2 data to be of length %u as indicated by CIGAR, but it is %u. E2=%.*s",
                    txt_name, dl->seq_len, value_len, value_len, value);

            dl->u2_data_start = value - vb->txt_data.data;
            dl->u2_data_len   = value_len;
            vb->txt_section_bytes[SEC_QUAL_DATA] += dl->seq_len + 1; // +1 for \t
        }
        // All other subfields - have their own dictionary
        else        
            seg_one_subfield ((VBlock *)vb, value, value_len, vb_line_i, dict_id, SEC_SAM_OPTNL_SF_B250,
                             (value_len) + 1); // +1 for \t
    }

    if (separator == '\n') 
        vb->txt_section_bytes[SEC_SAM_OPTNL_SF_B250]--; // \n is already accounted for in SEC_SAM_OPTIONAL_B250
}

// calculate the expected length of SEQ and QUAL from the CIGAR string
// A CIGAR looks something like: "109S19M23S", See: https://samtools.github.io/hts-specs/SAMv1.pdf 
uint32_t seg_sam_seq_len_from_cigar (const char *cigar, unsigned cigar_len)
{
    // ZIP case: if the CIGAR is "*", we later get the length from SEQ and store it as eg "151*". 
    // In PIZ it will be eg "151*" or "1*" if both SEQ and QUAL are "*"
    if (cigar_len == 1 && cigar[0] == '*') return 0;

    unsigned seq_len=0, n=0;

    for (unsigned i=0; i < cigar_len; i++) {
        char c = cigar[i];
        if (IS_DIGIT (c)) 
            n = n*10 + (c - '0');
        
        else if (c=='M' || c=='I' || c=='S' || c=='=' || c=='X' || c=='*') { // these "consume" sequences
            ASSERT (n, "Error: Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            seq_len += n;
            n = 0;
        }

        else if (c=='D' || c=='N' || c=='H' || c=='P') { // these don't consume sequence
            ASSERT (n, "Error: Invalid CIGAR in %s: operation %c not preceded by a number. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
            n = 0;
        }

        else ABORT ("Error: Invalid CIGAR in %s: invalid operation %c. CIGAR=%.*s", 
                    txt_name, c, cigar_len, cigar);
    }                          

    ASSERT (!n, "Error: Invalid CIGAR in %s: expecting it to end with an operation character. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    ASSERT (seq_len, "Error: Invalid CIGAR in %s: CIGAR implies 0-length SEQ. CIGAR=%.*s", 
            txt_name, cigar_len, cigar);

    return seq_len;
}

static void seg_sam_seq_qual_fields (VBlockSAM *vb, ZipDataLineSAM *dl, unsigned vb_line_i)
{
    const char *seq  = &vb->txt_data.data[dl->seq_data_start];
    const char *qual = &vb->txt_data.data[dl->qual_data_start];

    bool seq_is_available  = dl->seq_data_len  != 1 || *seq  != '*';
    bool qual_is_available = dl->qual_data_len != 1 || *qual != '*';

    // handle the case where CIGAR is "*" and we get the dl->seq_len directly from SEQ or QUAL
    if (!dl->seq_len) { // CIGAR is not available
        ASSERT (!seq_is_available || !qual_is_available || dl->seq_data_len==dl->qual_data_len, 
                "Bad line in %s: SEQ length is %u, QUAL length is %u, unexpectedly differ. SEQ=%.*s QUAL=%.*s", 
                txt_name, dl->seq_data_len, dl->qual_data_len, dl->seq_data_len, seq, dl->qual_data_len, qual);    

        dl->seq_len = MAX (dl->seq_data_len, dl->qual_data_len); // one or both might be not available and hence =1
    
        char new_cigar[20];
        sprintf (new_cigar, "%u*", dl->seq_len); // eg "151*". if both SEQ and QUAL are unavailable it will be "1*"
        unsigned new_cigar_len = strlen (new_cigar);

        seg_one_field (vb, new_cigar, new_cigar_len, vb_line_i, SAM_CIGAR); 
        vb->txt_section_bytes[SEC_SAM_CIGAR_B250] -= (new_cigar_len-1); // adjust - actual SAM length was only one (seg_one_field added new_cigar_len)
    } 
    else { // CIGAR is available - just check the seq and qual lengths
        ASSERT (!seq_is_available || dl->seq_data_len == dl->seq_len, 
                "Bad line in %s: according to CIGAR, expecting SEQ length to be %u but it is %u. SEQ=%.*s", 
                txt_name, dl->seq_len, dl->seq_data_len, dl->seq_data_len, seq);

        ASSERT (!qual_is_available || dl->qual_data_len == dl->seq_len, 
                "Bad line in %s: according to CIGAR, expecting QUAL length to be %u but it is %u. QUAL=%.*s", 
                txt_name, dl->seq_len, dl->qual_data_len, dl->qual_data_len, qual);    
    }

    // byte counts for --show-sections 
    vb->txt_section_bytes[SEC_SEQ_DATA]  += dl->seq_data_len  + 1; // +1 for terminating \t
    vb->txt_section_bytes[SEC_QUAL_DATA] += dl->qual_data_len + 1; 
}

const char *seg_sam_data_line (VBlock *vb_,   
                               const char *field_start_line,     // index in vb->txt_data where this line starts
                               uint32_t vb_line_i) // line within this vb (starting from 0)
{
    VBlockSAM *vb = (VBlockSAM *)vb_;
    ZipDataLineSAM *dl = DATA_LINE (vb, vb_line_i);

    const char *next_field, *field_start=field_start_line;
    unsigned field_len=0;
    char separator;
    bool has_13 = false; // does this line end in Windows-style \r\n rather than Unix-style \n

    int32_t len = AFTERENT (char, vb->txt_data) - field_start_line;

    // QNAME - We break down the QNAME into subfields separated by / and/or : - these are vendor-defined strings. Examples:
    // Illumina: <instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> for example "A00488:61:HMLGNDSXX:4:1101:15374:1031" see here: https://help.basespace.illumina.com/articles/descriptive/fastq-files/
    // PacBio BAM: {movieName}/{holeNumber}/{qStart}_{qEnd} see here: https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
    field_start = field_start_line;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "QNAME");
    
    seg_compound_field ((VBlockP)vb, &vb->mtf_ctx[SAM_QNAME], field_start, field_len, &vb->qname_mapper,
                        dict_id_sam_qname_sf (dict_id_make ("Q0NAME", 6)), 0, 
                        SEC_SAM_QNAME_B250, SEC_SAM_QNAME_SF_B250, vb_line_i);

    // FLAG
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "FLAG");
    seg_one_field (vb, field_start, field_len, vb_line_i, SAM_FLAG);

    // RNAME
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "RNAME");
    uint32_t rname_node_index = seg_one_field (vb, field_start, field_len, vb_line_i, SAM_RNAME);

    random_access_update_chrom (vb_, vb_line_i, rname_node_index);

    // POS - two options:
    // 1. if RNAME is the same as the previous line - store the delta in SAM_POS dictionary
    // 2. If its a different RNAME - add to vb->random_pos_data (not a dictionary)
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "POS");
    
    if (rname_node_index != vb->last_rname_node_index)  // different rname than prev line - store full pos in random
        vb->last_pos = seg_sam_add_to_random_pos_data (vb, SEC_SAM_RAND_POS_DATA, field_start, field_len, field_len+1, "POS");

    else // same rname - do a delta
        vb->last_pos = seg_pos_field (vb_, vb->last_pos, SAM_POS, SEC_SAM_POS_B250, field_start, field_len, vb_line_i, "POS");

    random_access_update_pos (vb_, vb->last_pos);

    vb->last_rname_node_index = rname_node_index;

    // MAPQ
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "MAPQ");
    seg_one_field (vb, field_start, field_len, vb_line_i, SAM_MAPQ);

    // CIGAR - if CIGAR is "*" and we wait to get the length from SEQ or QUAL
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "CIGAR");
    dl->seq_len = seg_sam_seq_len_from_cigar (field_start, field_len);
    if (dl->seq_len) seg_one_field (vb, field_start, field_len, vb_line_i, SAM_CIGAR); // not "*" - all good!

    // RNEXT - add to RNAME dictionary
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "RNEXT");
    uint32_t rnext_node_index = seg_one_field (vb, field_start, field_len, vb_line_i, SAM_RNAME);
    
    bool rnext_same_as_rname = (field_len == 1 && (field_start[0] == '=' || field_start[0] == '*')) || 
                                (rnext_node_index == rname_node_index);
    // PNEXT - 2.5 options:
    // 1. if RNEXT is the same as RNAME or unknown (i.e. either it is "=" or "*" or the string is the same) - store the delta
    //    vs. POS in SAM_PNEXT dictionary
    //    1A. If 1, but PNEXT is "0" and hence "unavailable" per SAM specification - store "*" in the SAM_PNEXT dictionary
    // 2. If RNEXT and RNAME are different - add to vb->random_pos_data (not a dictionary)
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "PNEXT");
        
    if (rnext_same_as_rname) { // RNAME and RNEXT are the same - store delta in dictionary
        if (field_len == 1 && *field_start == '0') // PNEXT is "unavailable"
            seg_one_field (vb, "*", 1, vb_line_i, SAM_PNEXT); 
        else
            seg_pos_field (vb_, vb->last_pos, SAM_PNEXT, SEC_SAM_PNEXT_B250, field_start, field_len, vb_line_i, "PNEXT");
    }

    else  // RNAME and RNEXT differ - store in random_pos
        seg_sam_add_to_random_pos_data (vb, SEC_SAM_RAND_POS_DATA, field_start, field_len, field_len+1, "PNEXT");

    // TLEN - replace a value that is the negative of the previous line with "*"
    field_start = next_field;
    next_field = seg_get_next_item (field_start, &len, false, true, false, vb_line_i, &field_len, &separator, &has_13, "TLEN");
    seg_sam_tlen_field (vb, field_start, field_len, vb_line_i);

    // SEQ & QUAL
    dl->seq_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, false, true, false, vb_line_i, &dl->seq_data_len, &separator, &has_13, "SEQ");

    dl->qual_data_start = next_field - vb->txt_data.data;
    next_field = seg_get_next_item (next_field, &len, true, true, false, vb_line_i, &dl->qual_data_len, &separator, &has_13, "QUAL");

    seg_sam_seq_qual_fields (vb, dl, vb_line_i); // also updates cigar if it is "*"

    // OPTIONAL fields - up to MAX_SUBFIELDS of them
    char oname[MAX_SUBFIELDS * 5 + 1]; // each name is 5 characters per SAM specification, eg "MC:Z:" ; +1 for # in case of \r
    unsigned oname_len;

    for (oname_len=0; oname_len < MAX_SUBFIELDS*5 && separator != '\n'; oname_len += 5) {
        field_start = next_field;
        next_field = seg_get_next_item (field_start, &len, true, true, false, vb_line_i, &field_len, &separator, &has_13, "OPTIONAL-subfield");
        seg_sam_optional_field (vb, dl, field_start, field_len, separator, vb_line_i);
        memcpy (&oname[oname_len], field_start, 5);
    }
    ASSERT (separator=='\n', "Error in %s: too many optional fields, limit is %u", txt_name, MAX_SUBFIELDS);

    // treat a Windows style \r (as part of a \r\n line termination) as part of the oname - add a # to represent it
    if (has_13) oname[oname_len++] = '#';

    // if oname is empty, we make it "*" (and adjust size accounting)
    if (!oname_len) oname[oname_len++] = '*';

    seg_one_field (vb, oname, oname_len, vb_line_i, SAM_OPTIONAL);

    // if we have no real OPTIONAL, seg_one_field adds +1 for the character (correct for '#' as it stands for \r, but incorrect for '*') and +1 for the wrongly presumsed \t
    if (oname_len == 1) 
        vb->txt_section_bytes[SEC_SAM_OPTIONAL_B250] -= (1 + (oname[0]=='*'));  

    return next_field;
}

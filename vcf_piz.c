// ------------------------------------------------------------------
//   vcf_piz.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <math.h>
#include "vcf_private.h"
#include "zfile.h"
#include "txtfile.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "endianness.h"
#include "sections.h"
#include "random_access.h"
#include "dict_id.h"
#include "strings.h"
#include "codec.h"
#include "reference.h"
#include "reconstruct.h"

// returns true if section is to be skipped reading / uncompressing
bool vcf_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (flag.drop_genotypes && // note: if all samples are filtered out with --samples then flag.drop_genotypes=true (set in samples_digest_vcf_header)
        (dict_id.num == dict_id_fields[VCF_FORMAT] || dict_id.num == dict_id_fields[VCF_SAMPLES] || dict_id_is_vcf_format_sf (dict_id)))
        return true;

    if (flag.gt_only && (dict_id_is_vcf_format_sf (dict_id) && dict_id.num != dict_id_FORMAT_GT))
        return true;

    return false;
}

CONTAINER_FILTER_FUNC (vcf_piz_filter)
{
    if (dict_id.num == dict_id_fields[VCF_SAMPLES]) {
        if (item < 0)  // filter for repeat
            return samples_am_i_included (rep); 

        else // filter for item
            if (flag.gt_only) return con->items[item].dict_id.num == dict_id_FORMAT_GT;
    }

    return true;    
}

// FORMAT - obey GT-only and drop-genotypes ; load haplotype line
SPECIAL_RECONSTRUCTOR (vcf_piz_special_FORMAT)
{
    VBlockVCF *vb_vcf = (VBlockVCF *)vb;

    if (flag.drop_genotypes) goto done;

    bool has_GT = (snip_len>=2 && snip[0]=='G' && snip[1] == 'T' && (snip_len==2 || snip[2] == ':'));
    
    if (reconstruct) {
        if (flag.gt_only) {
            if (has_GT)
                RECONSTRUCT ("GT\t", 3)
        }
        else 
            RECONSTRUCT (snip, snip_len);
    }

    // initialize haplotype stuff
    if (has_GT && !vb_vcf->ht_matrix_ctx) {

        ASSERT ((vb_vcf->ht_matrix_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_HT)), 
                "Error in vcf_piz_special_FORMAT vb_i=%u: cannot find GT_HT data", vb->vblock_i);

        // will exist in case of use of CODEC_HAPMAT but not CODEC_GTSHARK        
        vb_vcf->hapmat_index_ctx = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT_HT_INDEX);
        
        if (vb_vcf->hapmat_index_ctx)
            codec_hapmat_piz_calculate_columns (vb);
    }
done:
    return false; // no new value
}

// REFALT - reconstruct from reference and/or common SNPs
SPECIAL_RECONSTRUCTOR (vcf_piz_special_REFALT)
{
    if (!reconstruct) goto done;

    ASSERT (snip_len==2, "Error in vcf_piz_special_REFALT: expecting snip_len=2 but seeing %u", snip_len);

    // snip is 3 characters - REF, \t, ALT
    char ref_alt[3] = { 0, '\t', 0 };
    char ref_value = 0;
    
    PosType pos = vb->contexts[VCF_POS].last_value.i;

    if (snip[0] == '-' || snip[1] == '-') { 
        const Range *range = ref_piz_get_range (vb, pos, 1);
        ASSERT (range, "Error in vcf_piz_special_REFALT: failed to find range for chrom='%s' pos=%"PRId64, vb->chrom_name, pos);
        
        uint32_t idx = pos - range->first_pos;
        ASSERT (ref_is_nucleotide_set (range, idx), "Error in vcf_piz_special_REFALT: reference is not set: chrom=%.*s pos=%"PRId64, range->chrom_name_len, range->chrom_name, pos);
        ref_value = ref_get_nucleotide (range, idx);
    }

    // recover ref
    if (snip[0] == '-') 
        ref_alt[0] = ref_value;
    else 
        ref_alt[0] = snip[0];

    // recover alt
    if (snip[1] == '+') { // the alt has the most common value for a SNP
        if      (ref_alt[0] == 'A') ref_alt[2] = 'G';
        else if (ref_alt[0] == 'C') ref_alt[2] = 'T';
        else if (ref_alt[0] == 'G') ref_alt[2] = 'A';
        else if (ref_alt[0] == 'T') ref_alt[2] = 'C';
    }
    else if (snip[1] == '-')  // the alt has the reference value
        ref_alt[2] = ref_value;

    else // the alt is specified verbatim
        ref_alt[2] = snip[1];

    RECONSTRUCT (ref_alt, sizeof (ref_alt));

done:
    return false; // no new value
}   

SPECIAL_RECONSTRUCTOR (vcf_piz_special_AC)
{
    if (!reconstruct) goto done;
    
    bool is_an_before_ac = (bool)(snip[0] - '0');
    bool is_af_before_ac = (bool)(snip[1] - '0');

    Context *ctx_an = ctx_get_existing_ctx (vb, dict_id_INFO_AN);
    Context *ctx_af = ctx_get_existing_ctx (vb, dict_id_INFO_AF);

    uint32_t an = is_an_before_ac ? ctx_an->last_value.i : ctx_peek_next_int (vb, ctx_an);
    double   af = is_af_before_ac ? ctx_af->last_value.f : ctx_peek_next_float (vb, ctx_af);

    char ac_str[30];
    unsigned ac_str_len = str_int ((int64_t)round(an * af), ac_str);    

    RECONSTRUCT (ac_str, ac_str_len);

done:
    return false; // no new value
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_DS)
{
    if (!reconstruct) goto done;

    Context *ctx_gt = ctx_get_existing_ctx (vb, dict_id_FORMAT_GT);

    // we are guaranteed that if we have a special snip, then all values are either '0' or '1';
    char *gt = ENT (char, vb->txt_data, ctx_gt->last_txt);
    unsigned dosage=0;
    for (unsigned i=0; i < ctx_gt->last_txt_len; i+=2) 
        dosage += gt[i]-'0';

    char float_format[10];
    int32_t val;
    sscanf (snip, "%s %d", float_format, &val); // snip looks like eg: "%5.3f 50000"

    bufprintf (vb, &vb->txt_data, float_format, (double)val / 1000000 + dosage);

done:
    return false; // no new value
}

SPECIAL_RECONSTRUCTOR (vcf_piz_special_BaseCounts)
{
    Context *ctx_refalt = &vb->contexts[VCF_REFALT];
    ASSERTE (ctx_refalt->last_txt_len == 3, "Expecting ctx_refalt->last_txt_len=%u to be 3", ctx_refalt->last_txt_len);
    const char *refalt = ENT (const char, vb->txt_data, ctx_refalt->last_txt);

    uint32_t counts[4], sorted_counts[4] = {}; // counts of A, C, G, T

    new_value->i = 0;
    char *str = (char *)snip;

    for (unsigned i=0; i < 4; i++) {
        sorted_counts[i] = strtoul (str, &str, 10);
        str++; // skip comma seperator
        new_value->i += sorted_counts[i];
    }

    if (!reconstruct) goto done; // just return the new value

    ASSERTE (str - snip == snip_len + 1, "expecting (str-snip)=%d == (snip_len+1)=%u", (int)(str - snip), snip_len+1);

    unsigned ref_i = acgt_encode[(int)refalt[0]];
    unsigned alt_i = acgt_encode[(int)refalt[2]];
    
    counts[ref_i] = sorted_counts[0];
    counts[alt_i] = sorted_counts[1];
    
    unsigned sc_i=2;
    for (unsigned i=0; i <= 3; i++)
        if (ref_i != i && alt_i != i) counts[i] = sorted_counts[sc_i++];

    bufprintf (vb, &vb->txt_data, "%u,%u,%u,%u", counts[0], counts[1], counts[2], counts[3]);

done:
    return true; // has new value
}

// the case where SVLEN is minus the delta between END and POS
SPECIAL_RECONSTRUCTOR (vcf_piz_special_SVLEN)
{
    if (!reconstruct) goto done;

    int64_t value = -vb->contexts[VCF_POS].last_delta; // END is a alias of POS - they share the same data stream - so last_delta would be the delta between END and POS
    char str[30];
    unsigned str_len = str_int (value, str);
    RECONSTRUCT (str, str_len);

done:
    return false; // no new value
}


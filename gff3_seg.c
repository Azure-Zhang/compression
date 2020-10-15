// ------------------------------------------------------------------
//   seg_gff3.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "genozip.h"
#include "seg.h"
#include "vblock.h"
#include "context.h"
#include "random_access.h"
#include "file.h"
#include "strings.h"
#include "optimize.h"
#include "piz.h"
#include "dict_id.h"

#define MAX_ENST_ITEMS 10 // maximum number of items in an enst structure. this can be changed without impacting backward compatability.

// called from seg_all_data_lines
void gff3_seg_initialize (VBlock *vb)
{
    vb->contexts[GFF3_SEQID].inst = CTX_INST_NO_STONS; // needs b250 node_index for random access
    vb->contexts[GFF3_ATTRS].inst = CTX_INST_NO_STONS;
}

void gff3_seg_finalize (VBlockP vb)
{
    // top level snip
    Structured top_level = { 
        .repeats   = vb->lines.len,
        .flags     = STRUCTURED_TOPLEVEL,
        .num_items = 10,
        .items     = { { (DictId)dict_id_fields[GFF3_SEQID],  DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_SOURCE], DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_TYPE],   DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_START],  DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_END],    DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_SCORE],  DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_STRAND], DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_PHASE],  DID_I_NONE, "\t" },
                       { (DictId)dict_id_fields[GFF3_ATTRS],  DID_I_NONE, ""   },
                       { (DictId)dict_id_fields[GFF3_EOL],    DID_I_NONE, ""   } }
    };

    seg_structured_by_ctx (vb, &vb->contexts[GFF3_TOPLEVEL], &top_level, 0, 0, 0);
}


// returns length of next expected item, and 0 if unsuccessful
static unsigned gff3_seg_get_aofs_item_len (const char *str, unsigned len, bool is_last_item)
{
    unsigned i=0; for (; i < len; i++) {
        bool is_comma = str[i] == ',';
        bool is_space = str[i] == ' ';

        if ((is_last_item && is_comma) || (!is_last_item && is_space)) return i; // we reached the end of the item - return its length

        if ((is_last_item && is_space) || (!is_last_item && is_comma)) return 0; // invalid string - unexpected seperator
    }

    if (is_last_item) return i; // last item of last entry
    else return 0; // the string ended prematurely - this is not yet the last item
}

// a field that looks like: "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
// we have an array (2 entires in this example) of items (4 in this examples) - the entries are separated by comma and the items by space
// observed in Ensembel generated GVF: Variant_effect, sift_prediction, polyphen_prediction, variant_peptide
// The last item is treated as an ENST_ID (format: ENST00000399012) while the other items are regular dictionaries
// the names of the dictionaries are the same as the ctx, with the 2nd character replaced by 1,2,3...
// the field itself will contain the number of entries
static void gff3_seg_array_of_struct (VBlock *vb, Context *subfield_ctx, 
                                      Structured st, 
                                      const char *snip, unsigned snip_len)
{
    bool is_last_entry = false;

    // get ctx's
    Context *ctxs[MAX_ENST_ITEMS] = {}; // an array of length num_items_in_struct (pointer to start of sub-array in vb->contexts)
    for (unsigned i=0; i < st.num_items; i++) 
        ctxs[i] = mtf_get_ctx (vb, st.items[i].dict_id); 

    // set roll back point
    uint64_t saved_mtf_i_len[MAX_ENST_ITEMS], saved_local_len[MAX_ENST_ITEMS], saved_txt_len[MAX_ENST_ITEMS];
    for (unsigned item_i=0; item_i < st.num_items; item_i++) {
        saved_mtf_i_len[item_i] = ctxs[item_i]->mtf_i.len;
        saved_local_len[item_i] = ctxs[item_i]->local.len;
        saved_txt_len  [item_i] = ctxs[item_i]->txt_len;
    }
    const char *saved_snip = snip;
    unsigned saved_snip_len = snip_len;

    st.repeats = 0;

    while (snip_len) {
        
        for (unsigned item_i=0; item_i < st.num_items; item_i++) {
            bool is_last_item = (item_i == st.num_items-1);
            unsigned item_len = gff3_seg_get_aofs_item_len (snip, snip_len, is_last_item);
            if (!item_len) goto badly_formatted;

            if (!is_last_item)
                seg_by_dict_id (vb, snip, item_len, ctxs[item_i]->dict_id, item_len + (st.items[item_i].seperator[0] != 0) + (st.items[item_i].seperator[1] != 0));
            else {
                is_last_entry = (snip_len - item_len == 0);
                seg_id_field ((VBlockP)vb, (DictId)dict_id_ENSTid, snip, item_len, !is_last_entry);
            }
    
            snip     += item_len + 1 - is_last_entry; // 1 for either the , or the ' ' (except in the last item of the last entry)
            snip_len -= item_len + 1 - is_last_entry;
        }

        if (!is_last_entry && snip[-1]!=',') goto badly_formatted; // expecting a , after the end of all items in this entry
        
        st.repeats++;

        ASSSEG (st.repeats <= STRUCTURED_MAX_REPEATS, snip, "Error in gff3_seg_array_of_struct - exceeded maximum repeats allowed (%lu) while parsing %s",
                STRUCTURED_MAX_REPEATS, subfield_ctx->name);
    }

    // finally, the "structured" snip itself
    seg_structured_by_ctx ((VBlockP)vb, subfield_ctx, &st, NULL, 0, 0);

    return;

badly_formatted:
    // roll back all the changed data
    for (unsigned item_i=0; item_i < st.num_items ; item_i++) {
        ctxs[item_i]->mtf_i.len = saved_mtf_i_len[item_i];
        ctxs[item_i]->local.len = saved_local_len[item_i];
        ctxs[item_i]->txt_len   = saved_txt_len[item_i];
    }

    // now save the entire snip in the dictionary
    seg_by_dict_id (vb, saved_snip, saved_snip_len, subfield_ctx->dict_id, saved_snip_len); 
}                           

static bool gff3_seg_special_info_subfields (VBlockP vb, DictId dict_id, const char **this_value, unsigned *this_value_len, char *optimized_snip)
{
    // ID - this is a sequential number (at least in GRCh37/38)
    if (dict_id.num == dict_id_ATTR_ID) {
        Context *ctx = mtf_get_ctx (vb, dict_id);
        seg_pos_field ((VBlockP)vb, ctx->did_i, ctx->did_i, true, *this_value, *this_value_len, false);
        return false; // do not add to dictionary/b250 - we already did it
    }

    // Dbxref (example: "dbSNP_151:rs1307114892") - we divide to the non-numeric part which we store
    // in a dictionary and the numeric part which store in a NUMERICAL_ID_DATA section
    if (dict_id.num == dict_id_ATTR_Dbxref) {
        seg_id_field (vb, dict_id, *this_value, *this_value_len, false); // discard the const as seg_id_field modifies
        return false; // do not add to dictionary/b250 - we already did it
    }

    // subfields that are arrays of structs, for example:
    // "non_coding_transcript_variant 0 ncRNA ENST00000431238,intron_variant 0 primary_transcript ENST00000431238"
    if (dict_id.num == dict_id_ATTR_Variant_effect) {
        static const Structured Variant_effect = {
            .num_items   = 4, 
            .flags       = STRUCTURED_DROP_FINAL_ITEM_SEP,
            .repsep      = {0,0},
            .items       = { { .dict_id={.id="V0arEff" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="V1arEff" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="V2arEff" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="ENSTid"  }, .seperator = {','}, .did_i = DID_I_NONE } }
        };
        gff3_seg_array_of_struct (vb, mtf_get_ctx (vb, dict_id), Variant_effect, *this_value, *this_value_len);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_sift_prediction) {
        static const Structured sift_prediction = {
            .num_items   = 4, 
            .flags       = STRUCTURED_DROP_FINAL_ITEM_SEP,
            .repsep      = {0,0},
            .items       = { { .dict_id={.id="S0iftPr" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="S1iftPr" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="S2iftPr" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="ENSTid"  }, .seperator = {','}, .did_i = DID_I_NONE } }
        };
        gff3_seg_array_of_struct (vb, mtf_get_ctx (vb, dict_id), sift_prediction, *this_value, *this_value_len);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_polyphen_prediction) {
        static const Structured polyphen_prediction = {
            .num_items   = 4, 
            .flags       = STRUCTURED_DROP_FINAL_ITEM_SEP,
            .repsep      = {0,0},
            .items       = { { .dict_id={.id="P0olyPhP" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="P1olyPhP" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="P2olyPhP" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="ENSTid"   }, .seperator = {','}, .did_i = DID_I_NONE } }
        };
        gff3_seg_array_of_struct (vb, mtf_get_ctx (vb, dict_id), polyphen_prediction, *this_value, *this_value_len);
        return false;
    }

    if (dict_id.num == dict_id_ATTR_variant_peptide) {
        static const Structured variant_peptide = {
            .num_items   = 3, 
            .flags       = STRUCTURED_DROP_FINAL_ITEM_SEP,
            .repsep      = {0,0},
            .items       = { { .dict_id={.id="v0arPep" }, .seperator = {' '}, .did_i = DID_I_NONE }, // small v to differentiate from Variant_effect, so that dict_id to did_i mapper can map both
                             { .dict_id={.id="v1arPep" }, .seperator = {' '}, .did_i = DID_I_NONE },
                             { .dict_id={.id="ENSTid"  }, .seperator = {','}, .did_i = DID_I_NONE } }
        };
        gff3_seg_array_of_struct (vb, mtf_get_ctx (vb, dict_id), variant_peptide, *this_value, *this_value_len);
        return false;
    }

    // we store these 3 in one dictionary, as they are correlated and will compress better together
    if (dict_id.num == dict_id_ATTR_Variant_seq   ||
        dict_id.num == dict_id_ATTR_Reference_seq ||
        dict_id.num == dict_id_ATTR_ancestral_allele) {

        // note: all three are stored together in dict_id_ATTR_Reference_seq as they are correlated
        Context *ctx = mtf_get_ctx (vb, dict_id_ATTR_Reference_seq); 
        seg_add_to_local_text (vb, ctx, *this_value, *this_value_len, *this_value_len);
        
        return false; // do not add to dictionary/b250 - we already did it
    }

    // Optimize Variant_freq
    unsigned optimized_snip_len;
    if (flag_optimize_Vf && (dict_id.num == dict_id_ATTR_Variant_freq) &&
        optimize_float_2_sig_dig (*this_value, *this_value_len, 0, optimized_snip, &optimized_snip_len)) {
        
        vb->vb_data_size -= (int)(*this_value_len) - (int)optimized_snip_len;
        *this_value = optimized_snip;
        *this_value_len = optimized_snip_len;
        return true; // proceed with adding to dictionary/b250
    }

    return true; // all other cases -  procedue with adding to dictionary/b250
}

const char *gff3_seg_txt_line (VBlock *vb, const char *field_start_line, bool *has_13)     // index in vb->txt_data where this line starts
{
    const char *next_field=field_start_line, *field_start;
    unsigned field_len=0;
    char separator;

    int32_t len = &vb->txt_data.data[vb->txt_data.len] - field_start_line;

    GET_NEXT_ITEM ("SEQID");
    seg_chrom_field (vb, field_start, field_len);

    SEG_NEXT_ITEM (GFF3_SOURCE);
    SEG_NEXT_ITEM (GFF3_TYPE);

    GET_NEXT_ITEM ("START");
    seg_pos_field (vb, GFF3_START, GFF3_START, false, field_start, field_len, true);
    random_access_update_pos (vb, GFF3_START);

    GET_NEXT_ITEM ("END");
    seg_pos_field (vb, GFF3_END, GFF3_START, false, field_start, field_len, true);

    SEG_NEXT_ITEM (GFF3_SCORE);
    SEG_NEXT_ITEM (GFF3_STRAND);
    SEG_NEXT_ITEM (GFF3_PHASE);

    if (separator != '\n') {
        GET_LAST_ITEM (DTF(names)[GFF3_ATTRS] /* pointer to string to allow pointer comparison */); 
        seg_info_field (vb, gff3_seg_special_info_subfields, field_start, field_len);
    }
    else
        seg_by_did_i (vb, NULL, 0, GFF3_ATTRS, 0); // NULL=MISSING so previous \t is removed
             
    SEG_EOL (GFF3_EOL, false);

    return next_field;
}


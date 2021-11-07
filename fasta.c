// ------------------------------------------------------------------
//   fasta.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include "fasta_private.h"
#include "seg.h"
#include "context.h"
#include "file.h"
#include "piz.h"
#include "dict_id.h"
#include "random_access.h"
#include "strings.h"
#include "regions.h"
#include "codec.h"
#include "stats.h"
#include "reconstruct.h"
#include "kraken.h"
#include "reference.h"
#include "segconf.h"
#include "chrom.h"
#include "tokenizer.h"

#define dict_id_is_fasta_desc_sf dict_id_is_type_1
#define dict_id_fasta_desc_sf dict_id_type_1

typedef struct {
    uint32_t seq_data_start, seq_len; // regular fasta and make-reference: start & length within vb->txt_data
} ZipDataLineFASTA;
#define DATA_LINE(i) ENT (ZipDataLineFASTA, vb->lines, (i))

unsigned fasta_vb_size (DataType dt) 
{ 
    return dt == DT_REF && command == PIZ ? sizeof (VBlock) : sizeof (VBlockFASTA); 
}

unsigned fasta_vb_zip_dl_size (void) { return sizeof (ZipDataLineFASTA); }

void fasta_vb_release_vb (VBlockFASTA *vb)
{
    if (VB_DT(DT_REF) && command == PIZ) return; // this is actually a VBlock, not VBlockFASTA
    
    memset ((char *)vb + sizeof (VBlock), 0, sizeof (VBlockFASTA) - sizeof (VBlock)); // zero all data unique to VBlockFASTA
    CTX(FASTA_NONREF)->local.len = 0; // len might be is used even though buffer is not allocated (in make-ref)
}

void fasta_vb_destroy_vb (VBlockFASTA *vb) {}

// used by ref_make_create_range
void fasta_get_data_line (VBlockP vb_, uint32_t line_i, uint32_t *seq_data_start, uint32_t *seq_len)
{
    VBlockFASTA *vb = (VBlockFASTA *)vb_;

    *seq_data_start = DATA_LINE(line_i)->seq_data_start;
    *seq_len        = DATA_LINE(line_i)->seq_len;
}

//-------------------------
// TXTFILE stuff
//-------------------------

// returns true if txt_data[txt_i] is the end of a FASTA contig (= next char is '>' or end-of-file), false if not, 
// and -1 if more data (lower first_i) is needed 
static inline int fasta_is_end_of_contig (VBlock *vb, uint32_t first_i,
                                          int32_t txt_i) // index of a \n in txt_data
{
    ARRAY (char, txt, vb->txt_data);

    // if we're not at the end of the data - we can just look at the next character
    if (txt_i < vb->txt_data.len-1)
        return txt[txt_i+1] == '>';

    // if we're at the end of the line, we scan back to the previous \n and check if it NOT the >
    bool newline_run = true;
    for (int32_t i=txt_i-1; i >= first_i; i--) {
        if (txt[i] == '\n' && !newline_run) // ASSUMES NO COMMENT LINES (starting with ;), To do: fix this
            return txt[i+1] != '>'; // this row is a sequence row, not a description row
        
        else if (newline_run && txt[i] != '\n' && txt[i] != '\r')
            newline_run = false; // newline run ends when we encounter the first non newline
    }

    return -1; // we need more data (lower first_i)
}

// returns the length of the data at the end of vb->txt_data that will not be consumed by this VB is to be passed to the next VB
int32_t fasta_unconsumed (VBlockP vb, uint32_t first_i, int32_t *last_i)
{
    bool is_entire_vb = (first_i == 0 && *last_i == vb->txt_data.len-1);

    ASSERT (*last_i >= 0 && *last_i < vb->txt_data.len, "*last_i=%d is out of range [0,%"PRIu64"]", *last_i, vb->txt_data.len);

    ARRAY (char, txt, vb->txt_data);

    // case: reference file - we allow only one contig (or part of it) per VB - move second contig onwards to next vb
    // (note: first_i=0 when flag.make_reference)
    if (flag.make_reference) {
        bool data_found = false;
        for (uint32_t i=0; i < vb->txt_data.len; i++) {
            // just don't allow now-obsolete ';' rather than trying to disentangle comments from descriptions
            ASSINP (txt[i] != ';', "Error: %s contains a ';' character - this is not supported for reference files. Contig descriptions must begin with a >", txt_name);
        
            // if we've encountered a new DESC line after already seeing sequence data, move this DESC line and
            // everything following to the next VB
            if (data_found && txt[i]=='>' && txt[i-1]=='\n') 
                return vb->txt_data.len - i;

            if (!data_found && (txt[i] != '\n' && txt[i] != '\r')) data_found = true; // anything, except for empty lines, is considered data
        }
    }

    // we move the final partial line to the next vb (unless we are already moving more, due to a reference file)
    for (int32_t i=*last_i; i >= (int32_t)first_i; i--) {

        if (txt[i] == '\n') {

            // when compressing FASTA with a reference or --multifasta - an "end of line" means "end of contig" - 
            // i.e. one that the next character is >, or it is the end of the file
            // note: when compressing FASTA with a reference (eg long reads stored in a FASTA instead of a FASTQ), 
            // line cannot be too long - they must fit in a VB
            if ((flag.reference & REF_ZIP_LOADED) || flag.multifasta) {
                int is_end_of_contig = fasta_is_end_of_contig (vb, first_i, i);

                switch (is_end_of_contig) {
                    case true  : break;            // end of contig
                    case false : continue;         // not end of contig
                    default    : goto out_of_data; // need more data (lower first_i)
                }
            }

            // otherwise - tolerate a VB that ends part way through a SEQ
            else if (is_entire_vb && i+1 < vb->txt_data.len &&
                     txt[i+1] != ';' && txt[i+1] != '>') { // partial line isn't a Description or a Comment, hence its a Sequence
                ((VBlockFASTA *)vb)->vb_has_no_newline = true;
                return 0;                
            }
                
            *last_i = i;
            return vb->txt_data.len-1 - i;
        }
    }

out_of_data:
    // case: an entire FASTA VB without newlines (i.e. a very long sequential SEQ) - we accept a VB without newlines and deal with it in Seg
    if (is_entire_vb && !segconf.running) {
        ((VBlockFASTA *)vb)->vb_has_no_newline = true;
        return 0;
    }

    // case: a single BGZF block without newlines - test next block
    else 
        return -1; // cannot find end-of-line in the data starting first_i
}

//-------------------------
// SEG & ZIP stuff
//-------------------------

// called by main thread at the beginning of zipping this file
void fasta_zip_initialize (void)
{
}

// callback function for compress to get data of one line (called by codec_lzma_data_in_callback)
COMPRESSOR_CALLBACK (fasta_zip_seq)
{
    ZipDataLineFASTA *dl = DATA_LINE (vb_line_i);

    // note: maximum_len might be shorter than the data available if we're just sampling data in zip_assign_best_codec
    *line_data_len = MIN_(dl->seq_len, maximum_size);
    
    if (line_data) // if NULL, only length was requested
        *line_data = dl->seq_len ? ENT (char, vb->txt_data, dl->seq_data_start) : NULL;
}   

void fasta_seg_initialize (VBlock *vb)
{
    START_TIMER;

    ASSINP (vb->vblock_i > 1 || *FIRSTENT (char, vb->txt_data) == '>' || *FIRSTENT (char, vb->txt_data) == ';',
            "Error: expecting FASTA file %s to start with a '>' or a ';'", txt_name);

    CTX(FASTA_TOPLEVEL)->no_stons  = true; // keep in b250 so it can be eliminated as all_the_same
    CTX(FASTA_CONTIG)->flags.store = STORE_INDEX; // since v12
    CTX(FASTA_CONTIG)->no_stons    = true; // needs b250 node_index for reference
    CTX(FASTA_LINEMETA)->no_stons  = true; // avoid edge case where entire b250 is moved to local due to singletons, because fasta_reconstruct_vb iterates on ctx->b250

    if (kraken_is_loaded) {
        CTX(FASTA_TAXID)->flags.store    = STORE_INT;
        CTX(FASTA_TAXID)->no_stons       = true; // must be no_stons the SEC_COUNTS data needs to mirror the dictionary words
        CTX(FASTA_TAXID)->counts_section = true; 
    }

    // if this neocleotide FASTA of unrelated contigs, we're better off with ACGT        
    if (!flag.multifasta && segconf.seq_type == SQT_NUKE)
        codec_acgt_comp_init (VB);

    // if the contigs in this FASTA are related, let codec_assign_best_codec assign the bext codec 
    else 
        CTX(FASTA_NONREF)->ltype  = LT_SEQUENCE;

    if (flag.reference & REF_ZIP_LOADED) 
        CTX(FASTA_NONREF)->no_callback = true; // override callback if we are segmenting to a reference

    // in --stats, consolidate stats into FASTA_NONREF
    stats_set_consolidation (vb, FASTA_NONREF, 1, FASTA_NONREF_X);

    COPY_TIMER (seg_initialize);
}

void fasta_seg_finalize (VBlockP vb)
{
    // top level snip
    SmallContainer top_level = { 
        .repeats      = vb->lines.len,
        .is_toplevel  = true,
        .callback     = true,
        .nitems_lo    = 2,
        .items        = { { .dict_id = { _FASTA_LINEMETA }  },
                          { .dict_id = { _FASTA_EOL }, .translator = FASTA2PHYLIP_EOL } }
    };

    container_seg (vb, CTX(FASTA_TOPLEVEL), (ContainerP)&top_level, 0, 0, 0);

    // decide whether the sequences in this FASTA represent contigs (in which case we want a FASTA_CONTIG dictionary
    // and random access) or do they represent reads (in which case they are likely to numerous to be added to a dict)
    if (segconf.running) {
        uint64_t num_contigs_this_vb = CTX(FASTA_CONTIG)->nodes.len;
        ASSINP0 (num_contigs_this_vb, "Invalid FASTA file: no sequence description line");

        uint64_t avg_contig_size_this_vb = vb->txt_data.len / num_contigs_this_vb;
        uint64_t est_num_contigs_in_file = txtfile_get_seggable_size() / avg_contig_size_this_vb;

        // limit the number of contigs, to avoid the FASTA_CONTIG dictionary becoming too big. note this also
        // sets a limit for fasta-to-phylip translation
        #define MAX_CONTIGS_IN_FILE 1000000 
        segconf.fasta_has_contigs = num_contigs_this_vb == 1 || // the entire VB is a single contig
                                    est_num_contigs_in_file <  MAX_CONTIGS_IN_FILE; 

        // case: we've seen only characters that are both nucleotide and protein (as are A,C,G,T,N) - call it as nucleotide
        if (segconf.seq_type == SQT_NUKE_OR_AMINO) segconf.seq_type = SQT_NUKE;
    }
}

bool fasta_seg_is_small (ConstVBlockP vb, DictId dict_id)
{
    return dict_id.num == _FASTA_TOPLEVEL ||
           dict_id.num == _FASTA_DESC     ||
           dict_id.num == _FASTA_LINEMETA ||
           dict_id.num == _FASTA_TAXID    ||
           dict_id.num == _FASTA_EOL;
}

// description line - we segment it to its components
// note: we store the DESC container in its own ctx rather than just directly in LINEMETA, to make it easier to grep
static void fasta_seg_desc_line (VBlockFASTA *vb, const char *line_start, uint32_t line_len, bool *has_13)
{
    SAFE_NUL (&line_start[line_len]);
    
    // we store the contig name in a dictionary only (no b250), to be used if this fasta is used as a reference
    const char *chrom_name = line_start + 1;
    unsigned chrom_name_len = strcspn (line_start + 1, " \t\r\n");

    ASSSEG0 (chrom_name_len, line_start, "contig is missing a name");

    if (!flag.make_reference) {
        tokenizer_seg (VB, CTX(FASTA_DESC), line_start, line_len, sep_with_space, 0);
        
        char special_snip[100]; unsigned special_snip_len = sizeof (special_snip);
        seg_prepare_snip_other_do (SNIP_REDIRECTION, (DictId)_FASTA_DESC, false, 0, &special_snip[2], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_DESC;

        seg_by_did_i (VB, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    // case make_ref: add contig metadata (the rest of the line, except for the chrom_name)
    else {
        const char *md_start = chrom_name + chrom_name_len + strspn (&chrom_name[chrom_name_len], " \t");
        unsigned md_len = MIN_(strcspn (md_start, "\n\r"), REFCONTIG_MD_LEN-1);
        memcpy (vb->contig_metadata.str, md_start, md_len);
        vb->has_contig_metadata = true;        
    }

    // add contig to CONTIG dictionary (but not b250) and verify that its unique
    if (segconf.fasta_has_contigs || flag.make_reference || segconf.running) {
        bool is_new;
        chrom_seg_no_b250 (VB, STRa(chrom_name), &is_new);

        vb->ra_initialized = true;

        ASSINP (is_new || segconf.running, "Error: bad FASTA file - sequence \"%.*s\" appears more than once%s", chrom_name_len, chrom_name,
                flag.bind ? " (possibly in another FASTA being bound)" : 
                (flag.reference & REF_ZIP_LOADED) ? " (possibly the sequence size exceeds vblock size, try enlarging with --vblock)" : "");
    }

    vb->last_line = FASTA_LINE_DESC;    
    SAFE_RESTORE;
}

static void fast_seg_comment_line (VBlockFASTA *vb, const char *line_start, uint32_t line_len, bool *has_13)
{
    if (!flag.make_reference) {
        seg_add_to_local_text (VB, CTX(FASTA_COMMENT), line_start, line_len, line_len); 

        char special_snip[100]; unsigned special_snip_len = sizeof (special_snip);
        seg_prepare_snip_other_do (SNIP_OTHER_LOOKUP, (DictId)_FASTA_COMMENT, false, 0, &special_snip[2], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_COMMENT;

        seg_by_did_i (VB, special_snip, special_snip_len+2, FASTA_LINEMETA, 0);
        SEG_EOL (FASTA_EOL, true);
    }

    vb->last_line = FASTA_LINE_COMMENT;
}

// ZIP: main thread during segconf.running
static SeqType fasta_get_seq_type (STRp(seq))
{
    // we determine the type by characters that discriminate between protein and nucleotides, 
    // according to: https://www.bioinformatics.org/sms/iupac.html and https://en.wikipedia.org/wiki/FASTA_format
    // A,C,D,G,H,K,M,N,R,S,T,V,W,Y are can be either nucleoide or protein - in particular all of A,C,T,G,N can
    static bool uniq_amino[256]    = { ['E']=true, ['F']=true, ['I']=true, ['L']=true, ['P']=true, ['Q']=true, 
                                       ['X']=true, ['Z']=true,  // may be protein according to the FASTA_format page
                                       ['e']=true, ['f']=true, ['i']=true, ['l']=true, ['p']=true, ['q']=true, 
                                       ['x']=true, ['z']=true };
                                      
    static bool nuke_or_amino[256] = { ['A']=true, ['C']=true, ['D']=true, ['G']=true, ['H']=true, ['K']=true, ['M']=true, 
                                       ['N']=true, ['R']=true, ['S']=true, ['T']=true, ['V']=true, ['W']=true, ['Y']=true,
                                       ['U']=true, ['B']=true, // may be protein according to the FASTA_format page (in addition to standard nuke IUPACs)
                                       ['a']=true, ['c']=true, ['d']=true, ['g']=true, ['h']=true, ['k']=true, ['m']=true, 
                                       ['n']=true, ['r']=true, ['s']=true, ['t']=true, ['v']=true, ['w']=true, ['y']=true,
                                       ['u']=true, ['b']=true };
    
    bool evidence_of_amino=false, evidence_of_both=false;

    for (uint32_t i=0; i < seq_len; i++) {
        if (uniq_amino[(int)seq[i]])    evidence_of_amino = true;
        if (nuke_or_amino[(int)seq[i]]) evidence_of_both = true;

        segconf.seq_type_counter++;
    }

    if (evidence_of_amino) return SQT_AMINO;

    if (evidence_of_both) return segconf.seq_type_counter > 10000 ? SQT_NUKE : SQT_NUKE_OR_AMINO; // A,C,G,T,N are in both - call it as NUKE if we've seen enough without evidence of unique Amino characters
    
    return segconf.seq_type; // unchanged
}

static void fasta_seg_seq_line_do (VBlockFASTA *vb, uint32_t line_len, bool is_first_line_in_contig)
{
    Context *lm_ctx  = CTX(FASTA_LINEMETA);
    Context *seq_ctx = CTX(FASTA_NONREF);

    // line length is same as previous SEQ line
    if (!is_first_line_in_contig && ctx_has_value_in_line_(vb, lm_ctx) && line_len == lm_ctx->last_value.i) 
        seg_duplicate_last (VB, lm_ctx, 0);

    else { 
        char special_snip[100]; unsigned special_snip_len = sizeof (special_snip);
        seg_prepare_snip_other_do (SNIP_OTHER_LOOKUP, (DictId)_FASTA_NONREF, 
                                   true, (int32_t)line_len, &special_snip[3], &special_snip_len);

        special_snip[0] = SNIP_SPECIAL;
        special_snip[1] = FASTA_SPECIAL_SEQ;
        special_snip[2] = '0' + is_first_line_in_contig; 
        seg_by_ctx (VB, special_snip, 3 + special_snip_len, lm_ctx, 0);  // the payload of the special snip, is the OTHER_LOOKUP snip...
    }

    // note: we don't set value for first line, so that seg_duplicate_last doesn't copy it - since special_snip[2] is different 
    if (!is_first_line_in_contig) 
        ctx_set_last_value (VB, lm_ctx, (ValueType){ .i = line_len });

    seq_ctx->txt_len   += line_len;
    seq_ctx->local.len += line_len;
} 

static void fasta_seg_seq_line (VBlockFASTA *vb, STRp(line), 
                                bool is_last_line_vb_no_newline, bool is_last_line_in_contig, 
                                const bool *has_13)
{
    vb->lines_this_contig++;

    *DATA_LINE (vb->line_i) = (ZipDataLineFASTA){ .seq_data_start = ENTNUM (vb->txt_data, line),
                                                  .seq_len        = line_len };

    if (flag.make_reference)
        CTX(FASTA_NONREF)->local.len += line_len;

    // after last line in contig, we know if its SIMPLE or PBWT and seg all lines in this contig into FASTA_LINEMETA
    if (!flag.make_reference && is_last_line_in_contig) {

        ZipDataLineFASTA *dl = DATA_LINE (vb->line_i - vb->lines_this_contig + 1);

        for (int32_t i=0; i < vb->lines_this_contig; i++) {
            fasta_seg_seq_line_do (vb, dl[i].seq_len, i==0);

            // case last Seq line of VB, and this VB ends part-way through the Seq line (no newline)
            if (is_last_line_vb_no_newline && (i == vb->lines_this_contig - 1)) 
                seg_by_did_i (VB, *has_13 ? "\r" : "", *has_13, FASTA_EOL, *has_13); // an EOL without \n
            else 
                SEG_EOL (FASTA_EOL, true); 
        }        
        vb->lines_this_contig = 0;
    }

    if (segconf.running && (segconf.seq_type == SQT_UNKNOWN || segconf.seq_type == SQT_NUKE_OR_AMINO))
        segconf.seq_type = fasta_get_seq_type (STRa(line));

    vb->last_line = FASTA_LINE_SEQ;

    // case: this sequence is continuation from the previous VB - we don't yet know the chrom - we will update it,
    // and increment the min/max_pos relative to the beginning of the seq in the vb later, in random_access_finalize_entries
    if (!vb->chrom_name && !vb->ra_initialized) {
        random_access_update_chrom (VB, DC_PRIMARY, WORD_INDEX_NONE, 0, 0);
        vb->ra_initialized = true;
        vb->chrom_node_index = WORD_INDEX_NONE; // the chrom started in a previous VB, we don't yet know its index
    }

    random_access_increment_last_pos (VB, DC_PRIMARY, line_len); 
}

// Fasta format(s): https://en.wikipedia.org/wiki/FASTA_format
// concept: we segment each line separately, and for each line, we store an element in TEMPLATE about it. The
// Metadata elements are:
// > - description line - this (1) any line starting with > or (2) the first line starting with ; at the start 
//     of a file or after a sequence
//     the descrition line data is further segmented and stored in the DESC dictionary and D0SEC subfields
// ; - a comment line - any other line that starts with a ; or an empty line
//     the comment data (which can be empty for an empty line) is stored in a data buffer (not dictionary)
//     note: if a comment line is the first line in a VB - it will be segmented as a description. No harm done.
// 123 - a sequence line - any line that's not a description of sequence line - store its length
// these ^ are preceded by a 'Y' if the line has a Windows-style \r\n line ending or 'X' if not
const char *fasta_seg_txt_line (VBlockFASTA *vb, const char *line_start, uint32_t remaining_txt_len, bool *has_13) // index in vb->txt_data where this line starts
{
    // get entire line
    unsigned line_len;
    int32_t remaining_vb_txt_len = AFTERENT (char, vb->txt_data) - line_start;
    const char *next_field;
    
    next_field = seg_get_next_line (vb, line_start, &remaining_vb_txt_len, &line_len, !vb->vb_has_no_newline, has_13, "FASTA line");

    // case: description line - we segment it to its components
    if (*line_start == '>' || (*line_start == ';' && vb->last_line == FASTA_LINE_SEQ)) {
        fasta_seg_desc_line (vb, line_start, line_len, has_13);

        if (kraken_is_loaded) {
            unsigned qname_len = strcspn (line_start + 1, " \t\r\n"); // +1 to skip the '>' or ';'
            kraken_seg_taxid (VB, FASTA_TAXID, line_start + 1, qname_len, true);
        }
    }

    // case: comment line - stored in the comment buffer
    else if (*line_start == ';' || !line_len) 
        fast_seg_comment_line (vb, line_start, line_len, has_13);

    // case: sequence line
    else 
        fasta_seg_seq_line (vb, line_start, line_len, 
                            remaining_txt_len == line_len + *has_13, // true if this is the last line in the VB with no newline (but may or may not have \r)
                            !remaining_vb_txt_len || *next_field == '>' || *next_field == ';' || *next_field == '\n' || *next_field == '\r', // is_last_line_in_contig
                            has_13);

    return next_field;
}

//-------------------------
// PIZ stuff
//-------------------------

// Called by thread I/O to initialize for a new genozip file
bool fasta_piz_initialize (void)
{
    fasta_piz_initialize_contig_grepped_out (0,0,0); // initialize 

    return true; // proceed with PIZ
}

// returns true if section is to be skipped reading / uncompressing
bool fasta_piz_is_skip_section (VBlockP vb, SectionType st, DictId dict_id)
{
    if (!vb) return false; // we don't skip reading any SEC_DICT/SEC_COUNTS sections

    if (flag.reading_reference) return false;  // doesn't apply when using FASTA as a reference

    // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
    if (flag.header_only_fast && 
        (dict_id.num == _FASTA_NONREF || dict_id.num == _FASTA_NONREF_X || dict_id.num == _FASTA_COMMENT))
        return true;

    // when grepping by main thread - skipping all sections but DESC
    if ((flag.grep || flag.regions) && (vb->grep_stages == GS_TEST) && 
        dict_id.num != _FASTA_DESC && !dict_id_is_fasta_desc_sf (dict_id))
        return true;

    // if grepping, compute thread doesn't need to decompressed DESC again
    if ((flag.grep || flag.regions) && (vb->grep_stages == GS_UNCOMPRESS) && 
        (dict_id.num == _FASTA_DESC || dict_id_is_fasta_desc_sf (dict_id)))
        return true;

    // no need for the TAXID data if user didn't specify --taxid
    if (flag.kraken_taxid==TAXID_NONE && dict_id.num == _FASTA_TAXID)
        return true;

    return false;
}

// remove trailing newline before SEQ lines in case of --sequential. Note that there might be more than one newline
// in which case the subsequent newlines were part of empty COMMENT lines
static inline void fasta_piz_unreconstruct_trailing_newlines (VBlockFASTA *vb)
{
    char c;
    while ((c = *LASTENT (char, vb->txt_data)) == '\n' || c == '\r') 
        vb->txt_data.len--;

    // update final entries vb->lines to reflect the removal of the final newlines
    for (int32_t line_i = vb->line_i - vb->first_line; 
         line_i >= 0 && *ENT (char *, vb->lines, line_i) > AFTERENT (char, vb->txt_data); 
         line_i--)
        *ENT (char *, vb->lines, line_i) = AFTERENT (char, vb->txt_data);
}

// this is used for end-of-lines of a sequence line, that are not the last line of the sequence. we skip reconstructing
// the newline if the user selected --sequential
SPECIAL_RECONSTRUCTOR (fasta_piz_special_SEQ)
{
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;

    bool is_first_seq_line_in_this_contig = snip[0] - '0';

    // skip showing line if this contig is grepped - but consume it anyway
    if (fasta_vb->contig_grepped_out) vb->drop_curr_line = "grep";

    // --sequential - if this is NOT the first seq line in the contig, we delete the previous end-of-line
    else if (flag.sequential && !is_first_seq_line_in_this_contig) 
        fasta_piz_unreconstruct_trailing_newlines (fasta_vb);

    // in case of not showing the SEQ in the entire file - we can skip consuming it
    if (flag.header_only_fast) // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->drop_curr_line = "header_only_fast";     
    else 
        reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip+1, snip_len-1, true);    

    // case: --sequential, and this seq line is the last line in the vb, and it continues in the next vb
    if (  flag.sequential && // if we are asked for a sequential SEQ
          vb->line_i - vb->first_line == vb->lines.len-1 && // and this is the last line in this vb 
          !vb->drop_curr_line && 
          random_access_does_last_chrom_continue_in_next_vb (vb->vblock_i)) // and this sequence continues in the next VB 
        fasta_piz_unreconstruct_trailing_newlines (fasta_vb); // then: delete final newline if this VB ends with a 

    fasta_vb->last_line = FASTA_LINE_SEQ;

    return false; // no new value
}

SPECIAL_RECONSTRUCTOR (fasta_piz_special_COMMENT)
{
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;

    // skip showing comment line in case cases - but consume it anyway:
    if (  fasta_vb->contig_grepped_out || // 1. if this contig is grepped out
          flag.out_dt == DT_PHYLIP)       // 2. if we're outputting in Phylis format
        vb->drop_curr_line = "grep";

    // in case of not showing the COMMENT in the entire file (--header-only or this is a --reference) - we can skip consuming it
    if (flag.header_only_fast)  // note that flags_update_piz_one_file rewrites --header-only as flag.header_only_fast
        vb->drop_curr_line = "header_only_fast";     
    else 
        reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, true);    

    fasta_vb->last_line = FASTA_LINE_COMMENT;

    return false; // no new value
}

// this is called by piz_test_grep - it is called sequentially for all VBs by the main thread
// returns true if the last contig of the previous VB was grepped-in
bool fasta_piz_initialize_contig_grepped_out (VBlock *vb_, bool does_vb_have_any_desc, bool last_desc_in_this_vb_matches_grep)
{
    VBlockFASTA *vb = (VBlockFASTA *)vb_;

    // we pass the info from one VB to the next using this static variable
    static bool prev_vb_last_contig_grepped_out = false; 
    static uint32_t prev_vb_i = 0;

    if (!vb) { // vb=0 means initialize
        prev_vb_last_contig_grepped_out = false;
        prev_vb_i = 0;
        return 0;
    }

    // we're continuing the contig in the previous VB - until DESC is encountered
    vb->contig_grepped_out = prev_vb_last_contig_grepped_out || // last contig of previous VB had last_desc_in_this_vb_matches_grep
                             (prev_vb_i + 1 < vb->vblock_i);    // previous VB was skipped in piz_one_txt_file due to random_access_is_vb_included
    
    // update for use of next VB, IF this VB contains any DESC line, otherwise just carry forward the current value
    if (does_vb_have_any_desc) 
        prev_vb_last_contig_grepped_out = !last_desc_in_this_vb_matches_grep; 

    prev_vb_i = vb->vblock_i;
    
    return !vb->contig_grepped_out;
}

// Phylip format mandates exact 10 space-padded characters: http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html
static inline void fasta_piz_translate_desc_to_phylip (VBlock *vb, char *desc_start)
{
    uint32_t recon_len = AFTERENT (const char, vb->txt_data) - desc_start;
    *AFTERENT (char, vb->txt_data) = 0; // nul-terminate

    const char *chrom_name = desc_start + 1;
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");

    memmove (desc_start, chrom_name, MIN_(chrom_name_len, 10));
    if (chrom_name_len < 10) memcpy (desc_start + chrom_name_len, "          ", 10-chrom_name_len); // pad with spaces

    if (recon_len > 10) vb->txt_data.len -= recon_len - 10; // we do it this way to avoid signed problems
    else                vb->txt_data.len += 10 - recon_len;
}

// shorten DESC to the first white space
static inline void fasta_piz_desc_header_one (VBlock *vb, char *desc_start)
{
    uint32_t recon_len = AFTERENT (const char, vb->txt_data) - desc_start;
    *AFTERENT (char, vb->txt_data) = 0; // nul-terminate
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n");
    
    vb->txt_data.len -= recon_len - chrom_name_len -1;
}

SPECIAL_RECONSTRUCTOR (fasta_piz_special_DESC)
{
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;
    fasta_vb->contig_grepped_out = false;

    char *desc_start = AFTERENT (char, vb->txt_data);
    reconstruct_one_snip (vb, ctx, WORD_INDEX_NONE, snip, snip_len, true);    
    *AFTERENT (char, vb->txt_data) = 0; // for strstr and strcspn

    // if --grep: here we decide whether to show this contig or not
    if (flag.grep) 
        fasta_vb->contig_grepped_out = !strstr (desc_start, flag.grep);
    
    unsigned chrom_name_len = strcspn (desc_start + 1, " \t\r\n"); // +1 to skip the '>'
    
    // --taxid: grep out by Kraken taxid 
    if (flag.kraken_taxid) 
        fasta_vb->contig_grepped_out |= 
            (!kraken_is_loaded && !kraken_is_included_stored (vb, FASTA_TAXID, false)) ||
            ( kraken_is_loaded && !kraken_is_included_loaded (vb, desc_start + 1, chrom_name_len));
    
    vb->chrom_node_index = vb->last_index(CHROM) = ctx_search_for_word_index (CTX(CHROM), desc_start + 1, chrom_name_len);

    // note: this logic allows the to grep contigs even if --no-header 
    if (fasta_vb->contig_grepped_out)
        fasta_vb->drop_curr_line = "grep";     

    else if (flag.no_header)
        fasta_vb->drop_curr_line = "no_header";     

    if (flag.out_dt == DT_PHYLIP) 
        fasta_piz_translate_desc_to_phylip (vb, desc_start);

    if (flag.header_one)
        fasta_piz_desc_header_one (vb, desc_start);

    fasta_vb->last_line = FASTA_LINE_DESC;

    return false; // no new value
}

bool fasta_piz_read_one_vb (VBlock *vb, Section sl)
{ 
    // if we're grepping we we uncompress and reconstruct the DESC from the main thread, and terminate here if this VB is to be skipped
    if ((flag.grep || flag.regions) && !piz_test_grep (vb)) return false; 

    return true;
}

//---------------------------------------
// Multifasta -> PHYLIP translation stuff
//---------------------------------------

// create Phylip header line
TXTHEADER_TRANSLATOR (txtheader_fa2phy)
{
    // get length of contigs and error if they are not all the same length
    uint32_t contig_len = random_access_verify_all_contigs_same_length();

    bufprintf (comp_vb, txtheader_buf, "%"PRIu64" %u\n", ZCTX(FASTA_CONTIG)->word_list.len, contig_len);
}

// Translating FASTA->PHYLIP: drop EOL after DESC + don't allow multiple consecutive EOL
TRANSLATOR_FUNC (fasta_piz_fa2phy_EOL)
{
    VBlockFASTA *fasta_vb = (VBlockFASTA *)vb;

    if (fasta_vb->last_line == FASTA_LINE_DESC // line is a DESC
    ||  recon == vb->txt_data.data     // initial EOL - not allowed in Phylip)
    ||  recon[-1] == '\n')             // previous item was an EOL - remove this one then
    
        vb->txt_data.len -= recon_len;
    
    return 0;
}

CONTAINER_FILTER_FUNC (fasta_piz_filter)
{
    // currently unused, but specified in FASTA toplevel containers created up to v11, so the function needs to exist
    return true;
}

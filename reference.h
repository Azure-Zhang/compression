// ------------------------------------------------------------------
//   reference.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef REFERENCE_INCLUDED
#define REFERENCE_INCLUDED

#include <pthread.h>
#include "genozip.h"
#include "buffer.h"
#include "md5.h"
#include "bit_array.h"

// reference sequences - one per range of 1MB. ranges (chrom, pos) are mapped here with a hash function. In the rare case two unrelated ranges
// are mapped to the same entry - only the first range will have a reference, and the other ones will not. this will hurt the compression ratio,
// but not the correctness.
// Thread safety: "ref" is set atomically as the last set in initialization. If its set, then it is correct and will never change. If it appears to be
// not set yet, it is necesarry to lock the mutex and test again, and initialize if still not set.

typedef struct Range {
    BitArray ref;                // actual reference data - 2-bit array
    BitArray is_set;             // a 1-bit array - SEG: a pos is set if seg set this reference PIZ: is set if SEC_REF_IS_SET said so
    const char *chrom_name;
    unsigned chrom_name_len;
    WordIndex chrom;             // index to the chrom of the external reference
    uint32_t range_i;
    PosType first_pos, last_pos; // the range that includes all locii (note: in ZIP-INTERNAL it might include unset locii too)
    PosType gpos;                // position of this range in the "global position" 
    uint32_t copied_first_index, copied_len; // ZIP with REF_EXT_STORE: the subset of this range that was copied directly from the fasta file and doesn't need to be compressed
} Range;

#define ref_size(r) ((r) ? ((r)->last_pos - (r)->first_pos + 1) : 0)

#define REF_NUM_DENOVO_RANGES (1 << 20)
#define REF_NUM_DENOVO_SITES_PER_RANGE (1 << 20) // 1 Mbp

// ZIP ONLY: access range_i and index within range, for ranges configured for ZIP
#define range_i2pos(range_i) ((PosType)(range_i) * REF_NUM_DENOVO_SITES_PER_RANGE)
#define pos2range_i(pos)   ((uint32_t)((pos) / REF_NUM_DENOVO_SITES_PER_RANGE))
#define pos2range_idx(pos) ((uint32_t)((pos) % REF_NUM_DENOVO_SITES_PER_RANGE))

// locks
typedef struct { int32_t first_mutex, last_mutex; } RefLock;
#define REFLOCK_NONE ((RefLock){-1,-1})

extern RefLock ref_lock (PosType gpos_start, uint32_t seq_len);
extern RefLock ref_unlock (RefLock lock);
extern RefLock ref_lock_range (int32_t range_id);

typedef enum { RT_NONE,     // value of ranges.param if ranges is unallocated
               RT_MAKE_REF, // used in --make-ref one range per vb of fasta reference file - ranges in order of the fasta file
               RT_DENOVO,   // used in SAM with REF_INTERNAL - an large array of Range's, hashed by the chrom and pos
               RT_LOADED    // one Range per chrom (contig), overlayed on genome
              } RangesType;
              
extern void ref_initialize_ranges (RangesType ranges_type);
extern void ref_compress_ref (void);
extern void ref_load_stored_reference (void);
extern void ref_set_reference (const char *filename);
extern void ref_set_ref_file_info (Md5Hash md5, const char *fasta_name);
extern void ref_load_external_reference (bool display, bool is_last_file);
extern void ref_unload_reference (bool force_clean_all);
extern MemStats ref_memory_consumption (void);
extern const Range *ref_piz_get_range (VBlockP vb, PosType first_pos_needed, uint32_t num_nucleotides_needed);
extern void ref_consume_ref_fasta_global_area (void);
extern Range *ref_seg_get_locked_range (VBlockP vb, PosType pos, uint32_t seq_len, const char *field /* used for ASSSEG */, RefLock *lock);
extern void ref_print_subrange (const char *msg, const Range *r, PosType start_pos, PosType end_pos);
extern void ref_print_is_set (const Range *r, PosType around_pos);
extern const char *ref_get_cram_ref (void);
extern void ref_output_vb (VBlockP vb);
extern void ref_make_ref_init (void);
extern void ref_generate_reverse_complement_genome (void);

// contigs stuff
extern void ref_contigs_get (ConstBufferP *out_contig_dict, ConstBufferP *out_contigs);
extern void ref_contigs_verify_identical_chrom (const char *chrom_name, unsigned chrom_name_len, PosType last_pos);
extern void ref_contigs_sort_chroms (void);
extern void ref_contigs_load_contigs (void);

// alt chroms stuff
extern void ref_alt_chroms_load (void);
extern void ref_alt_chroms_compress (void);

#define ref_assert_nucleotide_available(range,pos) \
    ASSERT (/* piz w stored ref */ (flag_reference == REF_STORED && ref_is_nucleotide_set ((range), pos2range_idx(pos))) ||  \
            /* zip w ext ref    */ ((flag_reference == REF_EXTERNAL || flag_reference == REF_EXT_STORE) && ((pos) >= (range)->first_pos && (pos) <= (range)->last_pos)) || \
            /* zip internal ref */ flag_reference == REF_INTERNAL, \
        "Error in %s:%u: reference is not set: chrom=%.*s pos=%"PRId64, __FUNCTION__, __LINE__, (range)->chrom_name_len, (range)->chrom_name, (pos))

// note that the following work on idx and not pos! (idx is the index within the range)
#define ref_set_nucleotide(range,idx,value) { bit_array_assign (&(range)->ref, (idx) * 2,      acgt_encode[(uint8_t)value] & 1)       ;  \
                                              bit_array_assign (&(range)->ref, (idx) * 2 + 1, (acgt_encode[(uint8_t)value] & 2) >> 1) ; }

#define ref_is_nucleotide_set(range,idx) ((bool)bit_array_get (&(range)->is_set, (idx)))

#define ref_get_nucleotide(range,idx)   acgt_decode[(bit_array_get (&(range)->ref, (idx) * 2 + 1) << 1) | \
                                                     bit_array_get (&(range)->ref, (idx) * 2)]

// globals
extern const char *ref_filename;
extern Md5Hash ref_md5;
extern Buffer ref_stored_ra;

extern Range genome, genome_rev;
extern PosType genome_size;

#endif

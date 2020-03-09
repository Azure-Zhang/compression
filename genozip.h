// ------------------------------------------------------------------
//   genozip.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef GENOZIP_INCLUDED
#define GENOZIP_INCLUDED

#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER // Microsoft compiler
#include <inttypes.h>
#include <stdbool.h>
#include <unistd.h>
#else
#include "compatability/visual_c_stdint.h"
#include "compatability/visual_c_stdbool.h"
#include "compatability/visual_c_unistd.h"
#include "compatability/visual_c_misc_funcs.h"
#endif

#define GENOZIP_EXT ".genozip"

#define DEFAULT_MAX_MEMORY_PER_VB (128*1024*1024)

#define MAX_PLOIDY         100  // set to a reasonable 100 to avoid memory allocation explosion in case of an error in the VCF file
#if MAX_PLOIDY > 65535
#error "MAX_PLOIDY cannot go beyond 65535 as SectionHeaderVbHeader.ploidy and VariantBlock.ploidy are uint16_t"
#endif

#define MAX_SUBFIELDS      63   // maximum number of FORMAT subfield types (except for GT) and INFO subfield types that is supported in one GENOZIP file.

#define MAX_DICTS          (MAX_SUBFIELDS + MAX_SUBFIELDS + 8)   // dictionaries of subfields, infos and the 9 first fields (tabs) of the VCF file (+8 because REF and ALT are combined). 
#if MAX_DICTS > 255
#error "MAX_DICTS cannot go beyond 255 as SubfieldMapperZip and SubfieldInfoMapperPiz represent did_i as uint_8, and NIL=255"
#endif

// Note: the algorithm will use as many cores as it can - but there's no speed penalty for a higher MAX_COMPUTE_THREADS
// that number of cores - it will not be faster, but also not slower.
// However, each thread adds memory consumption approximately linearly

// the optimal number of compute threads is determined by the ratio between CPU time and I/O time. 
// For uncompress, 3 or more threads result in similar performance for HD and SSD, with 7 seaming to be about optimal (but not a big difference than 3). 
// Adding threads doesn't help. 2 or less threads result in significantly slower execution time. 
// Memory consumption is linear with the number of threads (each allocated a VB)

#define DEFAULT_MAX_THREADS 8 // maximum compute threads created - one I/O thread, and the rest of compute threads. This can be changed with command line option -@
 
#define MAX_32BIT_WINDOWS_MEMORY (1.7*1024*1024*1024) // 1.7GB - so Windows 32bit code doesn't explode at 2GB. TO DO - make this platform specific or get ulimits


// pointers used in header files - so we don't need to include the whole .h (and avoid cyclicity)
typedef struct variant_block_ *VariantBlockP;
typedef const struct variant_block_ *ConstVariantBlockP;
typedef struct file_ *FileP;
typedef const struct file_ *ConstFileP;
typedef struct buffer_ *BufferP;
typedef const struct buffer_ *ConstBufferP;
typedef struct mtfcontext_ *MtfContextP;
typedef struct mtfnode_ *MtfNodeP;

typedef enum { EXE_GENOZIP, EXE_GENOUNZIP, EXE_GENOLS, EXE_GENOCAT } ExeType;

// global parameters - set before any thread is created, and never change
extern uint32_t    global_num_samples, global_number_displayed_samples, global_max_memory_per_vb;
extern const char *global_cmd;            // set once in main()
extern int command;
extern ExeType exe_type;

// flags set by user's command line options
extern int flag_force, flag_quiet, flag_concat, flag_md5, flag_split, flag_show_alleles, flag_show_time, 
           flag_show_memory, flag_show_dict, flag_show_gt_nodes, flag_show_b250, flag_show_sections, flag_show_headers,
           flag_show_index, flag_show_gheader, flag_stdout, flag_replace, flag_show_content, flag_test, flag_regions,
           flag_samples, flag_drop_genotypes, flag_no_header, flag_header_only, flag_show_threads, flag_debug_memory;

// external vb - used when an operation is needed outside of the context of a specific variant block;
extern VariantBlockP evb;

// macros
#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b) )
#define MAX(a, b) (((a) > (b)) ? (a) : (b) )
#endif

// sanity checks
static inline void my_exit(void) { exit(1); }// an exit function so we can put a debugging break point when ASSERT exits
#define ASSERT(condition, format, ...)       { if (!(condition)) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); my_exit(); }}
#define ASSERT0(condition, string)           { if (!(condition)) { fprintf (stderr, "\n%s\n", string); my_exit(); }}
#define ASSERTW(condition, format, ...)      { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); }}
#define ASSERTW0(condition, string)          { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n%s\n", string); } }
#define RETURNW(condition, ret, format, ...) { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); return ret; }}
#define RETURNW0(condition, ret, string)     { if (!(condition) && !flag_quiet) { fprintf (stderr, "\n%s\n", string); return ret; } }
#define ABORT(format, ...)                   { fprintf (stderr, "\n"); fprintf (stderr, format, __VA_ARGS__); fprintf (stderr, "\n"); my_exit();}
#define ABORT0(string)                       { fprintf (stderr, "\n%s\n", string); my_exit();}

#endif
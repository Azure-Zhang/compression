// ------------------------------------------------------------------
//   bgzf.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
 
#include "genozip.h"

#define BGZF_MAX_BLOCK_SIZE 65536 // maximum block size of both compressed and uncompressed data of one block

// BGZF EOF marker is simply an empty block, see https://samtools.github.io/hts-specs/SAMv1.pdf section 4.1.2
#define BGZF_EOF_LEN 28
#define BGZF_EOF "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

// data type of vblock.bgzf_blocks
typedef struct BgzfBlockZip {
    uint32_t txt_index, txt_size;    // index of uncompressed block within vb->txt_data. The first block doesn't necessarily have index=0 bc there could be passed-down data
    uint32_t compressed_index, comp_size; // index within vb->compressed
    bool is_decompressed;
} BgzfBlockZip;

typedef struct BgzfBlockPiz {
    int32_t txt_index, txt_size; // index of uncompressed block within vb->txt_data. The first block index will be negative if there is passed-down unconsumed data
} BgzfBlockPiz;

//---------
// ZIP side
//---------

#define BGZF_BLOCK_GZIP_NOT_BGZIP -1
#define BGZF_BLOCK_IS_NOT_GZIP    -2
extern int32_t bgzf_read_block (FileP file, uint8_t *block, uint32_t *block_size, bool soft_fail);
extern void bgzf_uncompress_vb (VBlockP vb);
extern void bgzf_uncompress_one_block (VBlockP vb, BgzfBlockZip *bb);
extern void bgzf_compress_bgzf_section (void);

//---------
// PIZ side
//---------

extern bool bgzf_load_isizes (ConstSectionListEntryP sl_ent);
extern void bgzf_calculate_blocks_one_vb (VBlockP vb, uint32_t vb_txt_data_len);
extern void bgzf_compress_vb (VBlockP vb);
extern void bgzf_write_to_disk (VBlockP vb);


// ------------------------------------------------------------------
//   compresssor.h
//   Copyright (C) 2019-2022 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

#include "genozip.h"
#include "codec.h"

extern uint32_t comp_compress (VBlockP vb, BufferP z_data, SectionHeaderP header, 
                               rom uncompressed_data, // option 1 - compress contiguous data
                               LocalGetLineCB callback, rom name); // option 2 - compress data one line at a time

extern bool comp_compress_complex_codec (VBlockP vb, ContextP ctx, SectionHeaderP header, bool is_2nd_try,
                                         uint32_t *uncompressed_len, STRe (compressed), rom name);

extern void comp_uncompress (VBlockP vb, Codec codec, Codec sub_codec, uint8_t param,
                             STRp(compressed_data), Buffer *uncompressed_data, uint64_t uncompressed_len, rom name);

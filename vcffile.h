// ------------------------------------------------------------------
//   vcffile.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef VCFFILE_INCLUDED
#define VCFFILE_INCLUDED

#include "genozip.h"

//extern bool vcffile_get_line (VariantBlockP vb, unsigned line_i_in_file, bool skip_md5_vcf_header, BufferP line, const char *buf_name);

extern void vcffile_read_vcf_header (bool is_first_vcf);
extern void vcffile_read_variant_block (VariantBlockP vb);
extern void vcffile_write_one_variant_block (VariantBlockP vb);
extern unsigned vcffile_write_to_disk (ConstBufferP buf);

#endif
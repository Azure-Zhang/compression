// ------------------------------------------------------------------
//   stats.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#pragma once

extern void stats_compress (void); // ZIP
extern void stats_display (void);  // PIZ and ZIP
extern void stats_read_and_display (void);     // PIZ
extern void stats_add_txt_name (const char *fn);
extern void stats_set_consolidation (VBlockP vb, DidIType parent, unsigned num_deps, ...);
extern void stats_set_consolidation_(VBlockP vb, DidIType parent, unsigned num_deps, ContextP *dep_ctxs);

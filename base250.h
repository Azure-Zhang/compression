// ------------------------------------------------------------------
//   base250.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef BASE250_INCLUDED
#define BASE250_INCLUDED

#include "genozip.h"

// values 0 to 249 are used as numerals in base-250. 
// The remaining 6 values below are control characters, and can only appear in numerals[0].
#define BASE250_EMPTY_SF   250 // subfield declared in FORMAT is empty, terminating : present
#define BASE250_MISSING_SF 251 // subfield declared in FORMAT is missing at end of cell, no :
#define BASE250_ONE_UP     252 // value is one higher than previous value. used in B250_ENC_8
#define BASE250_2_NUMERALS 253 // this number has 2 numerals, starting from numerals[1]. Used in BASE250_ENCODING_V1
#define BASE250_MOST_FREQ  253 // this translates to 0 representing the most frequent value (according to vb_i=1 sorting). used in B250_ENC_16
#define BASE250_3_NUMERALS 254 // this number has 3 numerals. 
#define BASE250_4_NUMERALS 255 // this number has 4 numerals. 
#define MAX_BASE250_NUMERALS 5
typedef struct {
    uint32_t n;                                  // the number bering encoding
    uint8_t  numerals[2][MAX_BASE250_NUMERALS];  // encoded number. up to 5 numerals (first array for 8bit, 2nd for 16 bit encoding)
    uint8_t  num_numerals[2];                    // legal values - 1,2,3,4 (8 bit and 16 bit)
} Base250;

// these values go into the SectionHeaderBase250
typedef enum { B250_ENC_NONE=-1, B250_ENC_8=0, B250_ENC_16=1 } Base250Encoding; // B250_ENC_8/16 are used as indeces in arrays, so they need to be 0/1

// B250_ENC_8:  if n <= 249: one numeral which is n
//                         if n >= 250, first numeral is a code BASE250_2_NUMERALS, BASE250_3_NUMERALS or BASE250_4_NUMERALS
//                         and the next 2 to 4 numerals are the number in base-250 (each numeral 0 to 249), least significant first

// B250_ENC_16: if n==0 (representing the most frequent snip): one numeral which is BASE250_MOST_FREQ
//                         if 1 <= n < 250*250: 2 numerals, least significant first
//                         if n >= 250*250: first numeral is a code BASE250_3_NUMERALS or BASE250_4_NUMERALS
//                         and the next 3 or 4 numerals are the number in base-250 (each numeral 0 to 249), least significant first

extern Base250 base250_encode (uint32_t n);
extern uint32_t base250_decode (const uint8_t **str_p, Base250Encoding encoding); // decodes and advances str_p

//#define base250_len(data) (*(data) < BASE250_2_NUMERALS ? 1 : *(data) - BASE250_2_NUMERALS + 3) // number of bytes this base250 number consumes
// new econding:
unsigned base250_len (const uint8_t *data, Base250Encoding encoding);

#endif
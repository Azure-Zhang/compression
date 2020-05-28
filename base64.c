// ------------------------------------------------------------------
//   base64.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt
//
// inspired by: https://github.com/launchdarkly/c-client-sdk/blob/master/base64.c
//

#include "base64.h"

static const uint8_t base64_table[65] = 
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static const uint8_t dtable[256] = {
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 0-15
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 16-31
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 62  , 0x80, 0x80, 0x80, 63  , // ASCII 32-47
	52  , 53  , 54  , 55  , 56  , 57  , 58  , 59  , 60  , 61  , 0x80, 0x80, 0x80, 0   , 0x80, 0x80, // ASCII 48-63 (outside of base64_table '=' = 0) 
    0x80, 0   , 1   , 2   , 3   , 4   , 5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13  , 14  , // ASCII 64-79
	15  , 16  , 17  , 18  , 19  , 20  , 21  , 22  , 23  , 24  , 25  , 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 80-95
	0x80, 26  , 27  , 28  , 29  , 30  , 31  , 32  , 33  , 34  , 35  , 36  , 37  , 38  , 39  , 40  , // ASCII 96-111
	41  , 42  , 43  , 44  , 45  , 46  , 47  , 48  , 49  , 50  , 51  , 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 112-127
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, // ASCII 128-255
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
	0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80 };

// returns length of encoded (which is at most base64_sizeof)
// out must be allocated base64_sizeof bytes
unsigned base64_encode (const uint8_t *in, unsigned in_len, char *b64_str)
{
    ASSERT0 (in, "Error in base64_encode: in is NULL");

	const uint8_t *end = in + in_len;
	char *next = b64_str;
	while (end - in >= 3) {
		*next++ = base64_table[in[0] >> 2];
		*next++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
		*next++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
		*next++ = base64_table[  in[2] & 0x3f];
		in += 3;
	}

	if (end - in) {
		*next++ = base64_table[in[0] >> 2];
		if (end - in == 1) {
			*next++ = base64_table[(in[0] & 0x03) << 4];
			*next++ = '=';
		} else {
			*next++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
			*next++ = base64_table[ (in[1] & 0x0f) << 2];
		}
		*next++ = '=';
	}

    return next - b64_str;
}


void base64_decode (const char *b64_str, unsigned b64_str_len, uint8_t *out, unsigned *out_len /* out */)
{
	ASSERT (b64_str_len && !(b64_str_len % 4), "Error in base64_decode: bad base64 - expecting it to be a string with length divisable by 4 but its length is %u: %.*s",
            b64_str_len, b64_str_len, b64_str);

	unsigned pad=0;
	uint8_t block[4];
	for (unsigned i=0; i < b64_str_len; i++) {

		if (b64_str[i] == '=') pad++;
		block[i&3] = dtable[(unsigned)b64_str[i]];
		ASSERT (block[i&3] != 0x80, "Invalid character '%c' found in b64 string: %.*s", b64_str[i], b64_str_len, b64_str);

		if ((i&3) == 3) {
  			                 *out++ = (block[0] << 2) | (block[1] >> 4);
			if (pad <= 1)    *out++ = (block[1] << 4) | (block[2] >> 2);
			if (!pad)        *out++ = (block[2] << 6) |  block[3];
		}
	}	

	*out_len = (b64_str_len / 4) * 3 - pad;
}
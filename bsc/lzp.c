/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Lempel Ziv Prediction                                     */
/*-----------------------------------------------------------*/

/*--

This file is a part of bsc and/or libbsc, a program and a library for
lossless, block-sorting data compression.

   Copyright (c) 2009-2012 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information and file AUTHORS
for full list of contributors.

See also the bsc and libbsc web site:
  http://libbsc.com/ for more information.

--*/

// ------------------------------------------------------------------
//   This file was extensively modified to adapt it to genozip. All modifications:
//   Copyright (C) 2020-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "lzp.h"
#include "platform.h"
#include "libbsc.h"

#define LIBBSC_LZP_MATCH_FLAG 	0xf2

static inline int bsc_lzp_num_blocks(int n)
{
    if (n <       256 * 1024)   return 1;
    if (n <  4 * 1024 * 1024)   return 2;
    if (n < 16 * 1024 * 1024)   return 4;

    return 8;
}

static int bsc_lzp_encode_block (void *vb, const unsigned char * input, const unsigned char * inputEnd, unsigned char * output, unsigned char * outputEnd, int hashSize, int minLen)
{
    if (inputEnd - input < 16)
    {
        return LIBBSC_NOT_COMPRESSIBLE;
    }

    int * lookup = (int *)bsc_zero_malloc (vb, (int)(1 << hashSize) * sizeof(int));
    if (lookup)
    {
        unsigned int            mask        = (int)(1 << hashSize) - 1;
        const unsigned char *   inputStart  = input;
        const unsigned char *   outputStart = output;
        const unsigned char *   outputEOB   = outputEnd - 4;

        unsigned int context = 0;
        for (int i = 0; i < 4; ++i)
        {
            context = (context << 8) | (*output++ = *input++);
        }

        const unsigned char * heuristic      = input;
        const unsigned char * inputMinLenEnd = inputEnd - minLen - 8;
        while ((input < inputMinLenEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; 
            lookup[index] = (int)(input - inputStart);
            if (value > 0)
            {
                const unsigned char * reference = inputStart + value;
                if ((*(unsigned int *)(input + minLen - 4) == *(unsigned int *)(reference + minLen - 4)) && (*(unsigned int *)(input) == *(unsigned int *)(reference)))
                {
                    if ((heuristic > input) && (*(unsigned int *)heuristic != *(unsigned int *)(reference + (heuristic - input))))
                    {
                        goto LIBBSC_LZP_MATCH_NOT_FOUND;
                    }

                    int len = 4;
                    for (; input + len < inputMinLenEnd; len += 4)
                    {
                        if (*(unsigned int *)(input + len) != *(unsigned int *)(reference + len)) break;
                    }
                    if (len < minLen)
                    {
                        if (heuristic < input + len) heuristic = input + len;
                        goto LIBBSC_LZP_MATCH_NOT_FOUND;
                    }

                    if (input[len] == reference[len]) len++;
                    if (input[len] == reference[len]) len++;
                    if (input[len] == reference[len]) len++;

                    input += len; context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

                    *output++ = LIBBSC_LZP_MATCH_FLAG;

                    len -= minLen; while (len >= 254) { len -= 254; *output++ = 254; if (output >= outputEOB) break; }

                    *output++ = (unsigned char)(len);
                }
                else
                {
LIBBSC_LZP_MATCH_NOT_FOUND:{
                        unsigned char next = *output++ = *input++;  context = (context << 8) | next;
                        if (next == LIBBSC_LZP_MATCH_FLAG) *output++ = 255;
                    }
                }
            }
            else
            {
                context = (context << 8) | (*output++ = *input++);
            }
        }

        while ((input < inputEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);
            if (value > 0)
            {
                unsigned char next = *output++ = *input++; context = (context << 8) | next;
                if (next == LIBBSC_LZP_MATCH_FLAG) *output++ = 255;
            }
            else
            {
                context = (context << 8) | (*output++ = *input++);
            }
        }

        bsc_free (vb, lookup);

        return (output >= outputEOB) ? LIBBSC_NOT_COMPRESSIBLE : (int)(output - outputStart);
    }

    return LIBBSC_NOT_ENOUGH_MEMORY;
}

static int bsc_lzp_decode_block (void *vb, const unsigned char * input, const unsigned char * inputEnd, unsigned char * output, int hashSize, int minLen)
{
    if (inputEnd - input < 4)
    {
        return LIBBSC_UNEXPECTED_EOB;
    }

    int * lookup = (int *)bsc_zero_malloc (vb, (int)(1 << hashSize) * sizeof(int));
    if (lookup)
    {
        unsigned int            mask        = (int)(1 << hashSize) - 1;
        const unsigned char *   outputStart = output;

        unsigned int context = 0;
        for (int i = 0; i < 4; ++i)
        {
            context = (context << 8) | (*output++ = *input++);
        }

        while (input < inputEnd)
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(output - outputStart);
            if (*input == LIBBSC_LZP_MATCH_FLAG && value > 0)
            {
                input++;
                if (*input != 255)
                {
                    int len = minLen; while (1) { len += *input; if (*input++ != 254) break; }

                    const unsigned char * reference = outputStart + value;
                          unsigned char * outputEnd = output + len;

                    while (output < outputEnd) *output++ = *reference++;

                    context = output[-1] | (output[-2] << 8) | (output[-3] << 16) | (output[-4] << 24);
                }
                else
                {
                    input++; context = (context << 8) | (*output++ = LIBBSC_LZP_MATCH_FLAG);
                }
            }
            else
            {
                context = (context << 8) | (*output++ = *input++);
            }
        }

        bsc_free (vb, lookup);

        return (int)(output - outputStart);
    }

    return LIBBSC_NOT_ENOUGH_MEMORY;
}

static int bsc_lzp_compress_serial(void *vb, const unsigned char * input, unsigned char * output, int n, int hashSize, int minLen)
{
    if (bsc_lzp_num_blocks(n) == 1)
    {
        int result = bsc_lzp_encode_block(vb, input, input + n, output + 1, output + n - 1, hashSize, minLen);
        if (result >= LIBBSC_NO_ERROR) result = (output[0] = 1, result + 1);

        return result;
    }

    int nBlocks   = bsc_lzp_num_blocks(n);
    int chunkSize = n / nBlocks;
    int outputPtr = 1 + 8 * nBlocks;

    output[0] = nBlocks;
    for (int blockId = 0; blockId < nBlocks; ++blockId)
    {
        int inputStart  = blockId * chunkSize;
        int inputSize   = blockId != nBlocks - 1 ? chunkSize : n - inputStart;
        int outputSize  = inputSize; if (outputSize > n - outputPtr) outputSize = n - outputPtr;

        int result = bsc_lzp_encode_block(vb, input + inputStart, input + inputStart + inputSize, output + outputPtr, output + outputPtr + outputSize, hashSize, minLen);
        if (result < LIBBSC_NO_ERROR)
        {
            if (outputPtr + inputSize >= n) return LIBBSC_NOT_COMPRESSIBLE;
            result = inputSize; memcpy(output + outputPtr, input + inputStart, inputSize);
        }

        *(int *)(output + 1 + 8 * blockId + 0) = inputSize;
        *(int *)(output + 1 + 8 * blockId + 4) = result;

        outputPtr += result;
    }

    return outputPtr;
}


int bsc_lzp_compress(void *vb, const unsigned char * input, unsigned char * output, int n, int hashSize, int minLen, int features)
{
    return bsc_lzp_compress_serial (vb, input, output, n, hashSize, minLen);
}

int bsc_lzp_decompress(void *vb, const unsigned char * input, unsigned char * output, int n, int hashSize, int minLen, int features)
{
    int nBlocks = input[0];

    if (nBlocks == 1)
    {
        return bsc_lzp_decode_block(vb, input + 1, input + n, output, hashSize, minLen);
    }

    int decompressionResult[ALPHABET_SIZE];
    {
        for (int blockId = 0; blockId < nBlocks; ++blockId)
        {
            int inputPtr = 0;  for (int p = 0; p < blockId; ++p) inputPtr  += *(int *)(input + 1 + 8 * p + 4);
            int outputPtr = 0; for (int p = 0; p < blockId; ++p) outputPtr += *(int *)(input + 1 + 8 * p + 0);

            inputPtr += 1 + 8 * nBlocks;

            int inputSize  = *(int *)(input + 1 + 8 * blockId + 4);
            int outputSize = *(int *)(input + 1 + 8 * blockId + 0);

            if (inputSize != outputSize)
            {
                decompressionResult[blockId] = bsc_lzp_decode_block(vb, input + inputPtr, input + inputPtr + inputSize, output + outputPtr, hashSize, minLen);
            }
            else
            {
                decompressionResult[blockId] = inputSize; memcpy(output + outputPtr, input + inputPtr, inputSize);
            }
        }
    }

    int dataSize = 0, result = LIBBSC_NO_ERROR;
    for (int blockId = 0; blockId < nBlocks; ++blockId)
    {
        if (decompressionResult[blockId] < LIBBSC_NO_ERROR) result = decompressionResult[blockId];
        dataSize += decompressionResult[blockId];
    }

    return (result == LIBBSC_NO_ERROR) ? dataSize : result;
}

/*-----------------------------------------------------------*/
/* End                                               lzp.cpp */
/*-----------------------------------------------------------*/

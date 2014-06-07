///
/// @file   popcount.cpp
/// @brief  Fast algorithms to count the number of 1 bits in an array.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#if !defined(__STDC_CONSTANT_MACROS)
  #define __STDC_CONSTANT_MACROS
#endif

#include <stdint.h>

namespace {

const uint64_t start_masks[128] =
{
  UINT64_C(0xffffffffffffffff), UINT64_C(0xffffffffffffffff), 
  UINT64_C(0xfffffffffffffffe), UINT64_C(0xfffffffffffffffe), 
  UINT64_C(0xfffffffffffffffc), UINT64_C(0xfffffffffffffffc), 
  UINT64_C(0xfffffffffffffff8), UINT64_C(0xfffffffffffffff8), 
  UINT64_C(0xfffffffffffffff0), UINT64_C(0xfffffffffffffff0), 
  UINT64_C(0xffffffffffffffe0), UINT64_C(0xffffffffffffffe0), 
  UINT64_C(0xffffffffffffffc0), UINT64_C(0xffffffffffffffc0), 
  UINT64_C(0xffffffffffffff80), UINT64_C(0xffffffffffffff80), 
  UINT64_C(0xffffffffffffff00), UINT64_C(0xffffffffffffff00), 
  UINT64_C(0xfffffffffffffe00), UINT64_C(0xfffffffffffffe00), 
  UINT64_C(0xfffffffffffffc00), UINT64_C(0xfffffffffffffc00), 
  UINT64_C(0xfffffffffffff800), UINT64_C(0xfffffffffffff800), 
  UINT64_C(0xfffffffffffff000), UINT64_C(0xfffffffffffff000), 
  UINT64_C(0xffffffffffffe000), UINT64_C(0xffffffffffffe000), 
  UINT64_C(0xffffffffffffc000), UINT64_C(0xffffffffffffc000), 
  UINT64_C(0xffffffffffff8000), UINT64_C(0xffffffffffff8000), 
  UINT64_C(0xffffffffffff0000), UINT64_C(0xffffffffffff0000), 
  UINT64_C(0xfffffffffffe0000), UINT64_C(0xfffffffffffe0000), 
  UINT64_C(0xfffffffffffc0000), UINT64_C(0xfffffffffffc0000), 
  UINT64_C(0xfffffffffff80000), UINT64_C(0xfffffffffff80000), 
  UINT64_C(0xfffffffffff00000), UINT64_C(0xfffffffffff00000), 
  UINT64_C(0xffffffffffe00000), UINT64_C(0xffffffffffe00000), 
  UINT64_C(0xffffffffffc00000), UINT64_C(0xffffffffffc00000), 
  UINT64_C(0xffffffffff800000), UINT64_C(0xffffffffff800000), 
  UINT64_C(0xffffffffff000000), UINT64_C(0xffffffffff000000), 
  UINT64_C(0xfffffffffe000000), UINT64_C(0xfffffffffe000000), 
  UINT64_C(0xfffffffffc000000), UINT64_C(0xfffffffffc000000), 
  UINT64_C(0xfffffffff8000000), UINT64_C(0xfffffffff8000000), 
  UINT64_C(0xfffffffff0000000), UINT64_C(0xfffffffff0000000), 
  UINT64_C(0xffffffffe0000000), UINT64_C(0xffffffffe0000000), 
  UINT64_C(0xffffffffc0000000), UINT64_C(0xffffffffc0000000), 
  UINT64_C(0xffffffff80000000), UINT64_C(0xffffffff80000000), 
  UINT64_C(0xffffffff00000000), UINT64_C(0xffffffff00000000), 
  UINT64_C(0xfffffffe00000000), UINT64_C(0xfffffffe00000000), 
  UINT64_C(0xfffffffc00000000), UINT64_C(0xfffffffc00000000), 
  UINT64_C(0xfffffff800000000), UINT64_C(0xfffffff800000000), 
  UINT64_C(0xfffffff000000000), UINT64_C(0xfffffff000000000), 
  UINT64_C(0xffffffe000000000), UINT64_C(0xffffffe000000000), 
  UINT64_C(0xffffffc000000000), UINT64_C(0xffffffc000000000), 
  UINT64_C(0xffffff8000000000), UINT64_C(0xffffff8000000000), 
  UINT64_C(0xffffff0000000000), UINT64_C(0xffffff0000000000), 
  UINT64_C(0xfffffe0000000000), UINT64_C(0xfffffe0000000000), 
  UINT64_C(0xfffffc0000000000), UINT64_C(0xfffffc0000000000), 
  UINT64_C(0xfffff80000000000), UINT64_C(0xfffff80000000000), 
  UINT64_C(0xfffff00000000000), UINT64_C(0xfffff00000000000), 
  UINT64_C(0xffffe00000000000), UINT64_C(0xffffe00000000000), 
  UINT64_C(0xffffc00000000000), UINT64_C(0xffffc00000000000), 
  UINT64_C(0xffff800000000000), UINT64_C(0xffff800000000000), 
  UINT64_C(0xffff000000000000), UINT64_C(0xffff000000000000), 
  UINT64_C(0xfffe000000000000), UINT64_C(0xfffe000000000000), 
  UINT64_C(0xfffc000000000000), UINT64_C(0xfffc000000000000), 
  UINT64_C(0xfff8000000000000), UINT64_C(0xfff8000000000000), 
  UINT64_C(0xfff0000000000000), UINT64_C(0xfff0000000000000), 
  UINT64_C(0xffe0000000000000), UINT64_C(0xffe0000000000000), 
  UINT64_C(0xffc0000000000000), UINT64_C(0xffc0000000000000), 
  UINT64_C(0xff80000000000000), UINT64_C(0xff80000000000000), 
  UINT64_C(0xff00000000000000), UINT64_C(0xff00000000000000), 
  UINT64_C(0xfe00000000000000), UINT64_C(0xfe00000000000000), 
  UINT64_C(0xfc00000000000000), UINT64_C(0xfc00000000000000), 
  UINT64_C(0xf800000000000000), UINT64_C(0xf800000000000000), 
  UINT64_C(0xf000000000000000), UINT64_C(0xf000000000000000), 
  UINT64_C(0xe000000000000000), UINT64_C(0xe000000000000000), 
  UINT64_C(0xc000000000000000), UINT64_C(0xc000000000000000), 
  UINT64_C(0x8000000000000000), UINT64_C(0x8000000000000000)
};

const uint64_t stop_masks[128] =
{
  UINT64_C(0x1),                UINT64_C(0x1), 
  UINT64_C(0x3),                UINT64_C(0x3), 
  UINT64_C(0x7),                UINT64_C(0x7), 
  UINT64_C(0xf),                UINT64_C(0xf), 
  UINT64_C(0x1f),               UINT64_C(0x1f), 
  UINT64_C(0x3f),               UINT64_C(0x3f), 
  UINT64_C(0x7f),               UINT64_C(0x7f), 
  UINT64_C(0xff),               UINT64_C(0xff), 
  UINT64_C(0x1ff),              UINT64_C(0x1ff), 
  UINT64_C(0x3ff),              UINT64_C(0x3ff), 
  UINT64_C(0x7ff),              UINT64_C(0x7ff), 
  UINT64_C(0xfff),              UINT64_C(0xfff), 
  UINT64_C(0x1fff),             UINT64_C(0x1fff), 
  UINT64_C(0x3fff),             UINT64_C(0x3fff), 
  UINT64_C(0x7fff),             UINT64_C(0x7fff), 
  UINT64_C(0xffff),             UINT64_C(0xffff), 
  UINT64_C(0x1ffff),            UINT64_C(0x1ffff), 
  UINT64_C(0x3ffff),            UINT64_C(0x3ffff), 
  UINT64_C(0x7ffff),            UINT64_C(0x7ffff), 
  UINT64_C(0xfffff),            UINT64_C(0xfffff), 
  UINT64_C(0x1fffff),           UINT64_C(0x1fffff), 
  UINT64_C(0x3fffff),           UINT64_C(0x3fffff), 
  UINT64_C(0x7fffff),           UINT64_C(0x7fffff), 
  UINT64_C(0xffffff),           UINT64_C(0xffffff), 
  UINT64_C(0x1ffffff),          UINT64_C(0x1ffffff), 
  UINT64_C(0x3ffffff),          UINT64_C(0x3ffffff), 
  UINT64_C(0x7ffffff),          UINT64_C(0x7ffffff), 
  UINT64_C(0xfffffff),          UINT64_C(0xfffffff), 
  UINT64_C(0x1fffffff),         UINT64_C(0x1fffffff), 
  UINT64_C(0x3fffffff),         UINT64_C(0x3fffffff), 
  UINT64_C(0x7fffffff),         UINT64_C(0x7fffffff), 
  UINT64_C(0xffffffff),         UINT64_C(0xffffffff), 
  UINT64_C(0x1ffffffff),        UINT64_C(0x1ffffffff), 
  UINT64_C(0x3ffffffff),        UINT64_C(0x3ffffffff), 
  UINT64_C(0x7ffffffff),        UINT64_C(0x7ffffffff), 
  UINT64_C(0xfffffffff),        UINT64_C(0xfffffffff), 
  UINT64_C(0x1fffffffff),       UINT64_C(0x1fffffffff), 
  UINT64_C(0x3fffffffff),       UINT64_C(0x3fffffffff), 
  UINT64_C(0x7fffffffff),       UINT64_C(0x7fffffffff), 
  UINT64_C(0xffffffffff),       UINT64_C(0xffffffffff), 
  UINT64_C(0x1ffffffffff),      UINT64_C(0x1ffffffffff), 
  UINT64_C(0x3ffffffffff),      UINT64_C(0x3ffffffffff), 
  UINT64_C(0x7ffffffffff),      UINT64_C(0x7ffffffffff), 
  UINT64_C(0xfffffffffff),      UINT64_C(0xfffffffffff), 
  UINT64_C(0x1fffffffffff),     UINT64_C(0x1fffffffffff), 
  UINT64_C(0x3fffffffffff),     UINT64_C(0x3fffffffffff), 
  UINT64_C(0x7fffffffffff),     UINT64_C(0x7fffffffffff), 
  UINT64_C(0xffffffffffff),     UINT64_C(0xffffffffffff), 
  UINT64_C(0x1ffffffffffff),    UINT64_C(0x1ffffffffffff), 
  UINT64_C(0x3ffffffffffff),    UINT64_C(0x3ffffffffffff), 
  UINT64_C(0x7ffffffffffff),    UINT64_C(0x7ffffffffffff), 
  UINT64_C(0xfffffffffffff),    UINT64_C(0xfffffffffffff), 
  UINT64_C(0x1fffffffffffff),   UINT64_C(0x1fffffffffffff), 
  UINT64_C(0x3fffffffffffff),   UINT64_C(0x3fffffffffffff), 
  UINT64_C(0x7fffffffffffff),   UINT64_C(0x7fffffffffffff), 
  UINT64_C(0xffffffffffffff),   UINT64_C(0xffffffffffffff), 
  UINT64_C(0x1ffffffffffffff),  UINT64_C(0x1ffffffffffffff), 
  UINT64_C(0x3ffffffffffffff),  UINT64_C(0x3ffffffffffffff), 
  UINT64_C(0x7ffffffffffffff),  UINT64_C(0x7ffffffffffffff), 
  UINT64_C(0xfffffffffffffff),  UINT64_C(0xfffffffffffffff), 
  UINT64_C(0x1fffffffffffffff), UINT64_C(0x1fffffffffffffff), 
  UINT64_C(0x3fffffffffffffff), UINT64_C(0x3fffffffffffffff), 
  UINT64_C(0x7fffffffffffffff), UINT64_C(0x7fffffffffffffff), 
  UINT64_C(0xffffffffffffffff), UINT64_C(0xffffffffffffffff)
};

const uint64_t m1  = UINT64_C(0x5555555555555555);
const uint64_t m2  = UINT64_C(0x3333333333333333);
const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
const uint64_t h01 = UINT64_C(0x0101010101010101);

int64_t popcount(uint64_t bits)
{
  uint64_t bit_count = 0;
  uint64_t x;

  // http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
  x = bits;
  x =  x       - ((x >> 1)  & m1);
  x = (x & m2) + ((x >> 2)  & m2);
  x = (x       +  (x >> 4)) & m4;
  x = (x * h01) >> 56;
  bit_count += x;

  return bit_count;
}

int64_t popcount(const uint64_t* bits, int64_t start_idx, int64_t stop_idx)
{
  uint64_t bit_count = 0;
  uint64_t x;

  // http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
  for (int64_t i = start_idx; i <= stop_idx; i++)
  {
    x = bits[i];
    x =  x       - ((x >> 1)  & m1);
    x = (x & m2) + ((x >> 2)  & m2);
    x = (x       +  (x >> 4)) & m4;
    x = (x * h01) >> 56;
    bit_count += x;
  }

  return bit_count;
}

int64_t popcount_edges(const uint64_t* bits, int64_t start, int64_t stop)
{
  int64_t count = 0;
  int64_t start_idx = start >> 7;
  int64_t stop_idx = stop >> 7;
  int64_t start_mask = start_masks[start & 127];
  int64_t stop_mask = stop_masks[stop & 127];

  if (start_idx == stop_idx)
    count += popcount(bits[start_idx] & (start_mask & stop_mask));
  else
  {
    count += popcount(bits[start_idx] & start_mask);
    count += popcount(bits[stop_idx] & stop_mask);
  }

  return count;
}

} // namespace

namespace primecount {

int64_t popcount(const uint64_t* bits, int64_t start, int64_t stop, int64_t low)
{
  start += ~(start + (low & 1)) & 1;
  stop  -= ~(stop  + (low & 1)) & 1;

  if (start > stop)
    return 0;

  int64_t count = popcount_edges(bits, start, stop);

  int64_t start_idx = (start >> 7) + 1;
  int64_t stop_idx  = (stop  >> 7) - 1;

  if (start_idx <= stop_idx)
    count += ::popcount(bits, start_idx, stop_idx);

  return count;
}

} // namespace

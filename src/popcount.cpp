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

int64_t popcount(uint64_t bits)
{
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t h01 = UINT64_C(0x0101010101010101);

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
  const uint64_t m1  = UINT64_C(0x5555555555555555);
  const uint64_t m2  = UINT64_C(0x3333333333333333);
  const uint64_t m4  = UINT64_C(0x0F0F0F0F0F0F0F0F);
  const uint64_t h01 = UINT64_C(0x0101010101010101);

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
  int64_t start_idx = start >> 6;
  int64_t stop_idx = stop >> 6;
  int64_t start_mask = UINT64_C(0xffffffffffffffff) << (start & 63);
  int64_t stop_mask  = UINT64_C(0xffffffffffffffff) >> (63 - (stop & 63));

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

int64_t popcount(const uint64_t* bits, int64_t start, int64_t stop)
{
  if (start > stop)
    return 0;

  int64_t count = popcount_edges(bits, start, stop);

  int64_t start_idx = (start >> 6) + 1;
  int64_t stop_idx  = (stop  >> 6) - 1;

  if (start_idx <= stop_idx)
    count += ::popcount(bits, start_idx, stop_idx);

  return count;
}

} // namespace

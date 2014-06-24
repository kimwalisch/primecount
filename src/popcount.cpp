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
#include <popcount64.hpp>

namespace {

int64_t popcount_edges(const uint64_t* bits, int64_t start, int64_t stop)
{
  int64_t index1 = start >> 6;
  int64_t index2 = stop  >> 6;
  int64_t mask1 = UINT64_C(0xffffffffffffffff) << (start & 63);
  int64_t mask2 = UINT64_C(0xffffffffffffffff) >> (63 - (stop & 63));
  int64_t count = 0;

  if (index1 == index2)
    count += popcount64(bits[index1] & (mask1 & mask2));
  else
  {
    count += popcount64(bits[index1] & mask1);
    count += popcount64(bits[index2] & mask2);
  }

  return count;
}

int64_t popcount(const uint64_t* bits, int64_t start_idx, int64_t stop_idx)
{
  uint64_t bit_count = 0;

  for (int64_t i = start_idx; i <= stop_idx; i++)
    bit_count += popcount64(bits[i]);

  return bit_count;
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

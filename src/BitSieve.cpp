///
/// @file  BitSieve.cpp
/// @brief The BitSieve class is a bit array for use with
///        Eratosthenes-like prime sieving algorithms. BitSieve
///        assigns 64 numbers to the bits of an 8 byte word thus
///        reducing the memory usage by a factor of 8.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#if !defined(__STDC_CONSTANT_MACROS)
  #define __STDC_CONSTANT_MACROS
#endif

#include <BitSieve.hpp>
#include <popcount.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <vector>

namespace primecount {

const uint64_t BitSieve::unset_bit_[64] =
{
  ~(UINT64_C(1) <<  0), ~(UINT64_C(1) <<  1), ~(UINT64_C(1) <<  2),
  ~(UINT64_C(1) <<  3), ~(UINT64_C(1) <<  4), ~(UINT64_C(1) <<  5),
  ~(UINT64_C(1) <<  6), ~(UINT64_C(1) <<  7), ~(UINT64_C(1) <<  8),
  ~(UINT64_C(1) <<  9), ~(UINT64_C(1) << 10), ~(UINT64_C(1) << 11),
  ~(UINT64_C(1) << 12), ~(UINT64_C(1) << 13), ~(UINT64_C(1) << 14),
  ~(UINT64_C(1) << 15), ~(UINT64_C(1) << 16), ~(UINT64_C(1) << 17),
  ~(UINT64_C(1) << 18), ~(UINT64_C(1) << 19), ~(UINT64_C(1) << 20),
  ~(UINT64_C(1) << 21), ~(UINT64_C(1) << 22), ~(UINT64_C(1) << 23),
  ~(UINT64_C(1) << 24), ~(UINT64_C(1) << 25), ~(UINT64_C(1) << 26),
  ~(UINT64_C(1) << 27), ~(UINT64_C(1) << 28), ~(UINT64_C(1) << 29),
  ~(UINT64_C(1) << 30), ~(UINT64_C(1) << 31), ~(UINT64_C(1) << 32),
  ~(UINT64_C(1) << 33), ~(UINT64_C(1) << 34), ~(UINT64_C(1) << 35),
  ~(UINT64_C(1) << 36), ~(UINT64_C(1) << 37), ~(UINT64_C(1) << 38),
  ~(UINT64_C(1) << 39), ~(UINT64_C(1) << 40), ~(UINT64_C(1) << 41),
  ~(UINT64_C(1) << 42), ~(UINT64_C(1) << 43), ~(UINT64_C(1) << 44),
  ~(UINT64_C(1) << 45), ~(UINT64_C(1) << 46), ~(UINT64_C(1) << 47),
  ~(UINT64_C(1) << 48), ~(UINT64_C(1) << 49), ~(UINT64_C(1) << 50),
  ~(UINT64_C(1) << 51), ~(UINT64_C(1) << 52), ~(UINT64_C(1) << 53),
  ~(UINT64_C(1) << 54), ~(UINT64_C(1) << 55), ~(UINT64_C(1) << 56),
  ~(UINT64_C(1) << 57), ~(UINT64_C(1) << 58), ~(UINT64_C(1) << 59),
  ~(UINT64_C(1) << 60), ~(UINT64_C(1) << 61), ~(UINT64_C(1) << 62),
  ~(UINT64_C(1) << 63)
};

BitSieve::BitSieve(std::size_t size) :
  bits_(ceil_div(size, 64)),
  size_(size)
{ }

/// Set all bits to 1, except bits corresponding
/// to 0, 1 and even numbers > 2.
/// 
void BitSieve::fill(uint64_t low, uint64_t /* unused */)
{
  if (low % 2 == 0)
    std::fill(bits_.begin(), bits_.end(), UINT64_C(0xAAAAAAAAAAAAAAAA));
  else
    std::fill(bits_.begin(), bits_.end(), UINT64_C(0x5555555555555555));

  // correct 0, 1 and 2
  if (low <= 2)
  {
    uint64_t bit = 1 << (2 - low);
    bits_[0] &= ~(bit - 1);
    bits_[0] |= bit;
  }
}

/// Count the number of 1 bits inside [start, stop]
uint64_t BitSieve::count(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  assert(stop < size_);

  uint64_t start_idx = start / 64;
  uint64_t stop_idx = stop / 64;
  uint64_t m1 = UINT64_C(0xffffffffffffffff) << (start % 64);
  uint64_t m2 = UINT64_C(0xffffffffffffffff) >> (63 - stop % 64);
  uint64_t bit_count;

  if (start_idx == stop_idx)
  {
    uint64_t bits = bits_[start_idx] & (m1 & m2);
    bit_count = popcount64(bits);
  }
  else
  {
    uint64_t bits = bits_[start_idx] & m1;
    bit_count = popcount64(bits);

    for (uint64_t i = start_idx + 1; i < stop_idx; i++)
      bit_count += popcount64(bits_[i]);

    bits = bits_[stop_idx] & m2;
    bit_count += popcount64(bits);
  }

  return bit_count;
}

} // namespace

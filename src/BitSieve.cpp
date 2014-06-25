///
/// @file  BitSieve.cpp
/// @brief Bit array for prime sieving.
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
#include <algorithm>
#include <cassert>
#include <vector>

#include <BitSieve.hpp>
#include <popcount64.hpp>

namespace primecount {

const unsigned int BitSieve::unset_bit_[32] =
{
  ~(1u <<  0), ~(1u <<  1), ~(1u <<  2),
  ~(1u <<  3), ~(1u <<  4), ~(1u <<  5),
  ~(1u <<  6), ~(1u <<  7), ~(1u <<  8),
  ~(1u <<  9), ~(1u << 10), ~(1u << 11),
  ~(1u << 12), ~(1u << 13), ~(1u << 14),
  ~(1u << 15), ~(1u << 16), ~(1u << 17),
  ~(1u << 18), ~(1u << 19), ~(1u << 20),
  ~(1u << 21), ~(1u << 22), ~(1u << 23),
  ~(1u << 24), ~(1u << 25), ~(1u << 26),
  ~(1u << 27), ~(1u << 28), ~(1u << 29),
  ~(1u << 30), ~(1u << 31)
};

BitSieve::BitSieve(std::size_t size) :
  size_(size)
{
  std::size_t size32 = (size + 31) / 32;
  // align to 64-bit boundary
  if (size32 % 2 != 0)
    size32 += 2 - size32 % 2;
  bits_.resize(size32);
}

/// Set all bits to 1, except bits corresponding
/// to 0, 1 and even numbers > 2.
/// 
void BitSieve::memset(uint64_t low)
{
  unsigned mask = (low & 1) ? 0x55555555u : 0xAAAAAAAAu;
  std::fill(bits_.begin(), bits_.end(), mask);

  // correct 0, 1 and 2
  if (low <= 2)
  {
    uint32_t bitmask = ~((1u << (2 - low)) - 1);
    bits_[0] &= bitmask;
    bits_[0] |= 1 << (2 - low);
  }
}

uint64_t BitSieve::popcount(const uint64_t* bits, uint64_t start, uint64_t stop)
{
  if (start > stop)
    return 0;

  uint64_t bit_count = popcount_edges(bits, start, stop);
  uint64_t start_idx = (start >> 6) + 1;
  uint64_t limit = stop  >> 6;

  for (uint64_t i = start_idx; i < limit; i++)
    bit_count += popcount64(bits[i]);

  return bit_count;
}

uint64_t BitSieve::popcount_edges(const uint64_t* bits, uint64_t start, uint64_t stop)
{
  uint64_t index1 = start >> 6;
  uint64_t index2 = stop  >> 6;
  uint64_t mask1 = UINT64_C(0xffffffffffffffff) << (start & 63);
  uint64_t mask2 = UINT64_C(0xffffffffffffffff) >> (63 - (stop & 63));
  uint64_t count = 0;

  if (index1 == index2)
    count += popcount64(bits[index1] & (mask1 & mask2));
  else
  {
    count += popcount64(bits[index1] & mask1);
    count += popcount64(bits[index2] & mask2);
  }

  return count;
}

/// Count the number of 1 bits inside [start, stop]
uint64_t BitSieve::count(uint64_t start, uint64_t stop) const
{
  assert(stop < size_);
  const uint64_t* bits = reinterpret_cast<const uint64_t*>(bits_.data());
  return popcount(bits, start, stop);
}

} // namespace

///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n < x^(1/2). Since a PrimePi[n] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[n]
///        lookup table that uses only O(x^(1/4)) memory.
///
///        The SegmentedPiTable class is a compressed lookup table of
///        prime counts. Each bit of the lookup table corresponds to
///        an integer that is not divisible by 2, 3 and 5. The 8 bits
///        of each byte correspond to the offsets { 1, 7, 11, 13, 17,
///        19, 23, 29 }. Since our lookup table uses the uint64_t data
///        type, one array element (8 bytes) corresponds to an
///        interval of size 30 * 8 = 240.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <min.hpp>

#include <stdint.h>
#include <cassert>

namespace primecount {

void SegmentedPiTable::init(uint64_t low, uint64_t high)
{
  assert(low < high);
  assert(low % 240 == 0);
  uint64_t pi_low;

  if (low <= 5)
    pi_low = pi_tiny_[5];
  else if (low == high_)
    pi_low = operator[](low - 1);
  else
    pi_low = pi_noprint(low - 1, 1);

  low_ = low;
  high_ = high;
  uint64_t segment_size = high - low;
  uint64_t size = ceil_div(segment_size, 240);
  pi_.clear();
  pi_.resize(size);

  init_bits();
  init_count(pi_low);
}

/// Each thread computes PrimePi [low, high[
void SegmentedPiTable::init_bits()
{
  // Iterate over primes > 5
  primesieve::iterator it(max(low_, 5), high_);
  uint64_t prime = 0;

  // Each thread iterates over the primes
  // inside [low, high[ and initializes
  // the pi[x] lookup table.
  while ((prime = it.next_prime()) < high_)
  {
    uint64_t p = prime - low_;
    pi_[p / 240].bits |= set_bit_[p % 240];
  }
}

/// Each thread computes PrimePi [low, high[
void SegmentedPiTable::init_count(uint64_t pi_low)
{
  uint64_t high_idx = ceil_div((high_ - low_), 240);

  for (uint64_t i = 0; i < high_idx; i++)
  {
    pi_[i].count = pi_low;
    pi_low += popcnt64(pi_[i].bits);
  }
}

int64_t SegmentedPiTable::get_segment_size(uint64_t max_high, uint64_t x13, int threads)
{
  // CPU cache sizes per core
  uint64_t l1_cache_size = 32 << 10;
  uint64_t l2_cache_size = 512 << 10;
  uint64_t l3_cache_size = 16 << 20;
  uint64_t numbers_per_byte = 240 / sizeof(pi_t);
  uint64_t segment_size;

  if (threads == 1)
    segment_size = l2_cache_size * numbers_per_byte;
  else
  {
    // Segment size is tiny: x^(1/4)
    segment_size = isqrt(max_high);

    // Minimum segment size = 1 KiB
    uint64_t min_size = (1 << 10) * numbers_per_byte;
    // Most special leaves are below x^(1/3)
    uint64_t max_size = x13 / threads;
    max_size = max(segment_size, max_size);

    // Using a segment_size > x^(1/4) reduces thread scheduling
    // overhead. But we must ensure that the segment_size is
    // not larger than the CPU's cache (per core).
         if ((segment_size * 8) / sizeof(pi_t) <= l1_cache_size)
      segment_size *= 8;
    else if ((segment_size * 4) / sizeof(pi_t) <= l2_cache_size)
      segment_size *= 4;
    else if ((segment_size * 2) / sizeof(pi_t) <= l3_cache_size)
      segment_size *= 2;

    segment_size = in_between(min_size, segment_size, max_size);
  }

  segment_size = min(segment_size, max_high);
  segment_size += 240 - segment_size % 240;
  assert(segment_size % 240 == 0);

  return segment_size;
}

} // namespace

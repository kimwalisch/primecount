///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n < x^(1/2). Since a PrimePi[n] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[n]
///        lookup table that uses only O(y) memory.
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
#include <cstring>

namespace primecount {

SegmentedPiTable::SegmentedPiTable(uint64_t low,
                                   uint64_t segment_size)
{
  if (low <= 5)
    pi_low_ = pi_tiny_[5];
  else
    pi_low_ = pi_simple(low - 1, 1);

  segment_size_ = segment_size;
  pi_.resize(segment_size_ / 240);
  low_ = low;
  high_ = low + segment_size;
  init_bits(low_, high_);
  init_count(low_, high_);
}

/// Each thread computes PrimePi [start, stop[
void SegmentedPiTable::init_bits(uint64_t start, uint64_t stop)
{
  // Iterate over primes > 5
  primesieve::iterator it(max(start, 5), stop);
  uint64_t prime = 0;

  // Each thread iterates over the primes
  // inside [start, stop[ and initializes
  // the pi[x] lookup table.
  while ((prime = it.next_prime()) < stop)
  {
    uint64_t p = prime - low_;
    pi_[p / 240].bits |= set_bit_[p % 240];
  }
}

/// Each thread computes PrimePi [start, stop[
void SegmentedPiTable::init_count(uint64_t start, uint64_t stop)
{
  // First compute PrimePi[start - 1]
  uint64_t count = pi_low_;

  // Convert to array indexes
  uint64_t i = (start - low_) / 240;
  uint64_t stop_idx = ceil_div((stop - low_), 240);

  for (; i < stop_idx; i++)
  {
    pi_[i].count = count;
    count += popcnt64(pi_[i].bits);
  }
}

} // namespace

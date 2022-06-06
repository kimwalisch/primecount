///
/// @file  SegmentedPiTable.cpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[x] values with
///        x < x^(1/2). Since a PrimePi[x] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[x]
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
///        The algorithm of the easy special leaves and the usage of
///        the SegmentedPiTable are described in more detail in:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.md
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {

void SegmentedPiTable::init(uint64_t low, uint64_t high)
{
  ASSERT(low < high);
  ASSERT(low % 240 == 0);
  int threads = 1;
  uint64_t pi_low;

  // In order to make the threads completely independent
  // from each other each thread needs compute PrimePi[low]
  // at the start of each new segment. If only a single
  // thread is used we can compute PrimePi[low] in O(1) by
  // getting that value from the previous segment.
  if (low <= 5)
    pi_low = pi_tiny_[5];
  else if (low == high_)
    pi_low = operator[](low - 1);
  else
    pi_low = pi_noprint(low - 1, threads);

  low_ = low;
  high_ = high;
  uint64_t segment_size = high - low;
  uint64_t size = ceil_div(segment_size, 240);

  pi_.clear();
  pi_.resize(size);
  std::fill(pi_.begin(), pi_.end(), pi_t{0, 0});

  init_bits();
  init_count(pi_low);
}

/// Each thread computes PrimePi [low, high[
void SegmentedPiTable::init_bits()
{
  // Iterate over primes > 5
  uint64_t low = max(low_, 5);
  if (low + 1 >= high_)
    return;

  primesieve::iterator it(low, high_);
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
  uint64_t j = ceil_div((high_ - low_), 240);

  for (uint64_t i = 0; i < j; i++)
  {
    pi_[i].count = pi_low;
    pi_low += popcnt64(pi_[i].bits);
  }
}

} // namespace

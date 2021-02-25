///
/// @file  SegmentedPiTable.hpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n < x^(1/2). Since a PrimePi[n] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[n]
///        lookup table that uses only O(z) memory.
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

#ifndef SEGMENTEDPITABLE_HPP
#define SEGMENTEDPITABLE_HPP

#include <BitSieve240.hpp>
#include <popcnt.hpp>
#include <aligned_vector.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <cassert>

namespace primecount {

class SegmentedPiTable : public BitSieve240
{
public:
  SegmentedPiTable(uint64_t low,
                   uint64_t limit,
                   uint64_t segment_size,
                   int threads);

  void next();

  bool finished() const
  {
    return low_ >= max_high_;
  }

  int64_t low() const
  {
    return low_;
  }

  int64_t high() const
  {
    return high_;
  }

  /// Get number of primes <= n
  int64_t operator[](uint64_t n) const
  {
    assert(n >= low_);
    assert(n < high_);

    if_unlikely(n < pi_tiny_.size())
      return pi_tiny_[n];

    n -= low_;
    uint64_t bitmask = unset_larger_[n % 240];
    uint64_t prime_count = pi_[n / 240].prime_count;
    uint64_t bit_count = popcnt64(pi_[n / 240].bits & bitmask);
    return prime_count + bit_count;
  }

private:
  void init_bits(uint64_t start, uint64_t stop, uint64_t thread_num);
  void init_prime_count(uint64_t start, uint64_t stop, uint64_t thread_num);
  void init();

  struct pi_t
  {
    uint64_t prime_count;
    uint64_t bits;
  };

  pod_vector<pi_t> pi_;
  aligned_vector<uint64_t> counts_;
  uint64_t low_ = 0;
  uint64_t pi_low_ = 3;
  uint64_t high_;
  uint64_t max_high_;
  uint64_t segment_size_;
  int threads_;
};

} // namespace

#endif

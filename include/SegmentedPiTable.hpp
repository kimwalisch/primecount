///
/// @file  SegmentedPiTable.hpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[n] values with
///        n < x^(1/2). Since a PrimePi[n] lookup table of size
///        x^(1/2) would use too much memory we need a segmented
///        PrimePi[n] lookup table that uses only O(z) memory.
///
///        The SegmentedPiTable is based on the PiTable class which
///        is a compressed lookup table for prime counts. Each bit
///        in the lookup table corresponds to an odd integer and that
///        bit is set to 1 if the integer is a prime. PiTable uses
///        only (n / 8) bytes of memory and returns the number of
///        primes <= n in O(1) operations.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SEGMENTEDPITABLE_HPP
#define SEGMENTEDPITABLE_HPP

#include <popcnt.hpp>
#include <aligned_vector.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <array>
#include <cassert>

namespace primecount {

class SegmentedPiTable
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

    // Since we store only odd numbers in our lookup table,
    // we cannot store 2 which is the only even prime.
    // As a workaround we mark 1 as a prime (1st bit) and
    // add a check to return 0 for pi[1].
    if_unlikely(n == 1)
      return 0;

    n -= low_;
    uint64_t bitmask = unset_bits_[n % 128];
    uint64_t prime_count = pi_[n / 128].prime_count;
    uint64_t bit_count = popcnt64(pi_[n / 128].bits & bitmask);
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

  static const std::array<uint64_t, 128> unset_bits_;
  pod_vector<pi_t> pi_;
  aligned_vector<uint64_t> counts_;
  uint64_t low_ = 0;
  uint64_t pi_low_ = 0;
  uint64_t high_;
  uint64_t max_high_;
  uint64_t segment_size_;
  int threads_;
};

} // namespace

#endif

///
/// @file  SegmentedPiTable.hpp
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

#ifndef SEGMENTEDPITABLE_HPP
#define SEGMENTEDPITABLE_HPP

#include <BitSieve240.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {

class SegmentedPiTable : public BitSieve240
{
public:
  void init(uint64_t low, uint64_t high);

  int64_t low() const
  {
    return low_;
  }

  int64_t high() const
  {
    return high_;
  }

  static constexpr int64_t numbers_per_byte()
  {
    return 240 / sizeof(pi_t);
  }

  /// Make sure size % 240 == 0
  static int64_t get_segment_size(uint64_t size)
  {
    size = std::max<uint64_t>(240, size);

    if (size % 240)
      size += 240 - size % 240;

    return size;
  }

  /// Get number of primes <= x
  ALWAYS_INLINE int64_t operator[](uint64_t x) const
  {
    ASSERT(x >= low_);
    ASSERT(x < high_);

    if_unlikely(x < pi_tiny_.size())
      return pi_tiny_[x];

    x -= low_;
    uint64_t count = pi_[x / 240].count;
    uint64_t bits = pi_[x / 240].bits;
    uint64_t bitmask = unset_larger_[x % 240];
    return count + popcnt64(bits & bitmask);
  }

private:
  void init_bits();
  void init_count(uint64_t pi_low);

  struct pi_t
  {
    uint64_t count;
    uint64_t bits;
  };

  pod_vector<pi_t> pi_;
  uint64_t low_ = 0;
  uint64_t high_ = 0;
};

} // namespace

#endif

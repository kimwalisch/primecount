///
/// @file  SegmentedPiTable.hpp
/// @brief The A and C formulas in Xavier Gourdon's prime counting
///        algorithm require looking up PrimePi[x] values with
///        x < x^(1/2). Since a PrimePi[x] lookup table of size x^(1/2)
///        would use too much memory we need a segmented PrimePi[x]
///        lookup table that uses only O(x^(1/4)) memory.
///
///        The SegmentedPiTable class is a compressed lookup table of
///        prime counts. Since the size of SegmentedPiTable is very
///        small and will always fit into the CPU's cache, we don't
///        use a bit array with maximum compression because this adds
///        significant overhead. Instead we use a bit array where each
///        bit corresponds to an odd integer. This compression scheme
///        provides very fast access since the bit array index can be
///        calculated using a single right shift instruction.
///
///        The algorithm of the easy special leaves and the usage of
///        the SegmentedPiTable are described in more detail in:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SEGMENTEDPITABLE_HPP
#define SEGMENTEDPITABLE_HPP

#include <macros.hpp>
#include <Vector.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {

class SegmentedPiTable
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
    return 128 / sizeof(pi_t);
  }

  /// Make sure size % 128 == 0
  static int64_t align_segment_size(uint64_t size)
  {
    size = std::max<uint64_t>(128, size);

    if (size % 128)
      size += 128 - size % 128;

    return size;
  }

  /// Get number of primes <= x
  ALWAYS_INLINE int64_t operator[](uint64_t x) const
  {
    ASSERT(x >= low_);
    ASSERT(x < high_);

    // Workaround needed for prime 2 since
    // we are sieving with primes >= 3.
    if (x < 2)
      return 0;

    x -= low_;
    uint64_t count = pi_[x / 128].count;
    uint64_t bits = pi_[x / 128].bits;
    uint64_t bitmask = unset_larger_[x % 128];
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

  static const Array<uint64_t, 128> unset_larger_;
  Vector<pi_t> pi_;
  uint64_t low_ = 0;
  uint64_t high_ = 0;
};

} // namespace

#endif

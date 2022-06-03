///
/// @file  PiTable.hpp
/// @brief The PiTable class is a compressed lookup table of prime
///        counts. Each bit of the lookup table corresponds to an
///        integer that is not divisible by 2, 3 and 5. The 8 bits of
///        each byte correspond to the offsets { 1, 7, 11, 13, 17, 19,
///        23, 29 }. Since our lookup table uses the uint64_t data
///        type, one array element (8 bytes) corresponds to an
///        interval of size 30 * 8 = 240.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PITABLE_HPP
#define PITABLE_HPP

#include <BitSieve240.hpp>
#include <popcnt.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <cassert>

namespace primecount {

class PiTable : public BitSieve240
{
public:
  PiTable(uint64_t max_x, int threads);

  uint64_t size() const
  {
    return max_x_ + 1;
  }

  static int64_t max_cached()
  {
    return pi_cache_.size() * 240 - 1;
  }

  /// Get number of primes <= x
  ALWAYS_INLINE int64_t operator[](uint64_t x) const
  {
    assert(x <= max_x_);

    if_unlikely(x < pi_tiny_.size())
      return pi_tiny_[x];

    uint64_t count = pi_[x / 240].count;
    uint64_t bits = pi_[x / 240].bits;
    uint64_t bitmask = unset_larger_[x % 240];
    return count + popcnt64(bits & bitmask);
  }

  /// Get number of primes <= x
  static int64_t pi_cache(uint64_t x)
  {
    if_unlikely(x < pi_tiny_.size())
      return pi_tiny_[x];

    uint64_t count = pi_cache_[x / 240].count;
    uint64_t bits = pi_cache_[x / 240].bits;
    uint64_t bitmask = unset_larger_[x % 240];
    return count + popcnt64(bits & bitmask);
  }

private:
  struct pi_t
  {
    uint64_t count;
    uint64_t bits;
  };

  void init(uint64_t limit, uint64_t cache_limit, int threads);
  void init_bits(uint64_t low, uint64_t high, uint64_t thread_num);
  void init_count(uint64_t low, uint64_t high, uint64_t thread_num);
  static const pod_array<pi_t, 64> pi_cache_;
  pod_vector<pi_t> pi_;
  pod_vector<uint64_t> counts_;
  uint64_t max_x_;
};

} // namespace

#endif

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
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
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
  PiTable(uint64_t limit, int threads);

  uint64_t size() const
  {
    return limit_ + 1;
  }

  /// Get number of primes <= n
  ALWAYS_INLINE int64_t operator[](uint64_t n) const
  {
    assert(n <= limit_);

    if_unlikely(n < pi_tiny_.size())
      return pi_tiny_[n];

    uint64_t count = pi_[n / 240].count;
    uint64_t bits = pi_[n / 240].bits;
    uint64_t bitmask = unset_larger_[n % 240];
    return count + popcnt64(bits & bitmask);
  }

private:
  struct pi_t
  {
    uint64_t count;
    uint64_t bits;
  };

  void init_bits(uint64_t start, uint64_t stop, uint64_t thread_num);
  void init_count(uint64_t start, uint64_t stop, uint64_t thread_num);
  pod_vector<pi_t> pi_;
  pod_vector<uint64_t> counts_;
  uint64_t limit_;
};

} // namespace

#endif

///
/// @file  PiTable.hpp
/// @brief The PiTable class is a compressed lookup table for
///        prime counts. Each bit in the lookup table corresponds
///        to an odd integer and that bit is set to 1 if the
///        integer is a prime. PiTable uses only (n / 8) bytes of
///        memory and returns the number of primes <= n in O(1)
///        operations.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PITABLE_HPP
#define PITABLE_HPP

#include <popcnt.hpp>
#include <aligned_vector.hpp>
#include <noinline.hpp>
#include <pod_vector.hpp>
#include <unlikely.hpp>

#include <stdint.h>
#include <array>
#include <cassert>

namespace primecount {

class PiTable
{
public:
  NOINLINE PiTable(uint64_t limit, int threads);

  int64_t size() const
  {
    return limit_ + 1;
  }

  /// Get number of primes <= n
  int64_t operator[](uint64_t n) const
  {
    assert(n <= limit_);

    // Since we store only odd numbers in our lookup table,
    // we cannot store 2 which is the only even prime.
    // As a workaround we mark 1 as a prime (1st bit) and
    // add a check to return 0 for pi[1].
    if_unlikely(n == 1)
      return 0;

    uint64_t bitmask = unset_bits_[n % 128];
    uint64_t prime_count = pi_[n / 128].prime_count;
    uint64_t bit_count = popcnt64(pi_[n / 128].bits & bitmask);
    return prime_count + bit_count;
  }

private:
  struct pi_t
  {
    uint64_t prime_count;
    uint64_t bits;
  };

  void init_bits(uint64_t start, uint64_t stop, uint64_t thread_num);
  void init_prime_count(uint64_t start, uint64_t stop, uint64_t thread_num);
  static const std::array<uint64_t, 128> unset_bits_;
  pod_vector<pi_t> pi_;
  aligned_vector<uint64_t> counts_;
  uint64_t limit_;
};

} // namespace

#endif

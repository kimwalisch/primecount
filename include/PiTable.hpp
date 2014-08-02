///
/// @file    PiTable.hpp
/// @brief   The PiTable class is a compressed lookup table for prime
///          counts. It uses only n / 64 * 12 bytes of memory and
///          returns the number of primes below n using O(1) operations.
///
///          How it works:
///
///          The prime count is calculated in 2 steps:
///          1) pi(n - n % 64) is stored in pi_[n / 64].prime_count
///          2) The remaining primes inside ]n - n % 64, n]
///             are stored in pi_[n / 64].bits, each bit corresponds
///             to an integer i.e. bit[13] = n - n % 64 + 13. If the
///             bit is set the corresponding integer is a prime.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PITABLE_HPP
#define PITABLE_HPP

#include <popcount.hpp>

#include <stdint.h>
#include <cassert>
#include <vector>

namespace primecount {

class PiTable
{
public:
  PiTable(uint64_t max);

  /// Get the number of primes <= n.
  /// This implementation uses only 20 arithmetic operations.
  ///
  int64_t operator()(uint64_t n) const
  {
    assert(n <= max_);
    uint64_t bitmask = UINT64_C(0xffffffffffffffff) >> (63 - n % 64);
    return pi_[n / 64].prime_count + popcount64(pi_[n / 64].bits & bitmask);
  }

  int64_t size() const
  {
    return max_ + 1;
  }
private:
  struct PiData
  {
    PiData() : prime_count(0), bits(0) { }
    uint64_t prime_count;
    uint64_t bits;
  };

  std::vector<PiData> pi_;
  uint64_t max_;
};

} // namespace

#endif

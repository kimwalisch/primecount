///
/// @file   PiTable.hpp
/// @brief  The PiTable class is a compressed lookup table for prime
///         counts. It uses only n / 64 * 16 bytes of memory and
///         returns the number of primes below n using O(1) operations.
///
///         How it works:
///
///         1) pi_[n / 64].prime_count = pi(n - n % 64).
///
///         2) pi_[n / 64].bits corresponds to ]base, base + 64] with
///            base = n - n % 64, each bit represents one number and
///            if the bit is 1 the number is a prime. Using the
///            POPCNT instruction we can now easily count the primes
///            within ]n - n % 64, n] and add it to 1).
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

  /// @return  The number of primes <= n
  int64_t operator[](uint64_t n) const
  {
    assert(n <= max_);
    uint64_t bitmask = UINT64_C(0xffffffffffffffff) >> (63 - n % 64);
    return pi_[n / 64].prime_count + popcount_u64(pi_[n / 64].bits & bitmask);
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

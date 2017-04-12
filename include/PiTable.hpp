///
/// @file   PiTable.hpp
/// @brief  The PiTable class is a compressed lookup table for prime
///         counts. It uses only (n / 4) bytes of memory and
///         returns the number of primes below n in O(1) operations.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PITABLE_HPP
#define PITABLE_HPP

#include <popcnt.hpp>

#include <stdint.h>
#include <cassert>
#include <vector>

namespace primecount {

class PiTable
{
public:
  PiTable(uint64_t max);

  /// Get number of primes <= n
  int64_t operator[](uint64_t n) const
  {
    assert(n <= max_);
    uint64_t bitmask = 0xffffffffffffffffull >> (63 - n % 64);
    return pi_[n / 64].prime_count + popcnt64(pi_[n / 64].bits & bitmask);
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

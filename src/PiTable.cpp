///
/// @file  PiTable.cpp
/// @brief The PiTable class is a compressed lookup table for
///        prime counts. It uses only (n / 4) bytes of memory
///        and returns the number of primes <= n in O(1)
///        operations.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primesieve.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

PiTable::PiTable(uint64_t max) :
  max_(max)
{
  pi_.resize(max / 64 + 1);
  primesieve::iterator it(0, max);

  uint64_t pix = 0;
  uint64_t prime = 0;

  while ((prime = it.next_prime()) <= max)
    pi_[prime / 64].bits |= 1ull << (prime % 64);

  for (auto& i : pi_)
  {
    i.prime_count = pix;
    pix += popcnt64(i.bits);
  }
}

} // namespace

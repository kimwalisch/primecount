///
/// @file  PiTable.cpp
/// @see   PiTable.hpp for documentation.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <bit_sieve.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

PiTable::PiTable(uint64_t max) :
  max_(max)
{
  init();
}

void PiTable::init()
{
  bit_sieve sieve(max_ + 1);
  sieve.memset(0);

  // sieve of Eratosthenes
  for (uint64_t i = 2; i * i <= max_; i++)
    if (sieve[i])
      for (uint64_t j = i * i; j <= max_; j += i)
        sieve.unset(j);

  uint64_t index = ~0u;
  uint32_t pix = 1;
  pi_.resize(max_ / 64 + 1);

  // fill pi_[x] table
  for (uint64_t x = 2; x <= max_; x++)
  {
    if (x / 64 != index)
      pi_[x / 64].prime_count = pix;
    index = x / 64;

    // check whether x is a prime
    if (sieve[x])
    {
      // set bit corresponding to x to 1
      uint64_t one = 1;
      pi_[x / 64].bits |= one << (x % 64);
      pix++;
    }
  }
}

} // namespace

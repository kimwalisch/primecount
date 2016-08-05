///
/// @file  PiTable.cpp
/// @see   PiTable.hpp for documentation.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
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
    pi_[prime / 64].bits |= ((uint64_t) 1) << (prime % 64);

  for (uint64_t i = 0; i < pi_.size(); i++)
  {
    pi_[i].prime_count = pix;
    pix += popcount_u64(pi_[i].bits);
  }
}

} // namespace

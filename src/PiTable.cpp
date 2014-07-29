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
#include <primesieve.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

PiTable::PiTable(uint64_t max) :
  max_(max)
{
  pi_.resize(max_ / 64 + 1);
  primesieve::iterator iter(0, max_);

  uint64_t pix = 0;
  uint64_t prime = 0;
  uint64_t one = 1;

  while ((prime = iter.next_prime()) <= max_)
    pi_[prime / 64].bits |= one << (prime % 64);

  for (uint64_t i = 0; i < pi_.size(); i++)
  {
    pi_[i].prime_count = pix;
    pix += popcount64(pi_[i].bits);
  }
}

} // namespace

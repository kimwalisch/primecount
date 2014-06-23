///
/// @file  PiTable.cpp
/// @see   PiTable.hpp for documentation.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#if !defined(__STDC_CONSTANT_MACROS)
  #define __STDC_CONSTANT_MACROS
#endif

#include <stdint.h>
#include <vector>

#include <PiTable.hpp>
#include <bit_sieve.hpp>

namespace primecount {

const uint64_t PiTable::bitmasks_[128] =
{
  UINT64_C(0x0000000000000000), UINT64_C(0x0000000000000001),
  UINT64_C(0x0000000000000001), UINT64_C(0x0000000000000003),
  UINT64_C(0x0000000000000003), UINT64_C(0x0000000000000007),
  UINT64_C(0x0000000000000007), UINT64_C(0x000000000000000f),
  UINT64_C(0x000000000000000f), UINT64_C(0x000000000000001f),
  UINT64_C(0x000000000000001f), UINT64_C(0x000000000000003f),
  UINT64_C(0x000000000000003f), UINT64_C(0x000000000000007f),
  UINT64_C(0x000000000000007f), UINT64_C(0x00000000000000ff),
  UINT64_C(0x00000000000000ff), UINT64_C(0x00000000000001ff),
  UINT64_C(0x00000000000001ff), UINT64_C(0x00000000000003ff),
  UINT64_C(0x00000000000003ff), UINT64_C(0x00000000000007ff),
  UINT64_C(0x00000000000007ff), UINT64_C(0x0000000000000fff),
  UINT64_C(0x0000000000000fff), UINT64_C(0x0000000000001fff),
  UINT64_C(0x0000000000001fff), UINT64_C(0x0000000000003fff),
  UINT64_C(0x0000000000003fff), UINT64_C(0x0000000000007fff),
  UINT64_C(0x0000000000007fff), UINT64_C(0x000000000000ffff),
  UINT64_C(0x000000000000ffff), UINT64_C(0x000000000001ffff),
  UINT64_C(0x000000000001ffff), UINT64_C(0x000000000003ffff),
  UINT64_C(0x000000000003ffff), UINT64_C(0x000000000007ffff),
  UINT64_C(0x000000000007ffff), UINT64_C(0x00000000000fffff),
  UINT64_C(0x00000000000fffff), UINT64_C(0x00000000001fffff),
  UINT64_C(0x00000000001fffff), UINT64_C(0x00000000003fffff),
  UINT64_C(0x00000000003fffff), UINT64_C(0x00000000007fffff),
  UINT64_C(0x00000000007fffff), UINT64_C(0x0000000000ffffff),
  UINT64_C(0x0000000000ffffff), UINT64_C(0x0000000001ffffff),
  UINT64_C(0x0000000001ffffff), UINT64_C(0x0000000003ffffff),
  UINT64_C(0x0000000003ffffff), UINT64_C(0x0000000007ffffff),
  UINT64_C(0x0000000007ffffff), UINT64_C(0x000000000fffffff),
  UINT64_C(0x000000000fffffff), UINT64_C(0x000000001fffffff),
  UINT64_C(0x000000001fffffff), UINT64_C(0x000000003fffffff),
  UINT64_C(0x000000003fffffff), UINT64_C(0x000000007fffffff),
  UINT64_C(0x000000007fffffff), UINT64_C(0x00000000ffffffff),
  UINT64_C(0x00000000ffffffff), UINT64_C(0x00000001ffffffff),
  UINT64_C(0x00000001ffffffff), UINT64_C(0x00000003ffffffff),
  UINT64_C(0x00000003ffffffff), UINT64_C(0x00000007ffffffff),
  UINT64_C(0x00000007ffffffff), UINT64_C(0x0000000fffffffff),
  UINT64_C(0x0000000fffffffff), UINT64_C(0x0000001fffffffff),
  UINT64_C(0x0000001fffffffff), UINT64_C(0x0000003fffffffff),
  UINT64_C(0x0000003fffffffff), UINT64_C(0x0000007fffffffff),
  UINT64_C(0x0000007fffffffff), UINT64_C(0x000000ffffffffff),
  UINT64_C(0x000000ffffffffff), UINT64_C(0x000001ffffffffff),
  UINT64_C(0x000001ffffffffff), UINT64_C(0x000003ffffffffff),
  UINT64_C(0x000003ffffffffff), UINT64_C(0x000007ffffffffff),
  UINT64_C(0x000007ffffffffff), UINT64_C(0x00000fffffffffff),
  UINT64_C(0x00000fffffffffff), UINT64_C(0x00001fffffffffff),
  UINT64_C(0x00001fffffffffff), UINT64_C(0x00003fffffffffff),
  UINT64_C(0x00003fffffffffff), UINT64_C(0x00007fffffffffff),
  UINT64_C(0x00007fffffffffff), UINT64_C(0x0000ffffffffffff),
  UINT64_C(0x0000ffffffffffff), UINT64_C(0x0001ffffffffffff),
  UINT64_C(0x0001ffffffffffff), UINT64_C(0x0003ffffffffffff),
  UINT64_C(0x0003ffffffffffff), UINT64_C(0x0007ffffffffffff),
  UINT64_C(0x0007ffffffffffff), UINT64_C(0x000fffffffffffff),
  UINT64_C(0x000fffffffffffff), UINT64_C(0x001fffffffffffff),
  UINT64_C(0x001fffffffffffff), UINT64_C(0x003fffffffffffff),
  UINT64_C(0x003fffffffffffff), UINT64_C(0x007fffffffffffff),
  UINT64_C(0x007fffffffffffff), UINT64_C(0x00ffffffffffffff),
  UINT64_C(0x00ffffffffffffff), UINT64_C(0x01ffffffffffffff),
  UINT64_C(0x01ffffffffffffff), UINT64_C(0x03ffffffffffffff),
  UINT64_C(0x03ffffffffffffff), UINT64_C(0x07ffffffffffffff),
  UINT64_C(0x07ffffffffffffff), UINT64_C(0x0fffffffffffffff),
  UINT64_C(0x0fffffffffffffff), UINT64_C(0x1fffffffffffffff),
  UINT64_C(0x1fffffffffffffff), UINT64_C(0x3fffffffffffffff),
  UINT64_C(0x3fffffffffffffff), UINT64_C(0x7fffffffffffffff),
  UINT64_C(0x7fffffffffffffff), UINT64_C(0xffffffffffffffff)
};

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
  for (uint64_t i = 3; i * i <= max_; i += 2)
    if (sieve[i])
      for (uint64_t j = i * i; j <= max_; j += i * 2)
        sieve.unset(j);

  uint64_t index = ~0u;
  uint32_t pix = 1;
  pi_.resize(max_ / 128 + 1);

  // fill pi_[x] table
  for (uint64_t x = 3; x <= max_; x += 2)
  {
    if (x / 128 != index)
      pi_[x / 128].prime_count = pix;
    index = x / 128;

    // check whether x is a prime
    if (sieve[x])
    {
      // set bit corresponding to x to 1
      uint64_t one = 1;
      pi_[x / 128].bits |= one << ((x % 128) / 2);
      pix++;
    }
  }

  if (max_ > 0 && max_ % 128 == 0)
    pi_[max_ / 128].prime_count = pix;
}

} // namespace

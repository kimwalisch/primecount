///
/// @file  conditional_move.cpp
/// @brief For primecount's performance, it is important that the
///        compiler emits branchfree conditional selects for our
///        ARM64 CONDITIONAL_MOVE macro usage.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

extern const uint64_t unset_smaller[240];
extern const uint64_t unset_larger[240];

uint64_t sieve_count(uint64_t start, uint64_t stop, const uint64_t* sieve64)
{
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  /* Branchfree bitmask calculation: */
  /* if (start_idx == stop_idx) m1 = m1 & m2; */
  /* if (start_idx == stop_idx) m2 = 0; */
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2);
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0);

  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  uint64_t cnt = popcnt64(start_bits);
  cnt += popcnt64(stop_bits);

  for (uint64_t i = start_idx + 1; i < stop_idx; i++)
    cnt += popcnt64(sieve64[i]);

  return cnt;
}

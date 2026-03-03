///
/// @file  Sieve.cpp
/// @brief This file implements a highly optimized prime sieving
///        algorithm for computing the special leaves (sometimes named
///        hard special leaves) in the combinatorial prime counting
///        algorithms (e.g. Lagarias-Miller-Odlyzko, Deleglise-Rivat,
///        Gourdon).
///
///        The Sieve class contains a sieve of Eratosthenes
///        implementation with 30 numbers per byte i.e. the 8 bits of
///        each byte correspond to the offsets: { 1, 7, 11, 13, 17,
///        19, 23, 29 }. Unlike a traditional prime sieve this sieve
///        is designed for use in the combinatorial prime counting
///        algorithms: this sieve removes primes as well as multiples
///        of primes and it counts the number of elements that have
///        been crossed off for the first time in the sieve array.
///
///        Since there is a large number of leaves for which we have
///        to count the number of unsieved elements in the sieve
///        array, Lagarias-Miller-Odlyzko have suggested using a
///        binary indexed tree data structure (a.k.a. Fenwick tree) to
///        speedup counting. However using a binary indexed tree is
///        bad for performance as it causes many cache misses and
///        branch mispredictions. For this reason this implementation
///        instead uses a linear counter array whose elements contain
///        the total count of unsieved elements in a certain interval.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "Sieve.hpp"
#include "Sieve_arrays.hpp"
#include "Sieve_count_start_stop.hpp"
#include "Sieve_pre_sieve.hpp"

#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {

namespace {

constexpr uint8_t wheel_factor[8] = { 6, 4, 2, 4, 2, 4, 6, 2 };

constexpr uint8_t wheel_bits[8][8] = {
  { 0, 1, 2, 3, 4, 5, 6, 7 }, // p % 30 == 1
  { 1, 5, 4, 0, 7, 3, 2, 6 }, // p % 30 == 7
  { 2, 4, 0, 6, 1, 7, 3, 5 }, // p % 30 == 11
  { 3, 0, 6, 5, 2, 1, 7, 4 }, // p % 30 == 13
  { 4, 7, 1, 2, 5, 6, 0, 3 }, // p % 30 == 17
  { 5, 3, 7, 1, 6, 0, 4, 2 }, // p % 30 == 19
  { 6, 2, 3, 7, 0, 4, 5, 1 }, // p % 30 == 23
  { 7, 6, 5, 4, 3, 2, 1, 0 }  // p % 30 == 29
};

constexpr uint8_t wheel_corr[8][8] = {
  { 0, 0, 0, 0, 0, 0, 0, 1 }, // p % 30 == 1
  { 1, 1, 1, 0, 1, 1, 1, 1 }, // p % 30 == 7
  { 2, 2, 0, 2, 0, 2, 2, 1 }, // p % 30 == 11
  { 3, 1, 1, 2, 1, 1, 3, 1 }, // p % 30 == 13
  { 3, 3, 1, 2, 1, 3, 3, 1 }, // p % 30 == 17
  { 4, 2, 2, 2, 2, 2, 4, 1 }, // p % 30 == 19
  { 5, 3, 1, 4, 1, 3, 5, 1 }, // p % 30 == 23
  { 6, 4, 2, 4, 2, 4, 6, 1 }  // p % 30 == 29
};

} // namespace

Sieve::Sieve(uint64_t low,
             uint64_t segment_size,
             uint64_t primes_size)
{
  ASSERT(low % 30 == 0);
  ASSERT(segment_size % 240 == 0);

  start_ = low;
  segment_size = align_segment_size(segment_size);

  // sieve_size = segment_size / 30 as each byte corresponds
  // to 30 numbers i.e. the 8 bits correspond to the
  // offsets = {1, 7, 11, 13, 17, 19, 23, 29}.
  sieve_.resize(segment_size / 30);
  primeState_.reserve(primes_size);
}

/// The segment size is sieve.size() * 30 as each
/// byte corresponds to 30 numbers.
///
uint64_t Sieve::segment_size() const
{
  return sieve_.size() * 30;
}

/// segment_size must be a multiple of 240 as we
/// process 64-bit words (8 bytes) and each
/// byte contains 30 numbers.
///
uint64_t Sieve::align_segment_size(uint64_t size)
{
  size = max(size, 240);

  if (size % 240)
    size += 240 - size % 240;

  return size;
}

/// Use shorter sieve size for last segment
void Sieve::resize_sieve(uint64_t low, uint64_t high)
{
  uint64_t size = high - low;

  if (size < segment_size())
  {
    uint64_t last = size - 1;
    size = align_segment_size(size);
    sieve_.resize(size / 30);
    auto sieve64 = (uint64_t*) sieve_.data();
    sieve64[last / 240] &= unset_larger[last % 240];
  }
}

void Sieve::reset_counter()
{
  prev_stop_ = 0;
  count_ = 0;
}

void Sieve::init_counter(uint64_t low, uint64_t high)
{
  count_ = 0;
  total_count_ = 0;

  if (high > low)
    total_count_ = count(0, (high - 1) - low);
}

/// Add a sieving prime to the sieve.
/// Calculates the first multiple > start of prime that
/// is not divisible by 2, 3, 5 and its wheel index.
///
void Sieve::add(uint64_t prime, uint64_t i)
{
  if_unlikely(i > primeState_.size())
    primeState_.resize(i);

  // Find first multiple > start_
  ASSERT(start_ % 30 == 0);
  uint64_t quotient = start_ / prime + 1;
  uint64_t multiple = prime * quotient;

  // Find next multiple of prime that
  // is not divisible by 2, 3, 5
  uint64_t factor = wheel_init[quotient % 30].factor;
  multiple += prime * factor;

  ASSERT(multiple % 2 != 0);
  ASSERT(multiple % 3 != 0);
  ASSERT(multiple % 5 != 0);

  multiple = (multiple - start_) / 30;
  uint32_t multiple32 = (uint32_t) multiple;

  // Calculate wheel index of multiple
  uint32_t index = wheel_init[quotient % 30].index;
  index += wheel_init_offsets[prime % 30];
  primeState_.push_back({multiple32, index});
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Used for pre-sieving.
///
void Sieve::cross_off(uint64_t prime, uint64_t i)
{
  if (i >= primeState_.size())
    add(prime, i);

  prime /= 30;
  PrimeState& primeState = primeState_[i];
  uint64_t m = primeState.multiple;
  uint8_t* sieve = sieve_.data();
  uint64_t sieve_size = sieve_.size();

  #define CHECK_FINISHED(i) \
    if_unlikely(m >= sieve_size) \
    { \
      primeState.wheel_index = i; \
      primeState.multiple = (uint32_t) (m - sieve_size); \
      return; \
    }

  ASSERT(primeState.wheel_index <= 63);
  switch (primeState.wheel_index)
  {
    for (;;)
    {
      case 0: {
                uint64_t max_offset = m + prime * 28;
                uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                for (; m < limit; m += prime * 30 + 1)
                {
                  sieve[m + prime *  0] &= ~(1 << 0);
                  sieve[m + prime *  6] &= ~(1 << 1);
                  sieve[m + prime * 10] &= ~(1 << 2);
                  sieve[m + prime * 12] &= ~(1 << 3);
                  sieve[m + prime * 16] &= ~(1 << 4);
                  sieve[m + prime * 18] &= ~(1 << 5);
                  sieve[m + prime * 22] &= ~(1 << 6);
                  sieve[m + prime * 28] &= ~(1 << 7);
                }
              }
              CHECK_FINISHED(0); sieve[m] &= ~(1 << 0); m += prime * 6 + 0; FALLTHROUGH;
      case 1: CHECK_FINISHED(1); sieve[m] &= ~(1 << 1); m += prime * 4 + 0; FALLTHROUGH;
      case 2: CHECK_FINISHED(2); sieve[m] &= ~(1 << 2); m += prime * 2 + 0; FALLTHROUGH;
      case 3: CHECK_FINISHED(3); sieve[m] &= ~(1 << 3); m += prime * 4 + 0; FALLTHROUGH;
      case 4: CHECK_FINISHED(4); sieve[m] &= ~(1 << 4); m += prime * 2 + 0; FALLTHROUGH;
      case 5: CHECK_FINISHED(5); sieve[m] &= ~(1 << 5); m += prime * 4 + 0; FALLTHROUGH;
      case 6: CHECK_FINISHED(6); sieve[m] &= ~(1 << 6); m += prime * 6 + 0; FALLTHROUGH;
      case 7: CHECK_FINISHED(7); sieve[m] &= ~(1 << 7); m += prime * 2 + 1;
    }

    for (;;)
    {
      case  8: {
                 uint64_t max_offset = m + prime * 28 + 6;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 7)
                 {
                   sieve[m + prime *  0 + 0] &= ~(1 << 1);
                   sieve[m + prime *  6 + 1] &= ~(1 << 5);
                   sieve[m + prime * 10 + 2] &= ~(1 << 4);
                   sieve[m + prime * 12 + 3] &= ~(1 << 0);
                   sieve[m + prime * 16 + 3] &= ~(1 << 7);
                   sieve[m + prime * 18 + 4] &= ~(1 << 3);
                   sieve[m + prime * 22 + 5] &= ~(1 << 2);
                   sieve[m + prime * 28 + 6] &= ~(1 << 6);
                 }
               }
               CHECK_FINISHED( 8); sieve[m] &= ~(1 << 1); m += prime * 6 + 1; FALLTHROUGH;
      case  9: CHECK_FINISHED( 9); sieve[m] &= ~(1 << 5); m += prime * 4 + 1; FALLTHROUGH;
      case 10: CHECK_FINISHED(10); sieve[m] &= ~(1 << 4); m += prime * 2 + 1; FALLTHROUGH;
      case 11: CHECK_FINISHED(11); sieve[m] &= ~(1 << 0); m += prime * 4 + 0; FALLTHROUGH;
      case 12: CHECK_FINISHED(12); sieve[m] &= ~(1 << 7); m += prime * 2 + 1; FALLTHROUGH;
      case 13: CHECK_FINISHED(13); sieve[m] &= ~(1 << 3); m += prime * 4 + 1; FALLTHROUGH;
      case 14: CHECK_FINISHED(14); sieve[m] &= ~(1 << 2); m += prime * 6 + 1; FALLTHROUGH;
      case 15: CHECK_FINISHED(15); sieve[m] &= ~(1 << 6); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 16: {
                 uint64_t max_offset = m + prime * 28 + 10;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 11)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 2);
                   sieve[m + prime *  6 +  2] &= ~(1 << 4);
                   sieve[m + prime * 10 +  4] &= ~(1 << 0);
                   sieve[m + prime * 12 +  4] &= ~(1 << 6);
                   sieve[m + prime * 16 +  6] &= ~(1 << 1);
                   sieve[m + prime * 18 +  6] &= ~(1 << 7);
                   sieve[m + prime * 22 +  8] &= ~(1 << 3);
                   sieve[m + prime * 28 + 10] &= ~(1 << 5);
                 }
               }
               CHECK_FINISHED(16); sieve[m] &= ~(1 << 2); m += prime * 6 + 2; FALLTHROUGH;
      case 17: CHECK_FINISHED(17); sieve[m] &= ~(1 << 4); m += prime * 4 + 2; FALLTHROUGH;
      case 18: CHECK_FINISHED(18); sieve[m] &= ~(1 << 0); m += prime * 2 + 0; FALLTHROUGH;
      case 19: CHECK_FINISHED(19); sieve[m] &= ~(1 << 6); m += prime * 4 + 2; FALLTHROUGH;
      case 20: CHECK_FINISHED(20); sieve[m] &= ~(1 << 1); m += prime * 2 + 0; FALLTHROUGH;
      case 21: CHECK_FINISHED(21); sieve[m] &= ~(1 << 7); m += prime * 4 + 2; FALLTHROUGH;
      case 22: CHECK_FINISHED(22); sieve[m] &= ~(1 << 3); m += prime * 6 + 2; FALLTHROUGH;
      case 23: CHECK_FINISHED(23); sieve[m] &= ~(1 << 5); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 24: {
                 uint64_t max_offset = m + prime * 28 + 12;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 13)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 3);
                   sieve[m + prime *  6 +  3] &= ~(1 << 0);
                   sieve[m + prime * 10 +  4] &= ~(1 << 6);
                   sieve[m + prime * 12 +  5] &= ~(1 << 5);
                   sieve[m + prime * 16 +  7] &= ~(1 << 2);
                   sieve[m + prime * 18 +  8] &= ~(1 << 1);
                   sieve[m + prime * 22 +  9] &= ~(1 << 7);
                   sieve[m + prime * 28 + 12] &= ~(1 << 4);
                 }
               }
               CHECK_FINISHED(24); sieve[m] &= ~(1 << 3); m += prime * 6 + 3; FALLTHROUGH;
      case 25: CHECK_FINISHED(25); sieve[m] &= ~(1 << 0); m += prime * 4 + 1; FALLTHROUGH;
      case 26: CHECK_FINISHED(26); sieve[m] &= ~(1 << 6); m += prime * 2 + 1; FALLTHROUGH;
      case 27: CHECK_FINISHED(27); sieve[m] &= ~(1 << 5); m += prime * 4 + 2; FALLTHROUGH;
      case 28: CHECK_FINISHED(28); sieve[m] &= ~(1 << 2); m += prime * 2 + 1; FALLTHROUGH;
      case 29: CHECK_FINISHED(29); sieve[m] &= ~(1 << 1); m += prime * 4 + 1; FALLTHROUGH;
      case 30: CHECK_FINISHED(30); sieve[m] &= ~(1 << 7); m += prime * 6 + 3; FALLTHROUGH;
      case 31: CHECK_FINISHED(31); sieve[m] &= ~(1 << 4); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 32: {
                 uint64_t max_offset = m + prime * 28 + 16;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 17)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 4);
                   sieve[m + prime *  6 +  3] &= ~(1 << 7);
                   sieve[m + prime * 10 +  6] &= ~(1 << 1);
                   sieve[m + prime * 12 +  7] &= ~(1 << 2);
                   sieve[m + prime * 16 +  9] &= ~(1 << 5);
                   sieve[m + prime * 18 + 10] &= ~(1 << 6);
                   sieve[m + prime * 22 + 13] &= ~(1 << 0);
                   sieve[m + prime * 28 + 16] &= ~(1 << 3);
                 }
               }
               CHECK_FINISHED(32); sieve[m] &= ~(1 << 4); m += prime * 6 + 3; FALLTHROUGH;
      case 33: CHECK_FINISHED(33); sieve[m] &= ~(1 << 7); m += prime * 4 + 3; FALLTHROUGH;
      case 34: CHECK_FINISHED(34); sieve[m] &= ~(1 << 1); m += prime * 2 + 1; FALLTHROUGH;
      case 35: CHECK_FINISHED(35); sieve[m] &= ~(1 << 2); m += prime * 4 + 2; FALLTHROUGH;
      case 36: CHECK_FINISHED(36); sieve[m] &= ~(1 << 5); m += prime * 2 + 1; FALLTHROUGH;
      case 37: CHECK_FINISHED(37); sieve[m] &= ~(1 << 6); m += prime * 4 + 3; FALLTHROUGH;
      case 38: CHECK_FINISHED(38); sieve[m] &= ~(1 << 0); m += prime * 6 + 3; FALLTHROUGH;
      case 39: CHECK_FINISHED(39); sieve[m] &= ~(1 << 3); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 40: {
                 uint64_t max_offset = m + prime * 28 + 18;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 19)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 5);
                   sieve[m + prime *  6 +  4] &= ~(1 << 3);
                   sieve[m + prime * 10 +  6] &= ~(1 << 7);
                   sieve[m + prime * 12 +  8] &= ~(1 << 1);
                   sieve[m + prime * 16 + 10] &= ~(1 << 6);
                   sieve[m + prime * 18 + 12] &= ~(1 << 0);
                   sieve[m + prime * 22 + 14] &= ~(1 << 4);
                   sieve[m + prime * 28 + 18] &= ~(1 << 2);
                 }
               }
               CHECK_FINISHED(40); sieve[m] &= ~(1 << 5); m += prime * 6 + 4; FALLTHROUGH;
      case 41: CHECK_FINISHED(41); sieve[m] &= ~(1 << 3); m += prime * 4 + 2; FALLTHROUGH;
      case 42: CHECK_FINISHED(42); sieve[m] &= ~(1 << 7); m += prime * 2 + 2; FALLTHROUGH;
      case 43: CHECK_FINISHED(43); sieve[m] &= ~(1 << 1); m += prime * 4 + 2; FALLTHROUGH;
      case 44: CHECK_FINISHED(44); sieve[m] &= ~(1 << 6); m += prime * 2 + 2; FALLTHROUGH;
      case 45: CHECK_FINISHED(45); sieve[m] &= ~(1 << 0); m += prime * 4 + 2; FALLTHROUGH;
      case 46: CHECK_FINISHED(46); sieve[m] &= ~(1 << 4); m += prime * 6 + 4; FALLTHROUGH;
      case 47: CHECK_FINISHED(47); sieve[m] &= ~(1 << 2); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 48: {
                 uint64_t max_offset = m + prime * 28 + 22;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 23)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 6);
                   sieve[m + prime *  6 +  5] &= ~(1 << 2);
                   sieve[m + prime * 10 +  8] &= ~(1 << 3);
                   sieve[m + prime * 12 +  9] &= ~(1 << 7);
                   sieve[m + prime * 16 + 13] &= ~(1 << 0);
                   sieve[m + prime * 18 + 14] &= ~(1 << 4);
                   sieve[m + prime * 22 + 17] &= ~(1 << 5);
                   sieve[m + prime * 28 + 22] &= ~(1 << 1);
                 }
               }
               CHECK_FINISHED(48); sieve[m] &= ~(1 << 6); m += prime * 6 + 5; FALLTHROUGH;
      case 49: CHECK_FINISHED(49); sieve[m] &= ~(1 << 2); m += prime * 4 + 3; FALLTHROUGH;
      case 50: CHECK_FINISHED(50); sieve[m] &= ~(1 << 3); m += prime * 2 + 1; FALLTHROUGH;
      case 51: CHECK_FINISHED(51); sieve[m] &= ~(1 << 7); m += prime * 4 + 4; FALLTHROUGH;
      case 52: CHECK_FINISHED(52); sieve[m] &= ~(1 << 0); m += prime * 2 + 1; FALLTHROUGH;
      case 53: CHECK_FINISHED(53); sieve[m] &= ~(1 << 4); m += prime * 4 + 3; FALLTHROUGH;
      case 54: CHECK_FINISHED(54); sieve[m] &= ~(1 << 5); m += prime * 6 + 5; FALLTHROUGH;
      case 55: CHECK_FINISHED(55); sieve[m] &= ~(1 << 1); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 56: {
                 uint64_t max_offset = m + prime * 28 + 28;
                 uint64_t limit = std::max(max_offset, sieve_size) - max_offset;

                 for (; m < limit; m += prime * 30 + 29)
                 {
                   sieve[m + prime *  0 +  0] &= ~(1 << 7);
                   sieve[m + prime *  6 +  6] &= ~(1 << 6);
                   sieve[m + prime * 10 + 10] &= ~(1 << 5);
                   sieve[m + prime * 12 + 12] &= ~(1 << 4);
                   sieve[m + prime * 16 + 16] &= ~(1 << 3);
                   sieve[m + prime * 18 + 18] &= ~(1 << 2);
                   sieve[m + prime * 22 + 22] &= ~(1 << 1);
                   sieve[m + prime * 28 + 28] &= ~(1 << 0);
                 }
               }
               CHECK_FINISHED(56); sieve[m] &= ~(1 << 7); m += prime * 6 + 6; FALLTHROUGH;
      case 57: CHECK_FINISHED(57); sieve[m] &= ~(1 << 6); m += prime * 4 + 4; FALLTHROUGH;
      case 58: CHECK_FINISHED(58); sieve[m] &= ~(1 << 5); m += prime * 2 + 2; FALLTHROUGH;
      case 59: CHECK_FINISHED(59); sieve[m] &= ~(1 << 4); m += prime * 4 + 4; FALLTHROUGH;
      case 60: CHECK_FINISHED(60); sieve[m] &= ~(1 << 3); m += prime * 2 + 2; FALLTHROUGH;
      case 61: CHECK_FINISHED(61); sieve[m] &= ~(1 << 2); m += prime * 4 + 4; FALLTHROUGH;
      case 62: CHECK_FINISHED(62); sieve[m] &= ~(1 << 1); m += prime * 6 + 6; FALLTHROUGH;
      case 63: CHECK_FINISHED(63); sieve[m] &= ~(1 << 0); m += prime * 2 + 1;
    }

    default: UNREACHABLE;
  }

  #undef CHECK_FINISHED
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Also counts the number of elements
/// removed for the first time i.e. the count of sieved elements
/// whose least prime factor is the i-th prime.
///
void Sieve::cross_off_count(uint64_t prime, uint64_t i)
{
  if (i >= primeState_.size())
    add(prime, i);

  reset_counter();
  PrimeState& primeState = primeState_[i];
  ASSERT(primeState.wheel_index <= 63);
  uint32_t r = primeState.wheel_index >> 3;
  prime /= 30;

  const Array<uint64_t, 8> adv = {
    prime * wheel_factor[0] + wheel_corr[r][0],
    prime * wheel_factor[1] + wheel_corr[r][1],
    prime * wheel_factor[2] + wheel_corr[r][2],
    prime * wheel_factor[3] + wheel_corr[r][3],
    prime * wheel_factor[4] + wheel_corr[r][4],
    prime * wheel_factor[5] + wheel_corr[r][5],
    prime * wheel_factor[6] + wheel_corr[r][6],
    prime * wheel_factor[7] + wheel_corr[r][7]
  };

  const Array<uint8_t, 8> bitmask = {
    uint8_t(1u << wheel_bits[r][0]),
    uint8_t(1u << wheel_bits[r][1]),
    uint8_t(1u << wheel_bits[r][2]),
    uint8_t(1u << wheel_bits[r][3]),
    uint8_t(1u << wheel_bits[r][4]),
    uint8_t(1u << wheel_bits[r][5]),
    uint8_t(1u << wheel_bits[r][6]),
    uint8_t(1u << wheel_bits[r][7])
  };

  #define CHECK_FINISHED(i) \
    if_unlikely(m >= sieve_size) \
    { \
      primeState.wheel_index = (r << 3) + i; \
      primeState.multiple = (uint32_t) (m - sieve_size); \
      total_count_ = total_count; \
      return; \
    }

  #define COUNT_UNSET_BIT(i) \
    { \
      std::size_t b = sieve[m]; \
      std::size_t is_bit = (b & bitmask[i]) != 0; \
      sieve[m] = uint8_t(b & ~bitmask[i]); \
      total_count -= (uint64_t) is_bit; \
      m += adv[i]; \
    }

  uint64_t m = primeState.multiple;
  uint32_t s = primeState.wheel_index & 7;
  uint64_t sieve_size = sieve_.size();
  uint8_t* sieve = &sieve_[0];
  uint64_t total_count = total_count_;

  // Get ready for loop unrolling
  for (; s; s = (s + 1) & 7)
  {
    CHECK_FINISHED(s);
    COUNT_UNSET_BIT(s);
  }

  while (true)
  {
    CHECK_FINISHED(0); COUNT_UNSET_BIT(0);
    CHECK_FINISHED(1); COUNT_UNSET_BIT(1);
    CHECK_FINISHED(2); COUNT_UNSET_BIT(2);
    CHECK_FINISHED(3); COUNT_UNSET_BIT(3);
    CHECK_FINISHED(4); COUNT_UNSET_BIT(4);
    CHECK_FINISHED(5); COUNT_UNSET_BIT(5);
    CHECK_FINISHED(6); COUNT_UNSET_BIT(6);
    CHECK_FINISHED(7); COUNT_UNSET_BIT(7);
  }
}

} // namespace

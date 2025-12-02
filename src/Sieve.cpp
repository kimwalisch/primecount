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
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
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
#include <cmath>

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <cpu_supports_avx512_vpopcnt.hpp>
#endif

namespace primecount {

Sieve::Sieve(uint64_t low,
             uint64_t segment_size,
             uint64_t wheel_size)
{
  ASSERT(low % 30 == 0);
  ASSERT(segment_size % 240 == 0);

  start_ = low;
  segment_size = align_segment_size(segment_size);

  // sieve_size = segment_size / 30 as each byte corresponds
  // to 30 numbers i.e. the 8 bits correspond to the
  // offsets = {1, 7, 11, 13, 17, 19, 23, 29}.
  sieve_.resize(segment_size / 30);
  wheel_.reserve(wheel_size);
  allocate_counter(low);
}

/// Each element of the counter array contains the current
/// number of unsieved elements in the interval:
/// [i * counter_.dist, (i + 1) * counter_.dist[.
/// Ideally each element of the counter array should
/// represent an interval of size O(sqrt(average_leaf_dist)).
/// Also the counter distance should be adjusted regularly
/// whilst sieving as the distance between consecutive
/// leaves is very small ~ log(x) at the beginning of the
/// sieving algorithm but grows up to segment_size towards
/// the end of the algorithm.
///
void Sieve::allocate_counter(uint64_t low)
{
  double average_leaf_dist = std::sqrt(low);
  double counter_dist = std::sqrt(average_leaf_dist);

  // Here we balance counting with the counter array and
  // counting from the sieve array using the 64-bit POPCNT
  // instruction. Since the 64-bit POPCNT instructions
  // allows to count a distance of 240 using a single
  // instruction we slightly increase the counter distance
  // and slightly decrease the size of the counter array.
  uint64_t bytes_count_instruction = bytes_per_count_instruction();
  ASSERT(bytes_count_instruction >= sizeof(uint64_t));
  uint64_t dist_per_instruction = bytes_count_instruction * 30;
  counter_.dist = uint64_t(counter_dist * std::sqrt(dist_per_instruction));

  // Increasing the minimum counter distance decreases the
  // branch mispredictions (good) but on the other hand
  // increases the number of executed instructions (bad).
  // In my benchmarks setting the minimum amount of bytes to
  // bytes_count_instruction * 8 (or 16) performed best.
  uint64_t bytes = counter_.dist / 30;
  bytes = max(bytes, bytes_count_instruction * 8);
  bytes = next_power_of_2(bytes);

  // Make sure the counter (32-bit) doesn't overflow.
  // This can never happen since each counter array element
  // only counts the number of unsieved elements (1 bits) in
  // an interval of size: sieve_limit^(1/4) * sqrt(240).
  // Hence the max(counter value) = 2^18.
  ASSERT(bytes * 8 <= pstd::numeric_limits<uint32_t>::max());
  uint64_t counter_size = ceil_div(sieve_.size(), bytes);
  counter_.counter.resize(counter_size);
  counter_.dist = bytes * 30;
  counter_.log2_dist = ilog2(bytes);
}

string_view_t Sieve::count_algo_name()
{
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    if (cpu_supports_avx512_vpopcnt)
      return "Algorithm: AVX512 bit counting";
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_avx512_vpopcnt)
      return "Algorithm: ARM SVE bit counting";
  #endif

  return "Algorithm: " DEFAULT_SIEVE_COUNT_ALGO_NAME;
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
  counter_.i = 0;
  counter_.sum = 0;
  counter_.stop = counter_.dist;
}

void Sieve::init_counter(uint64_t low, uint64_t high)
{
  reset_counter();
  total_count_ = 0;

  uint64_t start = 0;
  uint64_t max_stop = (high - 1) - low;

  while (start <= max_stop)
  {
    uint64_t stop = start + counter_.dist - 1;
    stop = min(stop, max_stop);
    uint64_t cnt = count(start, stop);
    uint64_t byte_index = start / 30;
    uint64_t i = byte_index >> counter_.log2_dist;

    counter_[i] = (uint32_t) cnt;
    total_count_ += cnt;
    start += counter_.dist;
  }
}

/// Add a sieving prime to the sieve.
/// Calculates the first multiple > start of prime that
/// is not divisible by 2, 3, 5 and its wheel index.
///
void Sieve::add(uint64_t prime, uint64_t i)
{
  if_unlikely(i > wheel_.size())
    wheel_.resize(i);

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
  index += wheel_offsets[prime % 30];
  wheel_.push_back({multiple32, index});
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Used for pre-sieving.
///
void Sieve::cross_off(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime, i);

  prime /= 30;
  Wheel& wheel = wheel_[i];
  uint64_t m = wheel.multiple;
  uint8_t* sieve = sieve_.data();
  uint64_t sieve_size = sieve_.size();

  #define CHECK_FINISHED(wheel_index) \
    if_unlikely(m >= sieve_size) \
    { \
      wheel.index = wheel_index; \
      wheel.multiple = (uint32_t) (m - sieve_size); \
      return; \
    }

  ASSERT(wheel.index <= 63);

  switch (wheel.index)
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
  if (i >= wheel_.size())
    add(prime, i);

  reset_counter();
  Wheel& wheel = wheel_[i];
  prime /= 30;

  uint64_t m = wheel.multiple;
  uint64_t total_count = total_count_;
  uint64_t counter_log2_dist = counter_.log2_dist;
  uint64_t sieve_size = sieve_.size();
  uint32_t* counter = &counter_[0];
  uint8_t* sieve = &sieve_[0];

  #define CHECK_FINISHED(wheel_index) \
    if_unlikely(m >= sieve_size) \
    { \
      wheel.index = wheel_index; \
      wheel.multiple = (uint32_t) (m - sieve_size); \
      total_count_ = total_count; \
      return; \
    }

  #define COUNT_UNSET_BIT(bit_index) \
    { \
      std::size_t sieve_byte = sieve[m]; \
      std::size_t is_bit = (sieve_byte >> bit_index) & 1; \
      sieve[m] &= ~(1 << bit_index); \
      counter[m >> counter_log2_dist] -= (uint32_t) is_bit; \
      total_count -= (uint64_t) is_bit; \
    }

  ASSERT(wheel.index <= 63);

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: CHECK_FINISHED(0); COUNT_UNSET_BIT(0); m += prime * 6 + 0; FALLTHROUGH;
      case 1: CHECK_FINISHED(1); COUNT_UNSET_BIT(1); m += prime * 4 + 0; FALLTHROUGH;
      case 2: CHECK_FINISHED(2); COUNT_UNSET_BIT(2); m += prime * 2 + 0; FALLTHROUGH;
      case 3: CHECK_FINISHED(3); COUNT_UNSET_BIT(3); m += prime * 4 + 0; FALLTHROUGH;
      case 4: CHECK_FINISHED(4); COUNT_UNSET_BIT(4); m += prime * 2 + 0; FALLTHROUGH;
      case 5: CHECK_FINISHED(5); COUNT_UNSET_BIT(5); m += prime * 4 + 0; FALLTHROUGH;
      case 6: CHECK_FINISHED(6); COUNT_UNSET_BIT(6); m += prime * 6 + 0; FALLTHROUGH;
      case 7: CHECK_FINISHED(7); COUNT_UNSET_BIT(7); m += prime * 2 + 1;
    }

    for (;;)
    {
      case  8: CHECK_FINISHED( 8); COUNT_UNSET_BIT(1); m += prime * 6 + 1; FALLTHROUGH;
      case  9: CHECK_FINISHED( 9); COUNT_UNSET_BIT(5); m += prime * 4 + 1; FALLTHROUGH;
      case 10: CHECK_FINISHED(10); COUNT_UNSET_BIT(4); m += prime * 2 + 1; FALLTHROUGH;
      case 11: CHECK_FINISHED(11); COUNT_UNSET_BIT(0); m += prime * 4 + 0; FALLTHROUGH;
      case 12: CHECK_FINISHED(12); COUNT_UNSET_BIT(7); m += prime * 2 + 1; FALLTHROUGH;
      case 13: CHECK_FINISHED(13); COUNT_UNSET_BIT(3); m += prime * 4 + 1; FALLTHROUGH;
      case 14: CHECK_FINISHED(14); COUNT_UNSET_BIT(2); m += prime * 6 + 1; FALLTHROUGH;
      case 15: CHECK_FINISHED(15); COUNT_UNSET_BIT(6); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 16: CHECK_FINISHED(16); COUNT_UNSET_BIT(2); m += prime * 6 + 2; FALLTHROUGH;
      case 17: CHECK_FINISHED(17); COUNT_UNSET_BIT(4); m += prime * 4 + 2; FALLTHROUGH;
      case 18: CHECK_FINISHED(18); COUNT_UNSET_BIT(0); m += prime * 2 + 0; FALLTHROUGH;
      case 19: CHECK_FINISHED(19); COUNT_UNSET_BIT(6); m += prime * 4 + 2; FALLTHROUGH;
      case 20: CHECK_FINISHED(20); COUNT_UNSET_BIT(1); m += prime * 2 + 0; FALLTHROUGH;
      case 21: CHECK_FINISHED(21); COUNT_UNSET_BIT(7); m += prime * 4 + 2; FALLTHROUGH;
      case 22: CHECK_FINISHED(22); COUNT_UNSET_BIT(3); m += prime * 6 + 2; FALLTHROUGH;
      case 23: CHECK_FINISHED(23); COUNT_UNSET_BIT(5); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 24: CHECK_FINISHED(24); COUNT_UNSET_BIT(3); m += prime * 6 + 3; FALLTHROUGH;
      case 25: CHECK_FINISHED(25); COUNT_UNSET_BIT(0); m += prime * 4 + 1; FALLTHROUGH;
      case 26: CHECK_FINISHED(26); COUNT_UNSET_BIT(6); m += prime * 2 + 1; FALLTHROUGH;
      case 27: CHECK_FINISHED(27); COUNT_UNSET_BIT(5); m += prime * 4 + 2; FALLTHROUGH;
      case 28: CHECK_FINISHED(28); COUNT_UNSET_BIT(2); m += prime * 2 + 1; FALLTHROUGH;
      case 29: CHECK_FINISHED(29); COUNT_UNSET_BIT(1); m += prime * 4 + 1; FALLTHROUGH;
      case 30: CHECK_FINISHED(30); COUNT_UNSET_BIT(7); m += prime * 6 + 3; FALLTHROUGH;
      case 31: CHECK_FINISHED(31); COUNT_UNSET_BIT(4); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 32: CHECK_FINISHED(32); COUNT_UNSET_BIT(4); m += prime * 6 + 3; FALLTHROUGH;
      case 33: CHECK_FINISHED(33); COUNT_UNSET_BIT(7); m += prime * 4 + 3; FALLTHROUGH;
      case 34: CHECK_FINISHED(34); COUNT_UNSET_BIT(1); m += prime * 2 + 1; FALLTHROUGH;
      case 35: CHECK_FINISHED(35); COUNT_UNSET_BIT(2); m += prime * 4 + 2; FALLTHROUGH;
      case 36: CHECK_FINISHED(36); COUNT_UNSET_BIT(5); m += prime * 2 + 1; FALLTHROUGH;
      case 37: CHECK_FINISHED(37); COUNT_UNSET_BIT(6); m += prime * 4 + 3; FALLTHROUGH;
      case 38: CHECK_FINISHED(38); COUNT_UNSET_BIT(0); m += prime * 6 + 3; FALLTHROUGH;
      case 39: CHECK_FINISHED(39); COUNT_UNSET_BIT(3); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 40: CHECK_FINISHED(40); COUNT_UNSET_BIT(5); m += prime * 6 + 4; FALLTHROUGH;
      case 41: CHECK_FINISHED(41); COUNT_UNSET_BIT(3); m += prime * 4 + 2; FALLTHROUGH;
      case 42: CHECK_FINISHED(42); COUNT_UNSET_BIT(7); m += prime * 2 + 2; FALLTHROUGH;
      case 43: CHECK_FINISHED(43); COUNT_UNSET_BIT(1); m += prime * 4 + 2; FALLTHROUGH;
      case 44: CHECK_FINISHED(44); COUNT_UNSET_BIT(6); m += prime * 2 + 2; FALLTHROUGH;
      case 45: CHECK_FINISHED(45); COUNT_UNSET_BIT(0); m += prime * 4 + 2; FALLTHROUGH;
      case 46: CHECK_FINISHED(46); COUNT_UNSET_BIT(4); m += prime * 6 + 4; FALLTHROUGH;
      case 47: CHECK_FINISHED(47); COUNT_UNSET_BIT(2); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 48: CHECK_FINISHED(48); COUNT_UNSET_BIT(6); m += prime * 6 + 5; FALLTHROUGH;
      case 49: CHECK_FINISHED(49); COUNT_UNSET_BIT(2); m += prime * 4 + 3; FALLTHROUGH;
      case 50: CHECK_FINISHED(50); COUNT_UNSET_BIT(3); m += prime * 2 + 1; FALLTHROUGH;
      case 51: CHECK_FINISHED(51); COUNT_UNSET_BIT(7); m += prime * 4 + 4; FALLTHROUGH;
      case 52: CHECK_FINISHED(52); COUNT_UNSET_BIT(0); m += prime * 2 + 1; FALLTHROUGH;
      case 53: CHECK_FINISHED(53); COUNT_UNSET_BIT(4); m += prime * 4 + 3; FALLTHROUGH;
      case 54: CHECK_FINISHED(54); COUNT_UNSET_BIT(5); m += prime * 6 + 5; FALLTHROUGH;
      case 55: CHECK_FINISHED(55); COUNT_UNSET_BIT(1); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 56: CHECK_FINISHED(56); COUNT_UNSET_BIT(7); m += prime * 6 + 6; FALLTHROUGH;
      case 57: CHECK_FINISHED(57); COUNT_UNSET_BIT(6); m += prime * 4 + 4; FALLTHROUGH;
      case 58: CHECK_FINISHED(58); COUNT_UNSET_BIT(5); m += prime * 2 + 2; FALLTHROUGH;
      case 59: CHECK_FINISHED(59); COUNT_UNSET_BIT(4); m += prime * 4 + 4; FALLTHROUGH;
      case 60: CHECK_FINISHED(60); COUNT_UNSET_BIT(3); m += prime * 2 + 2; FALLTHROUGH;
      case 61: CHECK_FINISHED(61); COUNT_UNSET_BIT(2); m += prime * 4 + 4; FALLTHROUGH;
      case 62: CHECK_FINISHED(62); COUNT_UNSET_BIT(1); m += prime * 6 + 6; FALLTHROUGH;
      case 63: CHECK_FINISHED(63); COUNT_UNSET_BIT(0); m += prime * 2 + 1;
    }

    default: UNREACHABLE;
  }
}

} // namespace

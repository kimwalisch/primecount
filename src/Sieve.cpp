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
#include "primecount-internal.hpp"

#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>

namespace primecount {

Sieve::Sieve(uint64_t low,
             uint64_t limit,
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
  allocate_counter(low, limit);
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
void Sieve::allocate_counter(uint64_t low, uint64_t limit)
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

  // Very short counter intervals cause many branch
  // mispredictions in the counting phase and hurt
  // performance. Longer intervals reduce mispredictions
  // but increase counting work, so we choose a minimum
  // interval size based on benchmark data. For smaller
  // computations (x < 2^64) larger minimum intervals are
  // usually faster; for larger computations (x > 2^64)
  // smaller intervals tend to be faster.
  double log10_limit = std::log10(max(limit, 10));
  double exponent = 17.75 - 0.94 * log10_limit;
  int e = (int) std::round(exponent);
  int min_bytes_factor = 1 << in_between(4, e, 7);
  ASSERT(min_bytes_factor >= 16);
  ASSERT(min_bytes_factor <= 128);
  uint64_t bytes = counter_.dist / 30;
  bytes = max(bytes, bytes_count_instruction * min_bytes_factor);
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
  if_unlikely(i > primeState_.size())
    primeState_.resize(i);

  // Find first multiple > start_
  ASSERT(start_ % 30 == 0);
  uint64_t quotient = start_ / prime + 1;
  uint64_t multiple = prime * quotient;

  // Find next multiple of prime that
  // is not divisible by 2, 3, 5.
  multiple += prime * wheel_init_mul[quotient % 30];

  ASSERT(multiple % 2 != 0);
  ASSERT(multiple % 3 != 0);
  ASSERT(multiple % 5 != 0);

  multiple = (multiple - start_) / 30;
  uint32_t multiple32 = (uint32_t) multiple;
  uint8_t wheel_index = wheel_indexes[quotient % 30];
  uint8_t wheel_group = wheel_groups[prime % 30];
  primeState_.push_back({multiple32, wheel_group, wheel_index});
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
  ASSERT(primeState.wheel_group <= 7);
  ASSERT(primeState.wheel_index <= 7);
  uint64_t m = primeState.multiple;
  uint64_t g = primeState.wheel_group;
  uint64_t w = primeState.wheel_index;
  uint64_t sieve_size = sieve_.size();
  uint8_t* sieve = &sieve_[0];
  const uint8_t* bitmasks = wheel_bitmasks[g];

  const Array<uint64_t, 8> adv =
  {
    prime * wheel_mul[0] + wheel_corr[g][0],
    prime * wheel_mul[1] + wheel_corr[g][1],
    prime * wheel_mul[2] + wheel_corr[g][2],
    prime * wheel_mul[3] + wheel_corr[g][3],
    prime * wheel_mul[4] + wheel_corr[g][4],
    prime * wheel_mul[5] + wheel_corr[g][5],
    prime * wheel_mul[6] + wheel_corr[g][6],
    prime * wheel_mul[7] + wheel_corr[g][7]
  };

  #define CHECK_FINISHED(i) \
    if_unlikely(m >= sieve_size) \
    { \
      primeState.wheel_index = uint8_t(i); \
      primeState.multiple = uint32_t(m - sieve_size); \
      return; \
    }

  // Get ready for loop unrolling
  for (; w; w = (w + 1) & 7)
  {
    CHECK_FINISHED(w);
    sieve[m] &= bitmasks[w];
    m += adv[w];
  }

  while (true)
  {
    CHECK_FINISHED(0); sieve[m] &= bitmasks[0]; m += adv[0];
    CHECK_FINISHED(1); sieve[m] &= bitmasks[1]; m += adv[1];
    CHECK_FINISHED(2); sieve[m] &= bitmasks[2]; m += adv[2];
    CHECK_FINISHED(3); sieve[m] &= bitmasks[3]; m += adv[3];
    CHECK_FINISHED(4); sieve[m] &= bitmasks[4]; m += adv[4];
    CHECK_FINISHED(5); sieve[m] &= bitmasks[5]; m += adv[5];
    CHECK_FINISHED(6); sieve[m] &= bitmasks[6]; m += adv[6];
    CHECK_FINISHED(7); sieve[m] &= bitmasks[7]; m += adv[7];
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
  prime /= 30;
  PrimeState& primeState = primeState_[i];
  ASSERT(primeState.wheel_group <= 7);
  ASSERT(primeState.wheel_index <= 7);
  uint64_t m = primeState.multiple;
  uint64_t g = primeState.wheel_group;
  uint64_t w = primeState.wheel_index;

  uint64_t sieve_size = sieve_.size();
  uint64_t total_count = total_count_;
  uint64_t counter_log2_dist = counter_.log2_dist;
  uint8_t* sieve = &sieve_[0];
  uint32_t* counter = &counter_[0];
  const uint8_t* bitmasks = wheel_bitmasks[g];

  const Array<uint64_t, 8> adv =
  {
    prime * wheel_mul[0] + wheel_corr[g][0],
    prime * wheel_mul[1] + wheel_corr[g][1],
    prime * wheel_mul[2] + wheel_corr[g][2],
    prime * wheel_mul[3] + wheel_corr[g][3],
    prime * wheel_mul[4] + wheel_corr[g][4],
    prime * wheel_mul[5] + wheel_corr[g][5],
    prime * wheel_mul[6] + wheel_corr[g][6],
    prime * wheel_mul[7] + wheel_corr[g][7]
  };

  #define CHECK_FINISHED(i) \
    if_unlikely(m >= sieve_size) \
    { \
      primeState.wheel_index = uint8_t(i); \
      primeState.multiple = uint32_t(m - sieve_size); \
      total_count_ = total_count; \
      return; \
    }

  #define COUNT_UNSET_BIT(i) \
    { \
      uint8_t b = sieve[m]; \
      bool is_bit = (b & bitmasks[i]) != b; \
      sieve[m] = uint8_t(b & bitmasks[i]); \
      counter[m >> counter_log2_dist] -= is_bit; \
      total_count -= is_bit; \
      m += adv[i]; \
    }

  // Get ready for loop unrolling
  for (; w; w = (w + 1) & 7)
  {
    CHECK_FINISHED(w);
    COUNT_UNSET_BIT(w);
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

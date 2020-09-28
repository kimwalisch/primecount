///
/// @file  Sieve.cpp
/// @brief The Sieve class is a highly optimized sieve of
///        Eratosthenes implementation with 30 numbers per byte
///        i.e. the 8 bits of each byte correspond to the offsets
///        { 1, 7, 11, 13, 17, 19, 23, 29 }. This Sieve also
///        skips multiples of 2, 3, 5 using wheel factorization.
///
///        Unlike a traditional prime sieve this sieve is
///        designed for use in the combinatorial prime counting
///        algorithms: this sieve removes primes as well as
///        multiples of primes and it counts the number of
///        elements that have been crossed off for the first
///        time in the sieve array.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Sieve.hpp>
#include <SieveTables.hpp>
#include <imath.hpp>
#include <min.hpp>
#include <popcnt.hpp>
#include <unlikely.hpp>

#include <stdint.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

struct WheelInit
{
  uint8_t factor;
  uint8_t index;
};

/// Categorize sieving primes according to their modulo 30
/// congruence class { 1, 7, 11, 13, 17, 19, 23, 29 }.
///
const array<int, 30> wheel_offsets =
{
  0, 8 * 0, 0, 0, 0, 0,
  0, 8 * 1, 0, 0, 0, 8 * 2,
  0, 8 * 3, 0, 0, 0, 8 * 4,
  0, 8 * 5, 0, 0, 0, 8 * 6,
  0, 0,     0, 0, 0, 8 * 7
};

/// Used to calculate the first multiple > start of a
/// sieving prime that is coprime to 2, 3, 5.
///
const WheelInit wheel_init[30] =
{
  {1,  0}, {0,  0}, {5,  1}, {4,  1}, {3,  1},
  {2,  1}, {1,  1}, {0,  1}, {3,  2}, {2,  2},
  {1,  2}, {0,  2}, {1,  3}, {0,  3}, {3,  4},
  {2,  4}, {1,  4}, {0,  4}, {1,  5}, {0,  5},
  {3,  6}, {2,  6}, {1,  6}, {0,  6}, {5,  7},
  {4,  7}, {3,  7}, {2,  7}, {1,  7}, {0,  7}
};

} // namespace

namespace primecount {

Sieve::Sieve(uint64_t low,
             uint64_t segment_size, 
             uint64_t wheel_size)
{
  assert(low % 30 == 0);
  assert(segment_size % 240 == 0);

  start_ = low;
  segment_size = get_segment_size(segment_size);

  // sieve_size = segment_size / 30 as each byte corresponds
  // to 30 numbers i.e. the 8 bits correspond to the
  // offsets = {1, 7, 11, 13, 17, 19, 23, 29}.
  sieve_.resize(segment_size / 30);

  wheel_.reserve(wheel_size);
  wheel_.resize(4);
  allocate_counters(low);
}

/// Each element of the counters array contains the current
/// number of unsieved elements in the interval:
/// [i * counters_dist, (i + 1) * counters_dist[.
/// Ideally each element of the counters array should
/// represent an interval of size O(sqrt(average_leaf_dist)).
/// Also the counter distance should be adjusted regularly
/// whilst sieving as the distance between consecutive
/// leaves is very small ~ log(x) at the beginning of the
/// sieving algorithm but grows up to segment_size towards
/// the end of the algorithm.
///
void Sieve::allocate_counters(uint64_t low)
{
  double average_leaf_dist = sqrt((double) low);
  double counters_dist = sqrt(average_leaf_dist);

  // Here we balance counting with the counters array and
  // counting from the sieve array using the POPCNT
  // instruction. Since the POPCNT instructions allows to
  // count a distance of 240 using a single instruction we
  // slightly increase the counter distance and slightly
  // decrease the size of the counters array.
  double bits_sizet = numeric_limits<size_t>::digits;
  double popcnt_dist = (bits_sizet / 8) * 30;
  counters_dist_ = (uint64_t) (counters_dist * sqrt(popcnt_dist));

  // Each byte represents an interval of size 30
  uint64_t byte_dist = counters_dist_ / 30;
  byte_dist = max(byte_dist, 64);
  byte_dist = next_power_of_2(byte_dist);
  counters_dist_ = byte_dist * 30;
  counters_dist_log2_ = ilog2(byte_dist);

  uint64_t counters_size = ceil_div(sieve_.size(), byte_dist);
  counters_.resize(counters_size);
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
uint64_t Sieve::get_segment_size(uint64_t size)
{
  size = max<uint64_t>(240, size);

  if (size % 240)
    size += 240 - size % 240;

  return size;
}

void Sieve::reset_sieve(uint64_t low, uint64_t high)
{
  fill_n(sieve_.data(), sieve_.size(), 0xff);
  uint64_t size = high - low;

  if (size < segment_size())
  {
    uint64_t last = size - 1;
    size = get_segment_size(size);
    sieve_.resize(size / 30);
    auto sieve64 = (uint64_t*) sieve_.data();
    sieve64[last / 240] &= unset_larger[last % 240];
  }
}

void Sieve::reset_counters()
{
  prev_stop_ = 0;
  count_ = 0;
  counters_i_ = 0;
  counters_count_ = 0;
  counters_stop_ = counters_dist_;
}

void Sieve::init_counters(uint64_t low, uint64_t high)
{
  reset_counters();
  total_count_ = 0;

  uint64_t start = 0;
  uint64_t max_stop = (high - 1) - low;

  while (start <= max_stop)
  {
    uint64_t stop = start + counters_dist_ - 1;
    stop = min(stop, max_stop);
    uint64_t cnt = count(start, stop);
    uint64_t byte_index = start / 30;

    counters_[byte_index >> counters_dist_log2_] = cnt;

    total_count_ += cnt;
    start += counters_dist_;
  }
}

/// Count 1 bits inside [0, stop]
uint64_t Sieve::count(uint64_t stop)
{
  assert(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

  // Quickly count the number of unsieved elements (in
  // the sieve array) up to a value that is close to
  // the stop number i.e. (stop - value) <= counters_dist.
  // We do this using the counters array, each element
  // of the counters array contains the number of
  // unsieved elements in the interval:
  // [i * counters_dist, (i + 1) * counters_dist[.
  while (counters_stop_ <= stop)
  {
    start = counters_stop_;
    counters_stop_ += counters_dist_;
    counters_count_ += counters_[counters_i_++];
    count_ = counters_count_;
  }

  // Here the remaining distance is relatively small i.e.
  // (stop - start) <= counters_dist, hence we simply
  // count the remaining number of unsieved elements by
  // linearly iterating over the sieve array.
  count_ += count(start, stop);
  return count_;
}

/// Count 1 bits inside [start, stop]
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  assert(stop - start < segment_size());

  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];
  auto sieve64 = (uint64_t*) sieve_.data();

  if (start_idx == stop_idx)
    return popcnt64(sieve64[start_idx] & (m1 & m2));
  else
  {
    uint64_t cnt = popcnt64(sieve64[start_idx] & m1);
    for (uint64_t i = start_idx + 1; i < stop_idx; i++)
      cnt += popcnt64(sieve64[i]);
    cnt += popcnt64(sieve64[stop_idx] & m2);
    return cnt;
  }
}

/// Add a sieving prime to the sieve.
/// Calculates the first multiple > start of prime that
/// is not divisible by 2, 3, 5 and its wheel index.
///
void Sieve::add(uint64_t prime)
{
  assert(start_ % 30 == 0);

  // first multiple > start_
  uint64_t quotient = start_ / prime + 1;
  uint64_t multiple = prime * quotient;

  // find next multiple of prime that
  // is not divisible by 2, 3, 5
  uint64_t factor = wheel_init[quotient % 30].factor;
  multiple += prime * factor;
  multiple = (multiple - start_) / 30;
  uint32_t multiple32 = (uint32_t) multiple;

  // calculate wheel index of multiple
  uint32_t index = wheel_init[quotient % 30].index;
  index += wheel_offsets[prime % 30];
  wheel_.emplace_back(multiple32, index);
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Used for pre-sieving.
///
void Sieve::cross_off(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime);

  prime /= 30;
  Wheel& wheel = wheel_[i];
  uint64_t m = wheel.multiple;
  uint8_t* sieve = sieve_.data();
  uint64_t sieve_size = sieve_.size();

  #define check_finished(wheel_index) \
    if_unlikely(m >= sieve_size) \
    { \
      wheel.index = wheel_index; \
      wheel.multiple = (uint32_t) (m - sieve_size); \
      return; \
    }

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: check_finished(0);
      sieve[m] &= ~(1 << 0); m += prime * 6 + 0;
      case 1: check_finished(1);
      sieve[m] &= ~(1 << 1); m += prime * 4 + 0;
      case 2: check_finished(2);
      sieve[m] &= ~(1 << 2); m += prime * 2 + 0;
      case 3: check_finished(3);
      sieve[m] &= ~(1 << 3); m += prime * 4 + 0;
      case 4: check_finished(4);
      sieve[m] &= ~(1 << 4); m += prime * 2 + 0;
      case 5: check_finished(5);
      sieve[m] &= ~(1 << 5); m += prime * 4 + 0;
      case 6: check_finished(6);
      sieve[m] &= ~(1 << 6); m += prime * 6 + 0;
      case 7: check_finished(7);
      sieve[m] &= ~(1 << 7); m += prime * 2 + 1;

      while (m + prime * 28 < sieve_size)
      {
        sieve[m + prime *  0] &= ~(1 << 0);
        sieve[m + prime *  6] &= ~(1 << 1);
        sieve[m + prime * 10] &= ~(1 << 2);
        sieve[m + prime * 12] &= ~(1 << 3);
        sieve[m + prime * 16] &= ~(1 << 4);
        sieve[m + prime * 18] &= ~(1 << 5);
        sieve[m + prime * 22] &= ~(1 << 6);
        sieve[m + prime * 28] &= ~(1 << 7);
        m += prime * 30 + 1;
      }
    }

    for (;;)
    {
      case  8: check_finished( 8);
      sieve[m] &= ~(1 << 1); m += prime * 6 + 1;
      case  9: check_finished( 9);
      sieve[m] &= ~(1 << 5); m += prime * 4 + 1;
      case 10: check_finished(10);
      sieve[m] &= ~(1 << 4); m += prime * 2 + 1;
      case 11: check_finished(11);
      sieve[m] &= ~(1 << 0); m += prime * 4 + 0;
      case 12: check_finished(12);
      sieve[m] &= ~(1 << 7); m += prime * 2 + 1;
      case 13: check_finished(13);
      sieve[m] &= ~(1 << 3); m += prime * 4 + 1;
      case 14: check_finished(14);
      sieve[m] &= ~(1 << 2); m += prime * 6 + 1;
      case 15: check_finished(15);
      sieve[m] &= ~(1 << 6); m += prime * 2 + 1;

      while (m + prime * 28 + 6 < sieve_size)
      {
        sieve[m + prime *  0 + 0] &= ~(1 << 1);
        sieve[m + prime *  6 + 1] &= ~(1 << 5);
        sieve[m + prime * 10 + 2] &= ~(1 << 4);
        sieve[m + prime * 12 + 3] &= ~(1 << 0);
        sieve[m + prime * 16 + 3] &= ~(1 << 7);
        sieve[m + prime * 18 + 4] &= ~(1 << 3);
        sieve[m + prime * 22 + 5] &= ~(1 << 2);
        sieve[m + prime * 28 + 6] &= ~(1 << 6);
        m += prime * 30 + 7;
      }
    }

    for (;;)
    {
      case 16: check_finished(16);
      sieve[m] &= ~(1 << 2); m += prime * 6 + 2;
      case 17: check_finished(17);
      sieve[m] &= ~(1 << 4); m += prime * 4 + 2;
      case 18: check_finished(18);
      sieve[m] &= ~(1 << 0); m += prime * 2 + 0;
      case 19: check_finished(19);
      sieve[m] &= ~(1 << 6); m += prime * 4 + 2;
      case 20: check_finished(20);
      sieve[m] &= ~(1 << 1); m += prime * 2 + 0;
      case 21: check_finished(21);
      sieve[m] &= ~(1 << 7); m += prime * 4 + 2;
      case 22: check_finished(22);
      sieve[m] &= ~(1 << 3); m += prime * 6 + 2;
      case 23: check_finished(23);
      sieve[m] &= ~(1 << 5); m += prime * 2 + 1;

      while (m + prime * 28 + 10 < sieve_size)
      {
        sieve[m + prime *  0 +  0] &= ~(1 << 2);
        sieve[m + prime *  6 +  2] &= ~(1 << 4);
        sieve[m + prime * 10 +  4] &= ~(1 << 0);
        sieve[m + prime * 12 +  4] &= ~(1 << 6);
        sieve[m + prime * 16 +  6] &= ~(1 << 1);
        sieve[m + prime * 18 +  6] &= ~(1 << 7);
        sieve[m + prime * 22 +  8] &= ~(1 << 3);
        sieve[m + prime * 28 + 10] &= ~(1 << 5);
        m += prime * 30 + 11;
      }
    }

    for (;;)
    {
      case 24: check_finished(24);
      sieve[m] &= ~(1 << 3); m += prime * 6 + 3;
      case 25: check_finished(25);
      sieve[m] &= ~(1 << 0); m += prime * 4 + 1;
      case 26: check_finished(26);
      sieve[m] &= ~(1 << 6); m += prime * 2 + 1;
      case 27: check_finished(27);
      sieve[m] &= ~(1 << 5); m += prime * 4 + 2;
      case 28: check_finished(28);
      sieve[m] &= ~(1 << 2); m += prime * 2 + 1;
      case 29: check_finished(29);
      sieve[m] &= ~(1 << 1); m += prime * 4 + 1;
      case 30: check_finished(30);
      sieve[m] &= ~(1 << 7); m += prime * 6 + 3;
      case 31: check_finished(31);
      sieve[m] &= ~(1 << 4); m += prime * 2 + 1;

      while (m + prime * 28 + 12 < sieve_size)
      {
        sieve[m + prime *  0 +  0] &= ~(1 << 3);
        sieve[m + prime *  6 +  3] &= ~(1 << 0);
        sieve[m + prime * 10 +  4] &= ~(1 << 6);
        sieve[m + prime * 12 +  5] &= ~(1 << 5);
        sieve[m + prime * 16 +  7] &= ~(1 << 2);
        sieve[m + prime * 18 +  8] &= ~(1 << 1);
        sieve[m + prime * 22 +  9] &= ~(1 << 7);
        sieve[m + prime * 28 + 12] &= ~(1 << 4);
        m += prime * 30 + 13;
      }
    }

    for (;;)
    {
      case 32: check_finished(32);
      sieve[m] &= ~(1 << 4); m += prime * 6 + 3;
      case 33: check_finished(33);
      sieve[m] &= ~(1 << 7); m += prime * 4 + 3;
      case 34: check_finished(34);
      sieve[m] &= ~(1 << 1); m += prime * 2 + 1;
      case 35: check_finished(35);
      sieve[m] &= ~(1 << 2); m += prime * 4 + 2;
      case 36: check_finished(36);
      sieve[m] &= ~(1 << 5); m += prime * 2 + 1;
      case 37: check_finished(37);
      sieve[m] &= ~(1 << 6); m += prime * 4 + 3;
      case 38: check_finished(38);
      sieve[m] &= ~(1 << 0); m += prime * 6 + 3;
      case 39: check_finished(39);
      sieve[m] &= ~(1 << 3); m += prime * 2 + 1;

      while (m + prime * 28 + 16 < sieve_size)
      {
        sieve[m + prime *  0 +  0] &= ~(1 << 4);
        sieve[m + prime *  6 +  3] &= ~(1 << 7);
        sieve[m + prime * 10 +  6] &= ~(1 << 1);
        sieve[m + prime * 12 +  7] &= ~(1 << 2);
        sieve[m + prime * 16 +  9] &= ~(1 << 5);
        sieve[m + prime * 18 + 10] &= ~(1 << 6);
        sieve[m + prime * 22 + 13] &= ~(1 << 0);
        sieve[m + prime * 28 + 16] &= ~(1 << 3);
        m += prime * 30 + 17;
      }
    }

    for (;;)
    {
      case 40: check_finished(40);
      sieve[m] &= ~(1 << 5); m += prime * 6 + 4;
      case 41: check_finished(41);
      sieve[m] &= ~(1 << 3); m += prime * 4 + 2;
      case 42: check_finished(42);
      sieve[m] &= ~(1 << 7); m += prime * 2 + 2;
      case 43: check_finished(43);
      sieve[m] &= ~(1 << 1); m += prime * 4 + 2;
      case 44: check_finished(44);
      sieve[m] &= ~(1 << 6); m += prime * 2 + 2;
      case 45: check_finished(45);
      sieve[m] &= ~(1 << 0); m += prime * 4 + 2;
      case 46: check_finished(46);
      sieve[m] &= ~(1 << 4); m += prime * 6 + 4;
      case 47: check_finished(47);
      sieve[m] &= ~(1 << 2); m += prime * 2 + 1;

      while (m + prime * 28 + 18 < sieve_size)
      {
        sieve[m + prime *  0 +  0] &= ~(1 << 5);
        sieve[m + prime *  6 +  4] &= ~(1 << 3);
        sieve[m + prime * 10 +  6] &= ~(1 << 7);
        sieve[m + prime * 12 +  8] &= ~(1 << 1);
        sieve[m + prime * 16 + 10] &= ~(1 << 6);
        sieve[m + prime * 18 + 12] &= ~(1 << 0);
        sieve[m + prime * 22 + 14] &= ~(1 << 4);
        sieve[m + prime * 28 + 18] &= ~(1 << 2);
        m += prime * 30 + 19;
      }
    }

    for (;;)
    {
      case 48: check_finished(48);
      sieve[m] &= ~(1 << 6); m += prime * 6 + 5;
      case 49: check_finished(49);
      sieve[m] &= ~(1 << 2); m += prime * 4 + 3;
      case 50: check_finished(50);
      sieve[m] &= ~(1 << 3); m += prime * 2 + 1;
      case 51: check_finished(51);
      sieve[m] &= ~(1 << 7); m += prime * 4 + 4;
      case 52: check_finished(52);
      sieve[m] &= ~(1 << 0); m += prime * 2 + 1;
      case 53: check_finished(53);
      sieve[m] &= ~(1 << 4); m += prime * 4 + 3;
      case 54: check_finished(54);
      sieve[m] &= ~(1 << 5); m += prime * 6 + 5;
      case 55: check_finished(55);
      sieve[m] &= ~(1 << 1); m += prime * 2 + 1;

      while (m + prime * 28 + 22 < sieve_size)
      {
        sieve[m + prime *  0 +  0] &= ~(1 << 6);
        sieve[m + prime *  6 +  5] &= ~(1 << 2);
        sieve[m + prime * 10 +  8] &= ~(1 << 3);
        sieve[m + prime * 12 +  9] &= ~(1 << 7);
        sieve[m + prime * 16 + 13] &= ~(1 << 0);
        sieve[m + prime * 18 + 14] &= ~(1 << 4);
        sieve[m + prime * 22 + 17] &= ~(1 << 5);
        sieve[m + prime * 28 + 22] &= ~(1 << 1);
        m += prime * 30 + 23;
      }
    }

    for (;;)
    {
      case 56: check_finished(56);
      sieve[m] &= ~(1 << 7); m += prime * 6 + 6;
      case 57: check_finished(57);
      sieve[m] &= ~(1 << 6); m += prime * 4 + 4;
      case 58: check_finished(58);
      sieve[m] &= ~(1 << 5); m += prime * 2 + 2;
      case 59: check_finished(59);
      sieve[m] &= ~(1 << 4); m += prime * 4 + 4;
      case 60: check_finished(60);
      sieve[m] &= ~(1 << 3); m += prime * 2 + 2;
      case 61: check_finished(61);
      sieve[m] &= ~(1 << 2); m += prime * 4 + 4;
      case 62: check_finished(62);
      sieve[m] &= ~(1 << 1); m += prime * 6 + 6;
      case 63: check_finished(63);
      sieve[m] &= ~(1 << 0); m += prime * 2 + 1;

      while (m + prime * 28 + 28 < sieve_size)
      {
        sieve[m + prime *  0 +  0] &= ~(1 << 7);
        sieve[m + prime *  6 +  6] &= ~(1 << 6);
        sieve[m + prime * 10 + 10] &= ~(1 << 5);
        sieve[m + prime * 12 + 12] &= ~(1 << 4);
        sieve[m + prime * 16 + 16] &= ~(1 << 3);
        sieve[m + prime * 18 + 18] &= ~(1 << 2);
        sieve[m + prime * 22 + 22] &= ~(1 << 1);
        sieve[m + prime * 28 + 28] &= ~(1 << 0);
        m += prime * 30 + 29;
      }
    }
  }

  #undef check_finished
}

/// Remove the i-th prime and the multiples of the i-th prime
/// from the sieve array. Also counts the number of elements
/// removed for the first time i.e. the count of sieved elements
/// whose least prime factor is the i-th prime.
///
void Sieve::cross_off_count(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime);

  reset_counters();
  Wheel& wheel = wheel_[i];
  prime /= 30;

  uint64_t m = wheel.multiple;
  uint64_t total_count = total_count_;
  uint64_t counters_dist_log2 = counters_dist_log2_;
  uint64_t sieve_size = sieve_.size();
  uint64_t* counters = counters_.data();
  uint8_t* sieve = sieve_.data();

  #define count_and_unset_bit(bit_index) \
  { \
    auto is_bit = (sieve[m] >> bit_index) & 1; \
    counters[m >> counters_dist_log2] -= is_bit; \
    total_count -= is_bit; \
    sieve[m] &= ~(1 << bit_index); \
  }

  #define check_finished(wheel_index) \
    if_unlikely(m >= sieve_size) \
    { \
      wheel.index = wheel_index; \
      wheel.multiple = (uint32_t) (m - sieve_size); \
      total_count_ = total_count; \
      return; \
    }

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: check_finished(0);
      count_and_unset_bit(0); m += prime * 6 + 0;
      case 1: check_finished(1);
      count_and_unset_bit(1); m += prime * 4 + 0;
      case 2: check_finished(2);
      count_and_unset_bit(2); m += prime * 2 + 0;
      case 3: check_finished(3);
      count_and_unset_bit(3); m += prime * 4 + 0;
      case 4: check_finished(4);
      count_and_unset_bit(4); m += prime * 2 + 0;
      case 5: check_finished(5);
      count_and_unset_bit(5); m += prime * 4 + 0;
      case 6: check_finished(6);
      count_and_unset_bit(6); m += prime * 6 + 0;
      case 7: check_finished(7);
      count_and_unset_bit(7); m += prime * 2 + 1;
    }

    for (;;)
    {
      case  8: check_finished(8);
      count_and_unset_bit(1); m += prime * 6 + 1;
      case  9: check_finished(9);
      count_and_unset_bit(5); m += prime * 4 + 1;
      case 10: check_finished(10);
      count_and_unset_bit(4); m += prime * 2 + 1;
      case 11: check_finished(11);
      count_and_unset_bit(0); m += prime * 4 + 0;
      case 12: check_finished(12);
      count_and_unset_bit(7); m += prime * 2 + 1;
      case 13: check_finished(13);
      count_and_unset_bit(3); m += prime * 4 + 1;
      case 14: check_finished(14);
      count_and_unset_bit(2); m += prime * 6 + 1;
      case 15: check_finished(15);
      count_and_unset_bit(6); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 16: check_finished(16);
      count_and_unset_bit(2); m += prime * 6 + 2;
      case 17: check_finished(17);
      count_and_unset_bit(4); m += prime * 4 + 2;
      case 18: check_finished(18);
      count_and_unset_bit(0); m += prime * 2 + 0;
      case 19: check_finished(19);
      count_and_unset_bit(6); m += prime * 4 + 2;
      case 20: check_finished(20);
      count_and_unset_bit(1); m += prime * 2 + 0;
      case 21: check_finished(21);
      count_and_unset_bit(7); m += prime * 4 + 2;
      case 22: check_finished(22);
      count_and_unset_bit(3); m += prime * 6 + 2;
      case 23: check_finished(23);
      count_and_unset_bit(5); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 24: check_finished(24);
      count_and_unset_bit(3); m += prime * 6 + 3;
      case 25: check_finished(25);
      count_and_unset_bit(0); m += prime * 4 + 1;
      case 26: check_finished(26);
      count_and_unset_bit(6); m += prime * 2 + 1;
      case 27: check_finished(27);
      count_and_unset_bit(5); m += prime * 4 + 2;
      case 28: check_finished(28);
      count_and_unset_bit(2); m += prime * 2 + 1;
      case 29: check_finished(29);
      count_and_unset_bit(1); m += prime * 4 + 1;
      case 30: check_finished(30);
      count_and_unset_bit(7); m += prime * 6 + 3;
      case 31: check_finished(31);
      count_and_unset_bit(4); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 32: check_finished(32);
      count_and_unset_bit(4); m += prime * 6 + 3;
      case 33: check_finished(33);
      count_and_unset_bit(7); m += prime * 4 + 3;
      case 34: check_finished(34);
      count_and_unset_bit(1); m += prime * 2 + 1;
      case 35: check_finished(35);
      count_and_unset_bit(2); m += prime * 4 + 2;
      case 36: check_finished(36);
      count_and_unset_bit(5); m += prime * 2 + 1;
      case 37: check_finished(37);
      count_and_unset_bit(6); m += prime * 4 + 3;
      case 38: check_finished(38);
      count_and_unset_bit(0); m += prime * 6 + 3;
      case 39: check_finished(39);
      count_and_unset_bit(3); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 40: check_finished(40);
      count_and_unset_bit(5); m += prime * 6 + 4;
      case 41: check_finished(41);
      count_and_unset_bit(3); m += prime * 4 + 2;
      case 42: check_finished(42);
      count_and_unset_bit(7); m += prime * 2 + 2;
      case 43: check_finished(43);
      count_and_unset_bit(1); m += prime * 4 + 2;
      case 44: check_finished(44);
      count_and_unset_bit(6); m += prime * 2 + 2;
      case 45: check_finished(45);
      count_and_unset_bit(0); m += prime * 4 + 2;
      case 46: check_finished(46);
      count_and_unset_bit(4); m += prime * 6 + 4;
      case 47: check_finished(47);
      count_and_unset_bit(2); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 48: check_finished(48);
      count_and_unset_bit(6); m += prime * 6 + 5;
      case 49: check_finished(49);
      count_and_unset_bit(2); m += prime * 4 + 3;
      case 50: check_finished(50);
      count_and_unset_bit(3); m += prime * 2 + 1;
      case 51: check_finished(51);
      count_and_unset_bit(7); m += prime * 4 + 4;
      case 52: check_finished(52);
      count_and_unset_bit(0); m += prime * 2 + 1;
      case 53: check_finished(53);
      count_and_unset_bit(4); m += prime * 4 + 3;
      case 54: check_finished(54);
      count_and_unset_bit(5); m += prime * 6 + 5;
      case 55: check_finished(55);
      count_and_unset_bit(1); m += prime * 2 + 1;
    }

    for (;;)
    {
      case 56: check_finished(56);
      count_and_unset_bit(7); m += prime * 6 + 6;
      case 57: check_finished(57);
      count_and_unset_bit(6); m += prime * 4 + 4;
      case 58: check_finished(58);
      count_and_unset_bit(5); m += prime * 2 + 2;
      case 59: check_finished(59);
      count_and_unset_bit(4); m += prime * 4 + 4;
      case 60: check_finished(60);
      count_and_unset_bit(3); m += prime * 2 + 2;
      case 61: check_finished(61);
      count_and_unset_bit(2); m += prime * 4 + 4;
      case 62: check_finished(62);
      count_and_unset_bit(1); m += prime * 6 + 6;
      case 63: check_finished(63);
      count_and_unset_bit(0); m += prime * 2 + 1;
    }
  }
}

} // namespace

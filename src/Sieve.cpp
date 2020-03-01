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

#include <stdint.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <vector>

#define unset_bit(bit_index, i) \
  sieve[i] &= ~(1 << bit_index);

#define count_and_unset_bit(bit_index, i) \
  is_bit = (sieve[i] >> bit_index) & 1; \
  total_count -= is_bit; \
  counters[(i) >> counters_dist_log2] -= is_bit; \
  sieve[i] &= ~(1 << bit_index);

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
  sieve_size_ = segment_size / 30;
  sieve_ = new uint8_t[sieve_size_];
  deleter_.reset(sieve_);

  wheel_.reserve(wheel_size);
  wheel_.resize(4);
  allocate_counters(low);
}

/// Each element of the counters array contains the current
/// number of unsieved elements in the interval:
/// [i * counters_dist, (i + 1) * counters_dist[.
/// Ideally each element of the counters array should
/// represent an interval of size:
/// min(sqrt(average_leaf_dist), sqrt(segment_size)).
/// Also the counter distance should be adjusted whilst
/// sieving e.g. after each sieved segment. The distance
/// between consecutive leaves is very small (~ log(x))
/// at the beginning of the sieve algorithm but grows up
/// to segment_size towards the end of the sieve.
///
void Sieve::allocate_counters(uint64_t low)
{
  uint64_t average_leaf_dist = isqrt(low);
  counters_dist_ = isqrt(average_leaf_dist);

  // Each byte represents an interval of size 30
  uint64_t byte_dist = counters_dist_ / 30;
  byte_dist = max(byte_dist, 256);
  byte_dist = nearest_power_of_2(byte_dist);
  counters_dist_ = byte_dist * 30;
  counters_dist_log2_ = ilog2(byte_dist);

  uint64_t counters_size = ceil_div(sieve_size_, byte_dist);
  counters_.resize(counters_size);
}

/// The segment size is sieve_size * 30 as each
/// byte corresponds to 30 numbers.
///
uint64_t Sieve::segment_size() const
{
  return sieve_size_ * 30;
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
  fill_n(sieve_, sieve_size_, (uint8_t) 0xff);
  uint64_t size = high - low;

  if (size < segment_size())
  {
    uint64_t last = size - 1;
    size = get_segment_size(size);
    sieve_size_ = size / 30;
    auto sieve64 = (uint64_t*) sieve_;
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

  uint64_t max_stop = (high - 1) - low;

  for (uint64_t i = 0; i <= max_stop; i += counters_dist_)
  {
    uint64_t start = i;
    uint64_t stop = start + counters_dist_ - 1;
    stop = min(stop, max_stop);
    uint64_t cnt = count(start, stop);
    uint64_t byte_index = i / 30;

    counters_[byte_index >> counters_dist_log2_] = cnt;
    total_count_ += cnt;
  }
}

/// Count 1 bits inside [0, stop]
uint64_t Sieve::count(uint64_t stop)
{
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

  if (start > stop)
    return count_;

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
  assert(stop - start < segment_size());

  uint64_t bit_count = 0;
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];
  auto sieve64 = (uint64_t*) sieve_;

  if (start_idx == stop_idx)
    bit_count = popcnt64(sieve64[start_idx] & (m1 & m2));
  else
  {
    bit_count = popcnt64(sieve64[start_idx] & m1);
    bit_count += popcnt(&sieve64[start_idx + 1], stop_idx - (start_idx + 1));
    bit_count += popcnt64(sieve64[stop_idx] & m2);
  }

  return bit_count;
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

  Wheel& wheel = wheel_[i];
  prime /= 30;

  uint64_t m = wheel.multiple;
  uint64_t sieve_size = sieve_size_;
  uint8_t* sieve = sieve_;

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: if (m >= sieve_size) { wheel.index = 0; break; }
      unset_bit(0, m); m += prime * 6 + 0;
      case 1: if (m >= sieve_size) { wheel.index = 1; break; }
      unset_bit(1, m); m += prime * 4 + 0;
      case 2: if (m >= sieve_size) { wheel.index = 2; break; }
      unset_bit(2, m); m += prime * 2 + 0;
      case 3: if (m >= sieve_size) { wheel.index = 3; break; }
      unset_bit(3, m); m += prime * 4 + 0;
      case 4: if (m >= sieve_size) { wheel.index = 4; break; }
      unset_bit(4, m); m += prime * 2 + 0;
      case 5: if (m >= sieve_size) { wheel.index = 5; break; }
      unset_bit(5, m); m += prime * 4 + 0;
      case 6: if (m >= sieve_size) { wheel.index = 6; break; }
      unset_bit(6, m); m += prime * 6 + 0;
      case 7: if (m >= sieve_size) { wheel.index = 7; break; }
      unset_bit(7, m); m += prime * 2 + 1;

      while (m + prime * 28 < sieve_size)
      {
        unset_bit(0, m + prime *  0);
        unset_bit(1, m + prime *  6);
        unset_bit(2, m + prime * 10);
        unset_bit(3, m + prime * 12);
        unset_bit(4, m + prime * 16);
        unset_bit(5, m + prime * 18);
        unset_bit(6, m + prime * 22);
        unset_bit(7, m + prime * 28);
        m += prime * 30 + 1;
      }
    }
    break;

    for (;;)
    {
      case  8: if (m >= sieve_size) { wheel.index =  8; break; }
      unset_bit(1, m); m += prime * 6 + 1;
      case  9: if (m >= sieve_size) { wheel.index =  9; break; }
      unset_bit(5, m); m += prime * 4 + 1;
      case 10: if (m >= sieve_size) { wheel.index = 10; break; }
      unset_bit(4, m); m += prime * 2 + 1;
      case 11: if (m >= sieve_size) { wheel.index = 11; break; }
      unset_bit(0, m); m += prime * 4 + 0;
      case 12: if (m >= sieve_size) { wheel.index = 12; break; }
      unset_bit(7, m); m += prime * 2 + 1;
      case 13: if (m >= sieve_size) { wheel.index = 13; break; }
      unset_bit(3, m); m += prime * 4 + 1;
      case 14: if (m >= sieve_size) { wheel.index = 14; break; }
      unset_bit(2, m); m += prime * 6 + 1;
      case 15: if (m >= sieve_size) { wheel.index = 15; break; }
      unset_bit(6, m); m += prime * 2 + 1;

      while (m + prime * 28 + 6 < sieve_size)
      {
        unset_bit(1, m + prime *  0 + 0);
        unset_bit(5, m + prime *  6 + 1);
        unset_bit(4, m + prime * 10 + 2);
        unset_bit(0, m + prime * 12 + 3);
        unset_bit(7, m + prime * 16 + 3);
        unset_bit(3, m + prime * 18 + 4);
        unset_bit(2, m + prime * 22 + 5);
        unset_bit(6, m + prime * 28 + 6);
        m += prime * 30 + 7;
      }
    }
    break;

    for (;;)
    {
      case 16: if (m >= sieve_size) { wheel.index = 16; break; }
      unset_bit(2, m); m += prime * 6 + 2;
      case 17: if (m >= sieve_size) { wheel.index = 17; break; }
      unset_bit(4, m); m += prime * 4 + 2;
      case 18: if (m >= sieve_size) { wheel.index = 18; break; }
      unset_bit(0, m); m += prime * 2 + 0;
      case 19: if (m >= sieve_size) { wheel.index = 19; break; }
      unset_bit(6, m); m += prime * 4 + 2;
      case 20: if (m >= sieve_size) { wheel.index = 20; break; }
      unset_bit(1, m); m += prime * 2 + 0;
      case 21: if (m >= sieve_size) { wheel.index = 21; break; }
      unset_bit(7, m); m += prime * 4 + 2;
      case 22: if (m >= sieve_size) { wheel.index = 22; break; }
      unset_bit(3, m); m += prime * 6 + 2;
      case 23: if (m >= sieve_size) { wheel.index = 23; break; }
      unset_bit(5, m); m += prime * 2 + 1;

      while (m + prime * 28 + 10 < sieve_size)
      {
        unset_bit(2, m + prime *  0 +  0);
        unset_bit(4, m + prime *  6 +  2);
        unset_bit(0, m + prime * 10 +  4);
        unset_bit(6, m + prime * 12 +  4);
        unset_bit(1, m + prime * 16 +  6);
        unset_bit(7, m + prime * 18 +  6);
        unset_bit(3, m + prime * 22 +  8);
        unset_bit(5, m + prime * 28 + 10);
        m += prime * 30 + 11;
      }
    }
    break;

    for (;;)
    {
      case 24: if (m >= sieve_size) { wheel.index = 24; break; }
      unset_bit(3, m); m += prime * 6 + 3;
      case 25: if (m >= sieve_size) { wheel.index = 25; break; }
      unset_bit(0, m); m += prime * 4 + 1;
      case 26: if (m >= sieve_size) { wheel.index = 26; break; }
      unset_bit(6, m); m += prime * 2 + 1;
      case 27: if (m >= sieve_size) { wheel.index = 27; break; }
      unset_bit(5, m); m += prime * 4 + 2;
      case 28: if (m >= sieve_size) { wheel.index = 28; break; }
      unset_bit(2, m); m += prime * 2 + 1;
      case 29: if (m >= sieve_size) { wheel.index = 29; break; }
      unset_bit(1, m); m += prime * 4 + 1;
      case 30: if (m >= sieve_size) { wheel.index = 30; break; }
      unset_bit(7, m); m += prime * 6 + 3;
      case 31: if (m >= sieve_size) { wheel.index = 31; break; }
      unset_bit(4, m); m += prime * 2 + 1;

      while (m + prime * 28 + 12 < sieve_size)
      {
        unset_bit(3, m + prime *  0 +  0);
        unset_bit(0, m + prime *  6 +  3);
        unset_bit(6, m + prime * 10 +  4);
        unset_bit(5, m + prime * 12 +  5);
        unset_bit(2, m + prime * 16 +  7);
        unset_bit(1, m + prime * 18 +  8);
        unset_bit(7, m + prime * 22 +  9);
        unset_bit(4, m + prime * 28 + 12);
        m += prime * 30 + 13;
      }
    }
    break;

    for (;;)
    {
      case 32: if (m >= sieve_size) { wheel.index = 32; break; }
      unset_bit(4, m); m += prime * 6 + 3;
      case 33: if (m >= sieve_size) { wheel.index = 33; break; }
      unset_bit(7, m); m += prime * 4 + 3;
      case 34: if (m >= sieve_size) { wheel.index = 34; break; }
      unset_bit(1, m); m += prime * 2 + 1;
      case 35: if (m >= sieve_size) { wheel.index = 35; break; }
      unset_bit(2, m); m += prime * 4 + 2;
      case 36: if (m >= sieve_size) { wheel.index = 36; break; }
      unset_bit(5, m); m += prime * 2 + 1;
      case 37: if (m >= sieve_size) { wheel.index = 37; break; }
      unset_bit(6, m); m += prime * 4 + 3;
      case 38: if (m >= sieve_size) { wheel.index = 38; break; }
      unset_bit(0, m); m += prime * 6 + 3;
      case 39: if (m >= sieve_size) { wheel.index = 39; break; }
      unset_bit(3, m); m += prime * 2 + 1;

      while (m + prime * 28 + 16 < sieve_size)
      {
        unset_bit(4, m + prime *  0 +  0);
        unset_bit(7, m + prime *  6 +  3);
        unset_bit(1, m + prime * 10 +  6);
        unset_bit(2, m + prime * 12 +  7);
        unset_bit(5, m + prime * 16 +  9);
        unset_bit(6, m + prime * 18 + 10);
        unset_bit(0, m + prime * 22 + 13);
        unset_bit(3, m + prime * 28 + 16);
        m += prime * 30 + 17;
      }
    }
    break;

    for (;;)
    {
      case 40: if (m >= sieve_size) { wheel.index = 40; break; }
      unset_bit(5, m); m += prime * 6 + 4;
      case 41: if (m >= sieve_size) { wheel.index = 41; break; }
      unset_bit(3, m); m += prime * 4 + 2;
      case 42: if (m >= sieve_size) { wheel.index = 42; break; }
      unset_bit(7, m); m += prime * 2 + 2;
      case 43: if (m >= sieve_size) { wheel.index = 43; break; }
      unset_bit(1, m); m += prime * 4 + 2;
      case 44: if (m >= sieve_size) { wheel.index = 44; break; }
      unset_bit(6, m); m += prime * 2 + 2;
      case 45: if (m >= sieve_size) { wheel.index = 45; break; }
      unset_bit(0, m); m += prime * 4 + 2;
      case 46: if (m >= sieve_size) { wheel.index = 46; break; }
      unset_bit(4, m); m += prime * 6 + 4;
      case 47: if (m >= sieve_size) { wheel.index = 47; break; }
      unset_bit(2, m); m += prime * 2 + 1;

      while (m + prime * 28 + 18 < sieve_size)
      {
        unset_bit(5, m + prime *  0 +  0);
        unset_bit(3, m + prime *  6 +  4);
        unset_bit(7, m + prime * 10 +  6);
        unset_bit(1, m + prime * 12 +  8);
        unset_bit(6, m + prime * 16 + 10);
        unset_bit(0, m + prime * 18 + 12);
        unset_bit(4, m + prime * 22 + 14);
        unset_bit(2, m + prime * 28 + 18);
        m += prime * 30 + 19;
      }
    }
    break;

    for (;;)
    {
      case 48: if (m >= sieve_size) { wheel.index = 48; break; }
      unset_bit(6, m); m += prime * 6 + 5;
      case 49: if (m >= sieve_size) { wheel.index = 49; break; }
      unset_bit(2, m); m += prime * 4 + 3;
      case 50: if (m >= sieve_size) { wheel.index = 50; break; }
      unset_bit(3, m); m += prime * 2 + 1;
      case 51: if (m >= sieve_size) { wheel.index = 51; break; }
      unset_bit(7, m); m += prime * 4 + 4;
      case 52: if (m >= sieve_size) { wheel.index = 52; break; }
      unset_bit(0, m); m += prime * 2 + 1;
      case 53: if (m >= sieve_size) { wheel.index = 53; break; }
      unset_bit(4, m); m += prime * 4 + 3;
      case 54: if (m >= sieve_size) { wheel.index = 54; break; }
      unset_bit(5, m); m += prime * 6 + 5;
      case 55: if (m >= sieve_size) { wheel.index = 55; break; }
      unset_bit(1, m); m += prime * 2 + 1;

      while (m + prime * 28 + 22 < sieve_size)
      {
        unset_bit(6, m + prime *  0 +  0);
        unset_bit(2, m + prime *  6 +  5);
        unset_bit(3, m + prime * 10 +  8);
        unset_bit(7, m + prime * 12 +  9);
        unset_bit(0, m + prime * 16 + 13);
        unset_bit(4, m + prime * 18 + 14);
        unset_bit(5, m + prime * 22 + 17);
        unset_bit(1, m + prime * 28 + 22);
        m += prime * 30 + 23;
      }
    }
    break;

    for (;;)
    {
      case 56: if (m >= sieve_size) { wheel.index = 56; break; }
      unset_bit(7, m); m += prime * 6 + 6;
      case 57: if (m >= sieve_size) { wheel.index = 57; break; }
      unset_bit(6, m); m += prime * 4 + 4;
      case 58: if (m >= sieve_size) { wheel.index = 58; break; }
      unset_bit(5, m); m += prime * 2 + 2;
      case 59: if (m >= sieve_size) { wheel.index = 59; break; }
      unset_bit(4, m); m += prime * 4 + 4;
      case 60: if (m >= sieve_size) { wheel.index = 60; break; }
      unset_bit(3, m); m += prime * 2 + 2;
      case 61: if (m >= sieve_size) { wheel.index = 61; break; }
      unset_bit(2, m); m += prime * 4 + 4;
      case 62: if (m >= sieve_size) { wheel.index = 62; break; }
      unset_bit(1, m); m += prime * 6 + 6;
      case 63: if (m >= sieve_size) { wheel.index = 63; break; }
      unset_bit(0, m); m += prime * 2 + 1;

      while (m + prime * 28 + 28 < sieve_size)
      {
        unset_bit(7, m + prime *  0 +  0);
        unset_bit(6, m + prime *  6 +  6);
        unset_bit(5, m + prime * 10 + 10);
        unset_bit(4, m + prime * 12 + 12);
        unset_bit(3, m + prime * 16 + 16);
        unset_bit(2, m + prime * 18 + 18);
        unset_bit(1, m + prime * 22 + 22);
        unset_bit(0, m + prime * 28 + 28);
        m += prime * 30 + 29;
      }
    }
    break;
  }

  // update for the next segment
  wheel.multiple = (uint32_t) (m - sieve_size);
}

/// Remove the i-th prime and the multiples of the i-th
/// prime from the sieve array. Also counts the number
/// of elements removed for the first time i.e. the
/// count of sieved elements whose least prime factor
/// is the i-th prime.
///
void Sieve::cross_off_count(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime);

  Wheel& wheel = wheel_[i];
  prime /= 30;

  uint64_t is_bit = 0;
  uint64_t total_count = total_count_;
  uint64_t m = wheel.multiple;
  uint64_t sieve_size = sieve_size_;
  uint64_t counters_dist_log2 = counters_dist_log2_;
  uint64_t* counters = counters_.data();
  uint8_t* sieve = sieve_;

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: if (m >= sieve_size) { wheel.index = 0; break; }
      count_and_unset_bit(0, m); m += prime * 6 + 0;
      case 1: if (m >= sieve_size) { wheel.index = 1; break; }
      count_and_unset_bit(1, m); m += prime * 4 + 0;
      case 2: if (m >= sieve_size) { wheel.index = 2; break; }
      count_and_unset_bit(2, m); m += prime * 2 + 0;
      case 3: if (m >= sieve_size) { wheel.index = 3; break; }
      count_and_unset_bit(3, m); m += prime * 4 + 0;
      case 4: if (m >= sieve_size) { wheel.index = 4; break; }
      count_and_unset_bit(4, m); m += prime * 2 + 0;
      case 5: if (m >= sieve_size) { wheel.index = 5; break; }
      count_and_unset_bit(5, m); m += prime * 4 + 0;
      case 6: if (m >= sieve_size) { wheel.index = 6; break; }
      count_and_unset_bit(6, m); m += prime * 6 + 0;
      case 7: if (m >= sieve_size) { wheel.index = 7; break; }
      count_and_unset_bit(7, m); m += prime * 2 + 1;

      while (m + prime * 28 < sieve_size)
      {
        count_and_unset_bit(0, m + prime *  0);
        count_and_unset_bit(1, m + prime *  6);
        count_and_unset_bit(2, m + prime * 10);
        count_and_unset_bit(3, m + prime * 12);
        count_and_unset_bit(4, m + prime * 16);
        count_and_unset_bit(5, m + prime * 18);
        count_and_unset_bit(6, m + prime * 22);
        count_and_unset_bit(7, m + prime * 28);
        m += prime * 30 + 1;
      }
    }
    break;

    for (;;)
    {
      case  8: if (m >= sieve_size) { wheel.index =  8; break; }
      count_and_unset_bit(1, m); m += prime * 6 + 1;
      case  9: if (m >= sieve_size) { wheel.index =  9; break; }
      count_and_unset_bit(5, m); m += prime * 4 + 1;
      case 10: if (m >= sieve_size) { wheel.index = 10; break; }
      count_and_unset_bit(4, m); m += prime * 2 + 1;
      case 11: if (m >= sieve_size) { wheel.index = 11; break; }
      count_and_unset_bit(0, m); m += prime * 4 + 0;
      case 12: if (m >= sieve_size) { wheel.index = 12; break; }
      count_and_unset_bit(7, m); m += prime * 2 + 1;
      case 13: if (m >= sieve_size) { wheel.index = 13; break; }
      count_and_unset_bit(3, m); m += prime * 4 + 1;
      case 14: if (m >= sieve_size) { wheel.index = 14; break; }
      count_and_unset_bit(2, m); m += prime * 6 + 1;
      case 15: if (m >= sieve_size) { wheel.index = 15; break; }
      count_and_unset_bit(6, m); m += prime * 2 + 1;

      while (m + prime * 28 + 6 < sieve_size)
      {
        count_and_unset_bit(1, m + prime *  0 + 0);
        count_and_unset_bit(5, m + prime *  6 + 1);
        count_and_unset_bit(4, m + prime * 10 + 2);
        count_and_unset_bit(0, m + prime * 12 + 3);
        count_and_unset_bit(7, m + prime * 16 + 3);
        count_and_unset_bit(3, m + prime * 18 + 4);
        count_and_unset_bit(2, m + prime * 22 + 5);
        count_and_unset_bit(6, m + prime * 28 + 6);
        m += prime * 30 + 7;
      }
    }
    break;

    for (;;)
    {
      case 16: if (m >= sieve_size) { wheel.index = 16; break; }
      count_and_unset_bit(2, m); m += prime * 6 + 2;
      case 17: if (m >= sieve_size) { wheel.index = 17; break; }
      count_and_unset_bit(4, m); m += prime * 4 + 2;
      case 18: if (m >= sieve_size) { wheel.index = 18; break; }
      count_and_unset_bit(0, m); m += prime * 2 + 0;
      case 19: if (m >= sieve_size) { wheel.index = 19; break; }
      count_and_unset_bit(6, m); m += prime * 4 + 2;
      case 20: if (m >= sieve_size) { wheel.index = 20; break; }
      count_and_unset_bit(1, m); m += prime * 2 + 0;
      case 21: if (m >= sieve_size) { wheel.index = 21; break; }
      count_and_unset_bit(7, m); m += prime * 4 + 2;
      case 22: if (m >= sieve_size) { wheel.index = 22; break; }
      count_and_unset_bit(3, m); m += prime * 6 + 2;
      case 23: if (m >= sieve_size) { wheel.index = 23; break; }
      count_and_unset_bit(5, m); m += prime * 2 + 1;

      while (m + prime * 28 + 10 < sieve_size)
      {
        count_and_unset_bit(2, m + prime *  0 +  0);
        count_and_unset_bit(4, m + prime *  6 +  2);
        count_and_unset_bit(0, m + prime * 10 +  4);
        count_and_unset_bit(6, m + prime * 12 +  4);
        count_and_unset_bit(1, m + prime * 16 +  6);
        count_and_unset_bit(7, m + prime * 18 +  6);
        count_and_unset_bit(3, m + prime * 22 +  8);
        count_and_unset_bit(5, m + prime * 28 + 10);
        m += prime * 30 + 11;
      }
    }
    break;

    for (;;)
    {
      case 24: if (m >= sieve_size) { wheel.index = 24; break; }
      count_and_unset_bit(3, m); m += prime * 6 + 3;
      case 25: if (m >= sieve_size) { wheel.index = 25; break; }
      count_and_unset_bit(0, m); m += prime * 4 + 1;
      case 26: if (m >= sieve_size) { wheel.index = 26; break; }
      count_and_unset_bit(6, m); m += prime * 2 + 1;
      case 27: if (m >= sieve_size) { wheel.index = 27; break; }
      count_and_unset_bit(5, m); m += prime * 4 + 2;
      case 28: if (m >= sieve_size) { wheel.index = 28; break; }
      count_and_unset_bit(2, m); m += prime * 2 + 1;
      case 29: if (m >= sieve_size) { wheel.index = 29; break; }
      count_and_unset_bit(1, m); m += prime * 4 + 1;
      case 30: if (m >= sieve_size) { wheel.index = 30; break; }
      count_and_unset_bit(7, m); m += prime * 6 + 3;
      case 31: if (m >= sieve_size) { wheel.index = 31; break; }
      count_and_unset_bit(4, m); m += prime * 2 + 1;

      while (m + prime * 28 + 12 < sieve_size)
      {
        count_and_unset_bit(3, m + prime *  0 +  0);
        count_and_unset_bit(0, m + prime *  6 +  3);
        count_and_unset_bit(6, m + prime * 10 +  4);
        count_and_unset_bit(5, m + prime * 12 +  5);
        count_and_unset_bit(2, m + prime * 16 +  7);
        count_and_unset_bit(1, m + prime * 18 +  8);
        count_and_unset_bit(7, m + prime * 22 +  9);
        count_and_unset_bit(4, m + prime * 28 + 12);
        m += prime * 30 + 13;
      }
    }
    break;

    for (;;)
    {
      case 32: if (m >= sieve_size) { wheel.index = 32; break; }
      count_and_unset_bit(4, m); m += prime * 6 + 3;
      case 33: if (m >= sieve_size) { wheel.index = 33; break; }
      count_and_unset_bit(7, m); m += prime * 4 + 3;
      case 34: if (m >= sieve_size) { wheel.index = 34; break; }
      count_and_unset_bit(1, m); m += prime * 2 + 1;
      case 35: if (m >= sieve_size) { wheel.index = 35; break; }
      count_and_unset_bit(2, m); m += prime * 4 + 2;
      case 36: if (m >= sieve_size) { wheel.index = 36; break; }
      count_and_unset_bit(5, m); m += prime * 2 + 1;
      case 37: if (m >= sieve_size) { wheel.index = 37; break; }
      count_and_unset_bit(6, m); m += prime * 4 + 3;
      case 38: if (m >= sieve_size) { wheel.index = 38; break; }
      count_and_unset_bit(0, m); m += prime * 6 + 3;
      case 39: if (m >= sieve_size) { wheel.index = 39; break; }
      count_and_unset_bit(3, m); m += prime * 2 + 1;

      while (m + prime * 28 + 16 < sieve_size)
      {
        count_and_unset_bit(4, m + prime *  0 +  0);
        count_and_unset_bit(7, m + prime *  6 +  3);
        count_and_unset_bit(1, m + prime * 10 +  6);
        count_and_unset_bit(2, m + prime * 12 +  7);
        count_and_unset_bit(5, m + prime * 16 +  9);
        count_and_unset_bit(6, m + prime * 18 + 10);
        count_and_unset_bit(0, m + prime * 22 + 13);
        count_and_unset_bit(3, m + prime * 28 + 16);
        m += prime * 30 + 17;
      }
    }
    break;

    for (;;)
    {
      case 40: if (m >= sieve_size) { wheel.index = 40; break; }
      count_and_unset_bit(5, m); m += prime * 6 + 4;
      case 41: if (m >= sieve_size) { wheel.index = 41; break; }
      count_and_unset_bit(3, m); m += prime * 4 + 2;
      case 42: if (m >= sieve_size) { wheel.index = 42; break; }
      count_and_unset_bit(7, m); m += prime * 2 + 2;
      case 43: if (m >= sieve_size) { wheel.index = 43; break; }
      count_and_unset_bit(1, m); m += prime * 4 + 2;
      case 44: if (m >= sieve_size) { wheel.index = 44; break; }
      count_and_unset_bit(6, m); m += prime * 2 + 2;
      case 45: if (m >= sieve_size) { wheel.index = 45; break; }
      count_and_unset_bit(0, m); m += prime * 4 + 2;
      case 46: if (m >= sieve_size) { wheel.index = 46; break; }
      count_and_unset_bit(4, m); m += prime * 6 + 4;
      case 47: if (m >= sieve_size) { wheel.index = 47; break; }
      count_and_unset_bit(2, m); m += prime * 2 + 1;

      while (m + prime * 28 + 18 < sieve_size)
      {
        count_and_unset_bit(5, m + prime *  0 +  0);
        count_and_unset_bit(3, m + prime *  6 +  4);
        count_and_unset_bit(7, m + prime * 10 +  6);
        count_and_unset_bit(1, m + prime * 12 +  8);
        count_and_unset_bit(6, m + prime * 16 + 10);
        count_and_unset_bit(0, m + prime * 18 + 12);
        count_and_unset_bit(4, m + prime * 22 + 14);
        count_and_unset_bit(2, m + prime * 28 + 18);
        m += prime * 30 + 19;
      }
    }
    break;

    for (;;)
    {
      case 48: if (m >= sieve_size) { wheel.index = 48; break; }
      count_and_unset_bit(6, m); m += prime * 6 + 5;
      case 49: if (m >= sieve_size) { wheel.index = 49; break; }
      count_and_unset_bit(2, m); m += prime * 4 + 3;
      case 50: if (m >= sieve_size) { wheel.index = 50; break; }
      count_and_unset_bit(3, m); m += prime * 2 + 1;
      case 51: if (m >= sieve_size) { wheel.index = 51; break; }
      count_and_unset_bit(7, m); m += prime * 4 + 4;
      case 52: if (m >= sieve_size) { wheel.index = 52; break; }
      count_and_unset_bit(0, m); m += prime * 2 + 1;
      case 53: if (m >= sieve_size) { wheel.index = 53; break; }
      count_and_unset_bit(4, m); m += prime * 4 + 3;
      case 54: if (m >= sieve_size) { wheel.index = 54; break; }
      count_and_unset_bit(5, m); m += prime * 6 + 5;
      case 55: if (m >= sieve_size) { wheel.index = 55; break; }
      count_and_unset_bit(1, m); m += prime * 2 + 1;

      while (m + prime * 28 + 22 < sieve_size)
      {
        count_and_unset_bit(6, m + prime *  0 +  0);
        count_and_unset_bit(2, m + prime *  6 +  5);
        count_and_unset_bit(3, m + prime * 10 +  8);
        count_and_unset_bit(7, m + prime * 12 +  9);
        count_and_unset_bit(0, m + prime * 16 + 13);
        count_and_unset_bit(4, m + prime * 18 + 14);
        count_and_unset_bit(5, m + prime * 22 + 17);
        count_and_unset_bit(1, m + prime * 28 + 22);
        m += prime * 30 + 23;
      }
    }
    break;

    for (;;)
    {
      case 56: if (m >= sieve_size) { wheel.index = 56; break; }
      count_and_unset_bit(7, m); m += prime * 6 + 6;
      case 57: if (m >= sieve_size) { wheel.index = 57; break; }
      count_and_unset_bit(6, m); m += prime * 4 + 4;
      case 58: if (m >= sieve_size) { wheel.index = 58; break; }
      count_and_unset_bit(5, m); m += prime * 2 + 2;
      case 59: if (m >= sieve_size) { wheel.index = 59; break; }
      count_and_unset_bit(4, m); m += prime * 4 + 4;
      case 60: if (m >= sieve_size) { wheel.index = 60; break; }
      count_and_unset_bit(3, m); m += prime * 2 + 2;
      case 61: if (m >= sieve_size) { wheel.index = 61; break; }
      count_and_unset_bit(2, m); m += prime * 4 + 4;
      case 62: if (m >= sieve_size) { wheel.index = 62; break; }
      count_and_unset_bit(1, m); m += prime * 6 + 6;
      case 63: if (m >= sieve_size) { wheel.index = 63; break; }
      count_and_unset_bit(0, m); m += prime * 2 + 1;

      while (m + prime * 28 + 28 < sieve_size)
      {
        count_and_unset_bit(7, m + prime *  0 +  0);
        count_and_unset_bit(6, m + prime *  6 +  6);
        count_and_unset_bit(5, m + prime * 10 + 10);
        count_and_unset_bit(4, m + prime * 12 + 12);
        count_and_unset_bit(3, m + prime * 16 + 16);
        count_and_unset_bit(2, m + prime * 18 + 18);
        count_and_unset_bit(1, m + prime * 22 + 22);
        count_and_unset_bit(0, m + prime * 28 + 28);
        m += prime * 30 + 29;
      }
    }
    break;
  }

  // update for the next segment
  wheel.multiple = (uint32_t) (m - sieve_size);
  total_count_ = total_count;
  reset_counters();
}

} // namespace

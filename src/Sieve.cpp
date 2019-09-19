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
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Sieve.hpp>
#include <SieveTables.hpp>
#include <popcnt.hpp>

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

/// Unset the n-th bit
template <int n>
void unset_bit(byte_t* sieve)
{
  *sieve &= ~(1 << n);
}

/// Returns 1 if n-th bit is set, else 0
template <int n>
uint64_t is_bit(byte_t* sieve)
{
  return (*sieve >> n) & 1;
}

} // namespace

namespace primecount {

Sieve::Sieve(uint64_t start,
             uint64_t segment_size, 
             uint64_t wheel_size)
{
  assert(start % 30 == 0);
  assert(segment_size % 240 == 0);

  start_ = start;
  segment_size = get_segment_size(segment_size);

  // sieve_size = segment_size / 30 as each byte corresponds
  // to 30 numbers i.e. the 8 bits correspond to the
  // offsets = {1, 7, 11, 13, 17, 19, 23, 29}.
  sieve_size_ = segment_size / 30;
  sieve_ = new byte_t[sieve_size_];
  deleter_.reset(sieve_);

  wheel_.reserve(wheel_size);
  wheel_.resize(4);
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

/// Count 1 bits inside [start, stop]
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

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

void Sieve::reset_sieve(uint64_t low, uint64_t high)
{
  fill_n(sieve_, sieve_size_, (byte_t) 0xff);
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

/// Add a sieving prime to the sieve.
/// Calculates the first multiple > start_ of prime that
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
  // pointer to the byte with the 1st multiple
  byte_t* s = sieve_ + wheel.multiple;
  byte_t* sieve_end = sieve_ + sieve_size_;
  prime /= 30;

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: if (s >= sieve_end) { wheel.index = 0; break; }
      unset_bit<0>(s); s += prime * 6 + 0;
      case 1: if (s >= sieve_end) { wheel.index = 1; break; }
      unset_bit<1>(s); s += prime * 4 + 0;
      case 2: if (s >= sieve_end) { wheel.index = 2; break; }
      unset_bit<2>(s); s += prime * 2 + 0;
      case 3: if (s >= sieve_end) { wheel.index = 3; break; }
      unset_bit<3>(s); s += prime * 4 + 0;
      case 4: if (s >= sieve_end) { wheel.index = 4; break; }
      unset_bit<4>(s); s += prime * 2 + 0;
      case 5: if (s >= sieve_end) { wheel.index = 5; break; }
      unset_bit<5>(s); s += prime * 4 + 0;
      case 6: if (s >= sieve_end) { wheel.index = 6; break; }
      unset_bit<6>(s); s += prime * 6 + 0;
      case 7: if (s >= sieve_end) { wheel.index = 7; break; }
      unset_bit<7>(s); s += prime * 2 + 1;

      while (s + prime * 28 < sieve_end)
      {
        unset_bit<0>(s + prime *  0);
        unset_bit<1>(s + prime *  6);
        unset_bit<2>(s + prime * 10);
        unset_bit<3>(s + prime * 12);
        unset_bit<4>(s + prime * 16);
        unset_bit<5>(s + prime * 18);
        unset_bit<6>(s + prime * 22);
        unset_bit<7>(s + prime * 28);
        s += prime * 30 + 1;
      }
    }
    break;

    for (;;)
    {
      case  8: if (s >= sieve_end) { wheel.index =  8; break; }
      unset_bit<1>(s); s += prime * 6 + 1;
      case  9: if (s >= sieve_end) { wheel.index =  9; break; }
      unset_bit<5>(s); s += prime * 4 + 1;
      case 10: if (s >= sieve_end) { wheel.index = 10; break; }
      unset_bit<4>(s); s += prime * 2 + 1;
      case 11: if (s >= sieve_end) { wheel.index = 11; break; }
      unset_bit<0>(s); s += prime * 4 + 0;
      case 12: if (s >= sieve_end) { wheel.index = 12; break; }
      unset_bit<7>(s); s += prime * 2 + 1;
      case 13: if (s >= sieve_end) { wheel.index = 13; break; }
      unset_bit<3>(s); s += prime * 4 + 1;
      case 14: if (s >= sieve_end) { wheel.index = 14; break; }
      unset_bit<2>(s); s += prime * 6 + 1;
      case 15: if (s >= sieve_end) { wheel.index = 15; break; }
      unset_bit<6>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 6 < sieve_end)
      {
        unset_bit<1>(s + prime *  0 + 0);
        unset_bit<5>(s + prime *  6 + 1);
        unset_bit<4>(s + prime * 10 + 2);
        unset_bit<0>(s + prime * 12 + 3);
        unset_bit<7>(s + prime * 16 + 3);
        unset_bit<3>(s + prime * 18 + 4);
        unset_bit<2>(s + prime * 22 + 5);
        unset_bit<6>(s + prime * 28 + 6);
        s += prime * 30 + 7;
      }
    }
    break;

    for (;;)
    {
      case 16: if (s >= sieve_end) { wheel.index = 16; break; }
      unset_bit<2>(s); s += prime * 6 + 2;
      case 17: if (s >= sieve_end) { wheel.index = 17; break; }
      unset_bit<4>(s); s += prime * 4 + 2;
      case 18: if (s >= sieve_end) { wheel.index = 18; break; }
      unset_bit<0>(s); s += prime * 2 + 0;
      case 19: if (s >= sieve_end) { wheel.index = 19; break; }
      unset_bit<6>(s); s += prime * 4 + 2;
      case 20: if (s >= sieve_end) { wheel.index = 20; break; }
      unset_bit<1>(s); s += prime * 2 + 0;
      case 21: if (s >= sieve_end) { wheel.index = 21; break; }
      unset_bit<7>(s); s += prime * 4 + 2;
      case 22: if (s >= sieve_end) { wheel.index = 22; break; }
      unset_bit<3>(s); s += prime * 6 + 2;
      case 23: if (s >= sieve_end) { wheel.index = 23; break; }
      unset_bit<5>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 10 < sieve_end)
      {
        unset_bit<2>(s + prime *  0 +  0);
        unset_bit<4>(s + prime *  6 +  2);
        unset_bit<0>(s + prime * 10 +  4);
        unset_bit<6>(s + prime * 12 +  4);
        unset_bit<1>(s + prime * 16 +  6);
        unset_bit<7>(s + prime * 18 +  6);
        unset_bit<3>(s + prime * 22 +  8);
        unset_bit<5>(s + prime * 28 + 10);
        s += prime * 30 + 11;
      }
    }
    break;

    for (;;)
    {
      case 24: if (s >= sieve_end) { wheel.index = 24; break; }
      unset_bit<3>(s); s += prime * 6 + 3;
      case 25: if (s >= sieve_end) { wheel.index = 25; break; }
      unset_bit<0>(s); s += prime * 4 + 1;
      case 26: if (s >= sieve_end) { wheel.index = 26; break; }
      unset_bit<6>(s); s += prime * 2 + 1;
      case 27: if (s >= sieve_end) { wheel.index = 27; break; }
      unset_bit<5>(s); s += prime * 4 + 2;
      case 28: if (s >= sieve_end) { wheel.index = 28; break; }
      unset_bit<2>(s); s += prime * 2 + 1;
      case 29: if (s >= sieve_end) { wheel.index = 29; break; }
      unset_bit<1>(s); s += prime * 4 + 1;
      case 30: if (s >= sieve_end) { wheel.index = 30; break; }
      unset_bit<7>(s); s += prime * 6 + 3;
      case 31: if (s >= sieve_end) { wheel.index = 31; break; }
      unset_bit<4>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 12 < sieve_end)
      {
        unset_bit<3>(s + prime *  0 +  0);
        unset_bit<0>(s + prime *  6 +  3);
        unset_bit<6>(s + prime * 10 +  4);
        unset_bit<5>(s + prime * 12 +  5);
        unset_bit<2>(s + prime * 16 +  7);
        unset_bit<1>(s + prime * 18 +  8);
        unset_bit<7>(s + prime * 22 +  9);
        unset_bit<4>(s + prime * 28 + 12);
        s += prime * 30 + 13;
      }
    }
    break;

    for (;;)
    {
      case 32: if (s >= sieve_end) { wheel.index = 32; break; }
      unset_bit<4>(s); s += prime * 6 + 3;
      case 33: if (s >= sieve_end) { wheel.index = 33; break; }
      unset_bit<7>(s); s += prime * 4 + 3;
      case 34: if (s >= sieve_end) { wheel.index = 34; break; }
      unset_bit<1>(s); s += prime * 2 + 1;
      case 35: if (s >= sieve_end) { wheel.index = 35; break; }
      unset_bit<2>(s); s += prime * 4 + 2;
      case 36: if (s >= sieve_end) { wheel.index = 36; break; }
      unset_bit<5>(s); s += prime * 2 + 1;
      case 37: if (s >= sieve_end) { wheel.index = 37; break; }
      unset_bit<6>(s); s += prime * 4 + 3;
      case 38: if (s >= sieve_end) { wheel.index = 38; break; }
      unset_bit<0>(s); s += prime * 6 + 3;
      case 39: if (s >= sieve_end) { wheel.index = 39; break; }
      unset_bit<3>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 16 < sieve_end)
      {
        unset_bit<4>(s + prime *  0 +  0);
        unset_bit<7>(s + prime *  6 +  3);
        unset_bit<1>(s + prime * 10 +  6);
        unset_bit<2>(s + prime * 12 +  7);
        unset_bit<5>(s + prime * 16 +  9);
        unset_bit<6>(s + prime * 18 + 10);
        unset_bit<0>(s + prime * 22 + 13);
        unset_bit<3>(s + prime * 28 + 16);
        s += prime * 30 + 17;
      }
    }
    break;

    for (;;)
    {
      case 40: if (s >= sieve_end) { wheel.index = 40; break; }
      unset_bit<5>(s); s += prime * 6 + 4;
      case 41: if (s >= sieve_end) { wheel.index = 41; break; }
      unset_bit<3>(s); s += prime * 4 + 2;
      case 42: if (s >= sieve_end) { wheel.index = 42; break; }
      unset_bit<7>(s); s += prime * 2 + 2;
      case 43: if (s >= sieve_end) { wheel.index = 43; break; }
      unset_bit<1>(s); s += prime * 4 + 2;
      case 44: if (s >= sieve_end) { wheel.index = 44; break; }
      unset_bit<6>(s); s += prime * 2 + 2;
      case 45: if (s >= sieve_end) { wheel.index = 45; break; }
      unset_bit<0>(s); s += prime * 4 + 2;
      case 46: if (s >= sieve_end) { wheel.index = 46; break; }
      unset_bit<4>(s); s += prime * 6 + 4;
      case 47: if (s >= sieve_end) { wheel.index = 47; break; }
      unset_bit<2>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 18 < sieve_end)
      {
        unset_bit<5>(s + prime *  0 +  0);
        unset_bit<3>(s + prime *  6 +  4);
        unset_bit<7>(s + prime * 10 +  6);
        unset_bit<1>(s + prime * 12 +  8);
        unset_bit<6>(s + prime * 16 + 10);
        unset_bit<0>(s + prime * 18 + 12);
        unset_bit<4>(s + prime * 22 + 14);
        unset_bit<2>(s + prime * 28 + 18);
        s += prime * 30 + 19;
      }
    }
    break;

    for (;;)
    {
      case 48: if (s >= sieve_end) { wheel.index = 48; break; }
      unset_bit<6>(s); s += prime * 6 + 5;
      case 49: if (s >= sieve_end) { wheel.index = 49; break; }
      unset_bit<2>(s); s += prime * 4 + 3;
      case 50: if (s >= sieve_end) { wheel.index = 50; break; }
      unset_bit<3>(s); s += prime * 2 + 1;
      case 51: if (s >= sieve_end) { wheel.index = 51; break; }
      unset_bit<7>(s); s += prime * 4 + 4;
      case 52: if (s >= sieve_end) { wheel.index = 52; break; }
      unset_bit<0>(s); s += prime * 2 + 1;
      case 53: if (s >= sieve_end) { wheel.index = 53; break; }
      unset_bit<4>(s); s += prime * 4 + 3;
      case 54: if (s >= sieve_end) { wheel.index = 54; break; }
      unset_bit<5>(s); s += prime * 6 + 5;
      case 55: if (s >= sieve_end) { wheel.index = 55; break; }
      unset_bit<1>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 22 < sieve_end)
      {
        unset_bit<6>(s + prime *  0 +  0);
        unset_bit<2>(s + prime *  6 +  5);
        unset_bit<3>(s + prime * 10 +  8);
        unset_bit<7>(s + prime * 12 +  9);
        unset_bit<0>(s + prime * 16 + 13);
        unset_bit<4>(s + prime * 18 + 14);
        unset_bit<5>(s + prime * 22 + 17);
        unset_bit<1>(s + prime * 28 + 22);
        s += prime * 30 + 23;
      }
    }
    break;

    for (;;)
    {
      case 56: if (s >= sieve_end) { wheel.index = 56; break; }
      unset_bit<7>(s); s += prime * 6 + 6;
      case 57: if (s >= sieve_end) { wheel.index = 57; break; }
      unset_bit<6>(s); s += prime * 4 + 4;
      case 58: if (s >= sieve_end) { wheel.index = 58; break; }
      unset_bit<5>(s); s += prime * 2 + 2;
      case 59: if (s >= sieve_end) { wheel.index = 59; break; }
      unset_bit<4>(s); s += prime * 4 + 4;
      case 60: if (s >= sieve_end) { wheel.index = 60; break; }
      unset_bit<3>(s); s += prime * 2 + 2;
      case 61: if (s >= sieve_end) { wheel.index = 61; break; }
      unset_bit<2>(s); s += prime * 4 + 4;
      case 62: if (s >= sieve_end) { wheel.index = 62; break; }
      unset_bit<1>(s); s += prime * 6 + 6;
      case 63: if (s >= sieve_end) { wheel.index = 63; break; }
      unset_bit<0>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 28 < sieve_end)
      {
        unset_bit<7>(s + prime *  0 +  0);
        unset_bit<6>(s + prime *  6 +  6);
        unset_bit<5>(s + prime * 10 + 10);
        unset_bit<4>(s + prime * 12 + 12);
        unset_bit<3>(s + prime * 16 + 16);
        unset_bit<2>(s + prime * 18 + 18);
        unset_bit<1>(s + prime * 22 + 22);
        unset_bit<0>(s + prime * 28 + 28);
        s += prime * 30 + 29;
      }
    }
    break;
  }

  // update for the next segment
  wheel.multiple = (uint32_t) (s - sieve_end);
}

/// Remove the i-th prime and the multiples of the i-th
/// prime from the sieve array. Returns the count of
/// elements removed for the first time i.e. the count
/// of sieved elements whose least prime factor is the
/// i-th prime.
///
uint64_t Sieve::cross_off_count(uint64_t prime, uint64_t i)
{
  if (i >= wheel_.size())
    add(prime);

  uint64_t cnt = 0;
  Wheel& wheel = wheel_[i];
  // pointer to the byte with the 1st multiple
  byte_t* s = sieve_ + wheel.multiple;
  byte_t* sieve_end = sieve_ + sieve_size_;
  prime /= 30;

  switch (wheel.index)
  {
    for (;;)
    {
      case 0: if (s >= sieve_end) { wheel.index = 0; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 6 + 0;
      case 1: if (s >= sieve_end) { wheel.index = 1; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 4 + 0;
      case 2: if (s >= sieve_end) { wheel.index = 2; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 2 + 0;
      case 3: if (s >= sieve_end) { wheel.index = 3; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 4 + 0;
      case 4: if (s >= sieve_end) { wheel.index = 4; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 2 + 0;
      case 5: if (s >= sieve_end) { wheel.index = 5; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 4 + 0;
      case 6: if (s >= sieve_end) { wheel.index = 6; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 6 + 0;
      case 7: if (s >= sieve_end) { wheel.index = 7; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 2 + 1;

      while (s + prime * 28 < sieve_end)
      {
        cnt += is_bit<0>(s + prime *  0);
        unset_bit<0>(s + prime *  0);
        cnt += is_bit<1>(s + prime *  6);
        unset_bit<1>(s + prime *  6);
        cnt += is_bit<2>(s + prime * 10);
        unset_bit<2>(s + prime * 10);
        cnt += is_bit<3>(s + prime * 12);
        unset_bit<3>(s + prime * 12);
        cnt += is_bit<4>(s + prime * 16);
        unset_bit<4>(s + prime * 16);
        cnt += is_bit<5>(s + prime * 18);
        unset_bit<5>(s + prime * 18);
        cnt += is_bit<6>(s + prime * 22);
        unset_bit<6>(s + prime * 22);
        cnt += is_bit<7>(s + prime * 28);
        unset_bit<7>(s + prime * 28);
        s += prime * 30 + 1;
      }
    }
    break;

    for (;;)
    {
      case  8: if (s >= sieve_end) { wheel.index =  8; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 6 + 1;
      case  9: if (s >= sieve_end) { wheel.index =  9; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 4 + 1;
      case 10: if (s >= sieve_end) { wheel.index = 10; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 2 + 1;
      case 11: if (s >= sieve_end) { wheel.index = 11; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 4 + 0;
      case 12: if (s >= sieve_end) { wheel.index = 12; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 2 + 1;
      case 13: if (s >= sieve_end) { wheel.index = 13; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 4 + 1;
      case 14: if (s >= sieve_end) { wheel.index = 14; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 6 + 1;
      case 15: if (s >= sieve_end) { wheel.index = 15; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 6 < sieve_end)
      {
        cnt += is_bit<1>(s + prime *  0 + 0);
        unset_bit<1>(s + prime *  0 + 0);
        cnt += is_bit<5>(s + prime *  6 + 1);
        unset_bit<5>(s + prime *  6 + 1);
        cnt += is_bit<4>(s + prime * 10 + 2);
        unset_bit<4>(s + prime * 10 + 2);
        cnt += is_bit<0>(s + prime * 12 + 3);
        unset_bit<0>(s + prime * 12 + 3);
        cnt += is_bit<7>(s + prime * 16 + 3);
        unset_bit<7>(s + prime * 16 + 3);
        cnt += is_bit<3>(s + prime * 18 + 4);
        unset_bit<3>(s + prime * 18 + 4);
        cnt += is_bit<2>(s + prime * 22 + 5);
        unset_bit<2>(s + prime * 22 + 5);
        cnt += is_bit<6>(s + prime * 28 + 6);
        unset_bit<6>(s + prime * 28 + 6);
        s += prime * 30 + 7;
      }
    }
    break;

    for (;;)
    {
      case 16: if (s >= sieve_end) { wheel.index = 16; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 6 + 2;
      case 17: if (s >= sieve_end) { wheel.index = 17; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 4 + 2;
      case 18: if (s >= sieve_end) { wheel.index = 18; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 2 + 0;
      case 19: if (s >= sieve_end) { wheel.index = 19; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 4 + 2;
      case 20: if (s >= sieve_end) { wheel.index = 20; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 2 + 0;
      case 21: if (s >= sieve_end) { wheel.index = 21; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 4 + 2;
      case 22: if (s >= sieve_end) { wheel.index = 22; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 6 + 2;
      case 23: if (s >= sieve_end) { wheel.index = 23; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 10 < sieve_end)
      {
        cnt += is_bit<2>(s + prime *  0 +  0);
        unset_bit<2>(s + prime *  0 +  0);
        cnt += is_bit<4>(s + prime *  6 +  2);
        unset_bit<4>(s + prime *  6 +  2);
        cnt += is_bit<0>(s + prime * 10 +  4);
        unset_bit<0>(s + prime * 10 +  4);
        cnt += is_bit<6>(s + prime * 12 +  4);
        unset_bit<6>(s + prime * 12 +  4);
        cnt += is_bit<1>(s + prime * 16 +  6);
        unset_bit<1>(s + prime * 16 +  6);
        cnt += is_bit<7>(s + prime * 18 +  6);
        unset_bit<7>(s + prime * 18 +  6);
        cnt += is_bit<3>(s + prime * 22 +  8);
        unset_bit<3>(s + prime * 22 +  8);
        cnt += is_bit<5>(s + prime * 28 + 10);
        unset_bit<5>(s + prime * 28 + 10);
        s += prime * 30 + 11;
      }
    }
    break;

    for (;;)
    {
      case 24: if (s >= sieve_end) { wheel.index = 24; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 6 + 3;
      case 25: if (s >= sieve_end) { wheel.index = 25; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 4 + 1;
      case 26: if (s >= sieve_end) { wheel.index = 26; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 2 + 1;
      case 27: if (s >= sieve_end) { wheel.index = 27; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 4 + 2;
      case 28: if (s >= sieve_end) { wheel.index = 28; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 2 + 1;
      case 29: if (s >= sieve_end) { wheel.index = 29; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 4 + 1;
      case 30: if (s >= sieve_end) { wheel.index = 30; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 6 + 3;
      case 31: if (s >= sieve_end) { wheel.index = 31; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 12 < sieve_end)
      {
        cnt += is_bit<3>(s + prime *  0 +  0);
        unset_bit<3>(s + prime *  0 +  0);
        cnt += is_bit<0>(s + prime *  6 +  3);
        unset_bit<0>(s + prime *  6 +  3);
        cnt += is_bit<6>(s + prime * 10 +  4);
        unset_bit<6>(s + prime * 10 +  4);
        cnt += is_bit<5>(s + prime * 12 +  5);
        unset_bit<5>(s + prime * 12 +  5);
        cnt += is_bit<2>(s + prime * 16 +  7);
        unset_bit<2>(s + prime * 16 +  7);
        cnt += is_bit<1>(s + prime * 18 +  8);
        unset_bit<1>(s + prime * 18 +  8);
        cnt += is_bit<7>(s + prime * 22 +  9);
        unset_bit<7>(s + prime * 22 +  9);
        cnt += is_bit<4>(s + prime * 28 + 12);
        unset_bit<4>(s + prime * 28 + 12);
        s += prime * 30 + 13;
      }
    }
    break;

    for (;;)
    {
      case 32: if (s >= sieve_end) { wheel.index = 32; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 6 + 3;
      case 33: if (s >= sieve_end) { wheel.index = 33; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 4 + 3;
      case 34: if (s >= sieve_end) { wheel.index = 34; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 2 + 1;
      case 35: if (s >= sieve_end) { wheel.index = 35; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 4 + 2;
      case 36: if (s >= sieve_end) { wheel.index = 36; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 2 + 1;
      case 37: if (s >= sieve_end) { wheel.index = 37; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 4 + 3;
      case 38: if (s >= sieve_end) { wheel.index = 38; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 6 + 3;
      case 39: if (s >= sieve_end) { wheel.index = 39; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 16 < sieve_end)
      {
        cnt += is_bit<4>(s + prime *  0 +  0);
        unset_bit<4>(s + prime *  0 +  0);
        cnt += is_bit<7>(s + prime *  6 +  3);
        unset_bit<7>(s + prime *  6 +  3);
        cnt += is_bit<1>(s + prime * 10 +  6);
        unset_bit<1>(s + prime * 10 +  6);
        cnt += is_bit<2>(s + prime * 12 +  7);
        unset_bit<2>(s + prime * 12 +  7);
        cnt += is_bit<5>(s + prime * 16 +  9);
        unset_bit<5>(s + prime * 16 +  9);
        cnt += is_bit<6>(s + prime * 18 + 10);
        unset_bit<6>(s + prime * 18 + 10);
        cnt += is_bit<0>(s + prime * 22 + 13);
        unset_bit<0>(s + prime * 22 + 13);
        cnt += is_bit<3>(s + prime * 28 + 16);
        unset_bit<3>(s + prime * 28 + 16);
        s += prime * 30 + 17;
      }
    }
    break;

    for (;;)
    {
      case 40: if (s >= sieve_end) { wheel.index = 40; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 6 + 4;
      case 41: if (s >= sieve_end) { wheel.index = 41; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 4 + 2;
      case 42: if (s >= sieve_end) { wheel.index = 42; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 2 + 2;
      case 43: if (s >= sieve_end) { wheel.index = 43; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 4 + 2;
      case 44: if (s >= sieve_end) { wheel.index = 44; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 2 + 2;
      case 45: if (s >= sieve_end) { wheel.index = 45; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 4 + 2;
      case 46: if (s >= sieve_end) { wheel.index = 46; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 6 + 4;
      case 47: if (s >= sieve_end) { wheel.index = 47; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 18 < sieve_end)
      {
        cnt += is_bit<5>(s + prime *  0 +  0);
        unset_bit<5>(s + prime *  0 +  0);
        cnt += is_bit<3>(s + prime *  6 +  4);
        unset_bit<3>(s + prime *  6 +  4);
        cnt += is_bit<7>(s + prime * 10 +  6);
        unset_bit<7>(s + prime * 10 +  6);
        cnt += is_bit<1>(s + prime * 12 +  8);
        unset_bit<1>(s + prime * 12 +  8);
        cnt += is_bit<6>(s + prime * 16 + 10);
        unset_bit<6>(s + prime * 16 + 10);
        cnt += is_bit<0>(s + prime * 18 + 12);
        unset_bit<0>(s + prime * 18 + 12);
        cnt += is_bit<4>(s + prime * 22 + 14);
        unset_bit<4>(s + prime * 22 + 14);
        cnt += is_bit<2>(s + prime * 28 + 18);
        unset_bit<2>(s + prime * 28 + 18);
        s += prime * 30 + 19;
      }
    }
    break;

    for (;;)
    {
      case 48: if (s >= sieve_end) { wheel.index = 48; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 6 + 5;
      case 49: if (s >= sieve_end) { wheel.index = 49; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 4 + 3;
      case 50: if (s >= sieve_end) { wheel.index = 50; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 2 + 1;
      case 51: if (s >= sieve_end) { wheel.index = 51; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 4 + 4;
      case 52: if (s >= sieve_end) { wheel.index = 52; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 2 + 1;
      case 53: if (s >= sieve_end) { wheel.index = 53; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 4 + 3;
      case 54: if (s >= sieve_end) { wheel.index = 54; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 6 + 5;
      case 55: if (s >= sieve_end) { wheel.index = 55; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 22 < sieve_end)
      {
        cnt += is_bit<6>(s + prime *  0 +  0);
        unset_bit<6>(s + prime *  0 +  0);
        cnt += is_bit<2>(s + prime *  6 +  5);
        unset_bit<2>(s + prime *  6 +  5);
        cnt += is_bit<3>(s + prime * 10 +  8);
        unset_bit<3>(s + prime * 10 +  8);
        cnt += is_bit<7>(s + prime * 12 +  9);
        unset_bit<7>(s + prime * 12 +  9);
        cnt += is_bit<0>(s + prime * 16 + 13);
        unset_bit<0>(s + prime * 16 + 13);
        cnt += is_bit<4>(s + prime * 18 + 14);
        unset_bit<4>(s + prime * 18 + 14);
        cnt += is_bit<5>(s + prime * 22 + 17);
        unset_bit<5>(s + prime * 22 + 17);
        cnt += is_bit<1>(s + prime * 28 + 22);
        unset_bit<1>(s + prime * 28 + 22);
        s += prime * 30 + 23;
      }
    }
    break;

    for (;;)
    {
      case 56: if (s >= sieve_end) { wheel.index = 56; break; }
      cnt += is_bit<7>(s); unset_bit<7>(s); s += prime * 6 + 6;
      case 57: if (s >= sieve_end) { wheel.index = 57; break; }
      cnt += is_bit<6>(s); unset_bit<6>(s); s += prime * 4 + 4;
      case 58: if (s >= sieve_end) { wheel.index = 58; break; }
      cnt += is_bit<5>(s); unset_bit<5>(s); s += prime * 2 + 2;
      case 59: if (s >= sieve_end) { wheel.index = 59; break; }
      cnt += is_bit<4>(s); unset_bit<4>(s); s += prime * 4 + 4;
      case 60: if (s >= sieve_end) { wheel.index = 60; break; }
      cnt += is_bit<3>(s); unset_bit<3>(s); s += prime * 2 + 2;
      case 61: if (s >= sieve_end) { wheel.index = 61; break; }
      cnt += is_bit<2>(s); unset_bit<2>(s); s += prime * 4 + 4;
      case 62: if (s >= sieve_end) { wheel.index = 62; break; }
      cnt += is_bit<1>(s); unset_bit<1>(s); s += prime * 6 + 6;
      case 63: if (s >= sieve_end) { wheel.index = 63; break; }
      cnt += is_bit<0>(s); unset_bit<0>(s); s += prime * 2 + 1;

      while (s + prime * 28 + 28 < sieve_end)
      {
        cnt += is_bit<7>(s + prime *  0 +  0);
        unset_bit<7>(s + prime *  0 +  0);
        cnt += is_bit<6>(s + prime *  6 +  6);
        unset_bit<6>(s + prime *  6 +  6);
        cnt += is_bit<5>(s + prime * 10 + 10);
        unset_bit<5>(s + prime * 10 + 10);
        cnt += is_bit<4>(s + prime * 12 + 12);
        unset_bit<4>(s + prime * 12 + 12);
        cnt += is_bit<3>(s + prime * 16 + 16);
        unset_bit<3>(s + prime * 16 + 16);
        cnt += is_bit<2>(s + prime * 18 + 18);
        unset_bit<2>(s + prime * 18 + 18);
        cnt += is_bit<1>(s + prime * 22 + 22);
        unset_bit<1>(s + prime * 22 + 22);
        cnt += is_bit<0>(s + prime * 28 + 28);
        unset_bit<0>(s + prime * 28 + 28);
        s += prime * 30 + 29;
      }
    }
    break;
  }

  // update for the next segment
  wheel.multiple = (uint32_t) (s - sieve_end);

  return cnt;
}

} // namespace

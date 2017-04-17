///
/// @file  BitSieve.cpp
/// @brief The BitSieve class is a bit array for use with
///        Eratosthenes-like prime sieving algorithms. BitSieve
///        assigns 64 numbers to the bits of an 8 byte word thus
///        reducing the memory usage by a factor of 8.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BitSieve.hpp>
#include <popcnt.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

using namespace std;

namespace {

/// primes[1] = 2, primes[2] = 3, ...
const array<uint64_t, 10> primes = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

/// Bitmasks with multiples of the i-th prime set
const array<uint64_t, 10> masks =
{
  0x0000000000000000ull,
  0x5555555555555555ull, // 2
  0x9249249249249249ull, // 3
  0x1084210842108421ull, // 5
  0x8102040810204081ull, // 7
  0x0080100200400801ull, // 11
  0x0010008004002001ull, // 13
  0x0008000400020001ull, // 17
  0x0200004000080001ull, // 19
  0x0000400000800001ull  // 23
};

uint64_t unset_mask(uint64_t mask, uint64_t shift)
{
  return ~(mask << shift);
}

uint64_t fast_modulo(uint64_t x, uint64_t y)
{
  assert(x < y * 2);
  x = (x < y) ? x : x - y;
  return x;
}

} // namespace

namespace primecount {

const uint64_t BitSieve::unset_bit_[64] =
{
  ~(1ull <<  0), ~(1ull <<  1), ~(1ull <<  2), ~(1ull <<  3),
  ~(1ull <<  4), ~(1ull <<  5), ~(1ull <<  6), ~(1ull <<  7),
  ~(1ull <<  8), ~(1ull <<  9), ~(1ull << 10), ~(1ull << 11),
  ~(1ull << 12), ~(1ull << 13), ~(1ull << 14), ~(1ull << 15),
  ~(1ull << 16), ~(1ull << 17), ~(1ull << 18), ~(1ull << 19),
  ~(1ull << 20), ~(1ull << 21), ~(1ull << 22), ~(1ull << 23),
  ~(1ull << 24), ~(1ull << 25), ~(1ull << 26), ~(1ull << 27),
  ~(1ull << 28), ~(1ull << 29), ~(1ull << 30), ~(1ull << 31),
  ~(1ull << 32), ~(1ull << 33), ~(1ull << 34), ~(1ull << 35),
  ~(1ull << 36), ~(1ull << 37), ~(1ull << 38), ~(1ull << 39),
  ~(1ull << 40), ~(1ull << 41), ~(1ull << 42), ~(1ull << 43),
  ~(1ull << 44), ~(1ull << 45), ~(1ull << 46), ~(1ull << 47),
  ~(1ull << 48), ~(1ull << 49), ~(1ull << 50), ~(1ull << 51),
  ~(1ull << 52), ~(1ull << 53), ~(1ull << 54), ~(1ull << 55),
  ~(1ull << 56), ~(1ull << 57), ~(1ull << 58), ~(1ull << 59),
  ~(1ull << 60), ~(1ull << 61), ~(1ull << 62), ~(1ull << 63)
};

BitSieve::BitSieve(std::size_t size) :
  sieve_(ceil_div(size, 64)),
  size_(size)
{ }

/// Pre-sieve the multiples (>= low) of the first c primes.
/// @warning Also removes the first c primes.
/// @pre c < 10
///
void BitSieve::pre_sieve(uint64_t c, uint64_t low)
{
  assert(c < primes.size());

  if (sieve_.empty())
    return;

  // unset multiples of 2
  sieve_[0] = unset_mask(masks[1], low % 2);

  uint64_t last = 1;
  uint64_t sieved = last;
  uint64_t sieve_size = sieve_.size();

  // pre-sieve multiples of first c primes
  for (uint64_t i = 2; i <= c; i++)
  {
    uint64_t prime = primes[i];
    uint64_t end_copy = min(sieved * prime, sieve_size);

    // pre-sieve multiples of primes < i-th prime
    // by copying a small pre-sieved buffer
    while (last < end_copy)
    {
      uint64_t copy_words = min(sieved, sieve_size - last);
      copy(sieve_.begin(),
           sieve_.begin() + copy_words,
           sieve_.begin() + last);
      last += copy_words;
    }

    // calculate first multiple >= low of prime
    uint64_t multiple = ceil_div(low, prime) * prime;
    uint64_t shift = multiple - low;
    uint64_t next_shift = prime - 64 % prime;
    uint64_t mask = masks[i];

    // pre-sieve multiples of the i-th prime
    for (uint64_t j = 0; j < last; j++)
    {
      sieve_[j] &= unset_mask(mask, shift);
      shift = fast_modulo(shift + next_shift, prime);
    }

    sieved = last;
  }

  // fill up the rest of the sieve
  while (last < sieve_size)
  {
    uint64_t copy_words = min(sieved, sieve_size - last);
    copy(sieve_.begin(),
         sieve_.begin() + copy_words,
         sieve_.begin() + last);
    last += copy_words;
  }
}

/// Count the number of 1 bits inside [start, stop]
uint64_t BitSieve::count(uint64_t start,
                         uint64_t stop) const
{
  if (start > stop)
    return 0;

  assert(stop < size_);

  uint64_t start_idx = start / 64;
  uint64_t stop_idx = stop / 64;
  uint64_t m1 = 0xffffffffffffffffull << (start % 64);
  uint64_t m2 = 0xffffffffffffffffull >> (63 - stop % 64);
  uint64_t bit_count;

  if (start_idx == stop_idx)
    bit_count = popcnt64(sieve_[start_idx] & (m1 & m2));
  else
  {
    bit_count = popcnt64(sieve_[start_idx] & m1);
    bit_count += popcnt(&sieve_[start_idx + 1], stop_idx - (start_idx + 1));
    bit_count += popcnt64(sieve_[stop_idx] & m2);
  }

  return bit_count;
}

} // namespace

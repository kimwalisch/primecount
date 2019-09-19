///
/// @file  Sieve.hpp
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

#ifndef SIEVE_HPP
#define SIEVE_HPP

#include <stdint.h>
#include <cassert>
#include <memory>
#include <vector>

namespace primecount {

using byte_t = uint8_t;

struct Wheel
{
  Wheel()
    : multiple(0),
      index(0)
  { }
  Wheel(uint32_t m, uint32_t i)
    : multiple(m),
      index(i)
  { }
  uint32_t multiple;
  uint32_t index;
};

class Sieve
{
public:
  Sieve(uint64_t start, uint64_t segment_size, uint64_t wheel_size);
  static uint64_t get_segment_size(uint64_t size);
  uint64_t segment_size() const;
  void cross_off(uint64_t prime, uint64_t i);
  uint64_t cross_off_count(uint64_t prime, uint64_t i);

  template <typename T>
  void pre_sieve(const std::vector<T>& primes, uint64_t c, uint64_t low, uint64_t high)
  {
    assert(c < primes.size());
    reset_sieve(low, high);

    for (uint64_t i = 4; i <= c; i++)
      cross_off(primes[i], i);
  }

  /// Count 1 bits inside [start, stop]
  uint64_t count(uint64_t start, uint64_t stop) const;

  /// Count 1 bits inside [0, stop]
  uint64_t count(uint64_t stop) const
  {
    return count(0, stop);
  }

  /// Count 1 bits inside [start, stop].
  /// This method counts either forwards or backwards 
  /// depending on what's faster.
  ///
  uint64_t count(uint64_t start,
                 uint64_t stop,
                 uint64_t low,
                 uint64_t high,
                 uint64_t count_0_start,
                 uint64_t count_low_high) const
  {
    if (start > stop)
      return 0;

    if (stop - start < ((high - 1) - low) - stop)
      return count(start, stop);
    else
      return count_low_high - count_0_start - count(stop + 1, (high - 1) - low);
  }

private:
  void reset_sieve(uint64_t low, uint64_t high);
  void add(uint64_t prime);

  uint64_t start_;
  uint64_t sieve_size_;
  byte_t* sieve_;
  std::unique_ptr<byte_t[]> deleter_;
  std::vector<Wheel> wheel_;
};

} // namespace

#endif

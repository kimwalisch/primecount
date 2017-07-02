///
/// @file  Sieve.hpp
/// @brief The Sieve class is a highly optimized sieve of
///        Eratosthenes implementation with 30 numbers per
///        byte i.e. the 8 bits of each byte correspond to
///        the offsets { 1, 7, 11, 13, 17, 19, 23, 29 }.
///        The Sieve also skips multiples of 2, 3, 5 using
///        wheel factorization.
///
///        Unlike a traditional prime sieve this sieve is
///        designed for use in the combinatorial prime
///        counting algorithms: This sieve removes primes
///        as well as multiples of primes and it counts
///        the number of elements that have been crossed
///        off for the first time in the sieve array.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_HPP
#define SIEVE_HPP

#include <stdint.h>
#include <vector>

namespace primecount {

typedef uint8_t byte_t;

struct Wheel
{
  Wheel() { }
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
  Sieve(uint64_t start, 
      uint64_t stop, 
      uint64_t segment_size, 
      uint64_t wheel_size);

  void reset();
  uint64_t distance() const;
  static uint64_t get_segment_size(uint64_t size);
  uint64_t cross_off(uint64_t i, uint64_t prime);

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
  void set_sieve_size(uint64_t segment_size);
  void add_wheel(uint64_t prime);

  uint64_t low_;
  uint64_t high_;
  uint64_t start_;
  uint64_t stop_;
  bool first_segment_;

  std::vector<byte_t> sieve_;
  std::vector<Wheel> wheel_;
};

} // namespace

#endif

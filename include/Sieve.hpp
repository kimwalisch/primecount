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
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_HPP
#define SIEVE_HPP

#include <noinline.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <cassert>
#include <vector>

namespace primecount {

class Sieve
{
public:
  NOINLINE Sieve(uint64_t low, uint64_t segment_size, uint64_t wheel_size);
  NOINLINE void cross_off(uint64_t prime, uint64_t i);
  NOINLINE void cross_off_count(uint64_t prime, uint64_t i);
  static uint64_t get_segment_size(uint64_t size);
  uint64_t count(uint64_t start, uint64_t stop) const;
  uint64_t count(uint64_t stop);

  uint64_t get_total_count() const
  {
    return total_count_;
  }

  template <typename T>
  void pre_sieve(const std::vector<T>& primes, uint64_t c, uint64_t low, uint64_t high)
  {
    assert(c < primes.size());
    reset_sieve(low, high);

    for (uint64_t i = 4; i <= c; i++)
      cross_off(primes[i], i);

    init_counters(low, high);
  }

private:
  void add(uint64_t prime);
  void allocate_counters(uint64_t low);
  void reset_counters();
  void reset_sieve(uint64_t low, uint64_t high);
  void init_counters(uint64_t low, uint64_t high);
  uint64_t segment_size() const;

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

  uint64_t start_ = 0;
  uint64_t prev_stop_ = 0;
  uint64_t count_ = 0;
  uint64_t total_count_ = 0;
  uint64_t counters_i_ = 0;
  uint64_t counters_count_ = 0;
  uint64_t counters_dist_ = 0;
  uint64_t counters_dist_log2_ = 0;
  uint64_t counters_stop_ = 0;
  pod_vector<uint8_t> sieve_;
  pod_vector<uint64_t> counters_;
  std::vector<Wheel> wheel_;
};

} // namespace

#endif

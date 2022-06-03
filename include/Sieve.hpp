///
/// @file  Sieve.hpp
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
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.md
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_HPP
#define SIEVE_HPP

#include <pod_vector.hpp>
#include <stdint.h>

namespace primecount {

class Sieve
{
public:
  Sieve(uint64_t low, uint64_t segment_size, uint64_t wheel_size);
  void cross_off(uint64_t prime, uint64_t i);
  void cross_off_count(uint64_t prime, uint64_t i);
  static uint64_t get_segment_size(uint64_t size);
  uint64_t count(uint64_t start, uint64_t stop) const;
  uint64_t count(uint64_t stop);

  uint64_t get_total_count() const
  {
    return total_count_;
  }

  template <typename T>
  void pre_sieve(const pod_vector<T>& primes, uint64_t c, uint64_t low, uint64_t high)
  {
    reset_sieve(low, high);

    for (uint64_t i = 4; i <= c; i++)
      cross_off(primes[i], i);

    init_counter(low, high);
  }

private:
  void add(uint64_t prime);
  void allocate_counter(uint64_t low);
  void init_counter(uint64_t low, uint64_t high);
  void reset_counter();
  void reset_sieve(uint64_t low, uint64_t high);
  uint64_t segment_size() const;

  struct Wheel
  {
    uint32_t multiple;
    uint32_t index;

    Wheel() = default;
    Wheel(uint32_t m, uint32_t i)
      : multiple(m),
        index(i)
    { }
  };

  struct Counter
  {
    uint64_t stop = 0;
    uint64_t dist = 0;
    uint64_t log2_dist = 0;
    uint64_t sum = 0;
    uint64_t i = 0;
    pod_vector<uint32_t> counter;

    uint32_t& operator[](std::size_t pos)
    {
      return counter[pos];
    }
  
    uint32_t operator[](std::size_t pos) const
    {
      return counter[pos];
    }
  };

  uint64_t start_ = 0;
  uint64_t prev_stop_ = 0;
  uint64_t count_ = 0;
  uint64_t total_count_ = 0;
  pod_vector<uint8_t> sieve_;
  pod_vector<Wheel> wheel_;
  Counter counter_;
};

} // namespace

#endif

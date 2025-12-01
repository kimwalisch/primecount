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
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_HPP
#define SIEVE_HPP

#include <cpu_arch_macros.hpp>
#include <print.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primecount {

class Sieve
{
public:
  Sieve(uint64_t low, uint64_t segment_size, uint64_t wheel_size);
  uint64_t count(uint64_t stop);
  uint64_t count(uint64_t start, uint64_t stop) const;
  void init_counter(uint64_t low, uint64_t high);
  void cross_off(uint64_t prime, uint64_t i);
  void cross_off_count(uint64_t prime, uint64_t i);
  static uint64_t align_segment_size(uint64_t size);
  static string_view_t count_algo_name();

  uint64_t get_total_count() const
  {
    return total_count_;
  }

  template <typename T>
  void pre_sieve(const Vector<T>& primes, uint64_t c, uint64_t low, uint64_t high)
  {
    uint64_t primePi = pre_sieve(c, low);
    resize_sieve(low, high);

    for (uint64_t i = primePi + 1; i <= c; i++)
      cross_off(primes[i], i);
  }

#if defined(ENABLE_PORTABLE_POPCNT64)
  /// Count 1 bits inside [0, stop]
  uint64_t count_popcnt64(uint64_t stop);

  /// Count 1 bits inside [start, stop]
  uint64_t count_popcnt64(uint64_t start, uint64_t stop) const;
#endif

#if defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

  /// Count 1 bits inside [0, stop]
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    __attribute__ ((target ("avx512f,avx512vpopcntdq")))
  #endif
  uint64_t count_avx512(uint64_t stop);

  /// Count 1 bits inside [start, stop]
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    __attribute__ ((target ("avx512f,avx512vpopcntdq")))
  #endif
  uint64_t count_avx512(uint64_t start, uint64_t stop) const;

#elif defined(ENABLE_ARM_SVE) || \
      defined(ENABLE_MULTIARCH_ARM_SVE)

  /// Count 1 bits inside [0, stop]
  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    __attribute__ ((target ("arch=armv8-a+sve")))
  #endif
  uint64_t count_arm_sve(uint64_t stop);

  /// Count 1 bits inside [start, stop]
  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    __attribute__ ((target ("arch=armv8-a+sve")))
  #endif
  uint64_t count_arm_sve(uint64_t start, uint64_t stop) const;

#endif

private:
  void add(uint64_t prime, uint64_t i);
  void allocate_counter(uint64_t low);
  void reset_counter();
  void resize_sieve(uint64_t low, uint64_t high);
  uint64_t pre_sieve(uint64_t c, uint64_t low);
  uint64_t segment_size() const;
  static const Array<uint64_t, 240> unset_smaller;
  static const Array<uint64_t, 240> unset_larger;

  struct Wheel
  {
    uint32_t multiple;
    uint32_t index;
  };

  struct Counter
  {
    uint64_t stop = 0;
    uint64_t dist = 0;
    uint64_t log2_dist = 0;
    uint64_t sum = 0;
    uint64_t i = 0;
    Vector<uint32_t> counter;

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
  Vector<uint8_t> sieve_;
  Vector<Wheel> wheel_;
  Counter counter_;
};

} // namespace

// These forced inline Sieve::count(stop) methods
// must be defined in a header file that is
// visible to all translation units that use them.
#include "Sieve_count_stop.hpp"

#endif

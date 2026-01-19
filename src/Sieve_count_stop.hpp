///
/// @file  Sieve_count_stop.hpp
/// @brief Highly optimized methods to count the number of 1 bits
///        in the sieve array in the special leaves algorithm
///        (used in the combinatorial prime counting algorithms
///        e.g. Lagarias-Miller-Odlyzko, Deleglise-Rivat, Gourdon).
///
///        The Sieve::count(stop) methods defined in this file are
///        called very frequently. Therefore, these methods have all
///        been annotated using ALWAYS_INLINE.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_COUNT_STOP_HPP
#define SIEVE_COUNT_STOP_HPP

#include <Sieve.hpp>
#include <Sieve_count_simd.hpp>
#include <cpu_arch_macros.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

namespace primecount {

/// Count 1 bits inside [0, stop].
/// This method is safe to run on any CPU without runtime
/// CPUID checks. In most cases (e.g. when compiled
/// without -march=native) this will call the portable
/// but slow count_popcnt64() method.
///
ALWAYS_INLINE uint64_t Sieve::count(uint64_t stop)
{
  #if defined(ENABLE_AVX512_VPOPCNT)
    #define DEFAULT_SIEVE_COUNT_ALGO_NAME "AVX512 bit counting"
    return count_avx512(stop);
  #elif defined(ENABLE_ARM_SVE)
    #define DEFAULT_SIEVE_COUNT_ALGO_NAME "ARM SVE bit counting"
    return count_arm_sve(stop);
  #else
    #define DEFAULT_SIEVE_COUNT_ALGO_NAME "POPCNT64 bit counting"
    return count_popcnt64(stop);
  #endif
}

#if defined(ENABLE_PORTABLE_POPCNT64)

/// Count 1 bits inside [0, stop]
ALWAYS_INLINE uint64_t Sieve::count_popcnt64(uint64_t stop)
{
  ASSERT(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

  if (start > stop)
    return count_;

  // Quickly count the number of unsieved elements (in
  // the sieve array) up to a value that is close to
  // the stop number i.e. (stop - start) < counter_.dist.
  // We do this using the counter array, each element
  // of the counter array contains the number of
  // unsieved elements in the interval:
  // [i * counter_.dist, (i + 1) * counter_.dist[.
  while (counter_.stop <= stop)
  {
    start = counter_.stop;
    counter_.stop += counter_.dist;
    counter_.sum += counter_[counter_.i++];
    count_ = counter_.sum;
  }

  // Here the remaining distance is relatively small i.e.
  // (stop - start) < counter_.dist, hence we simply
  // count the remaining number of unsieved elements by
  // linearly iterating over the sieve array.
  SIEVE_COUNT_POPCNT64(start, stop);
  count_ += cnt;

  return count_;
}

#endif

#if defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

/// Count 1 bits inside [0, stop]
#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  __attribute__ ((target ("avx512f,avx512vpopcntdq")))
#endif
ALWAYS_INLINE uint64_t Sieve::count_avx512(uint64_t stop)
{
  ASSERT(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

  if (start > stop)
    return count_;

  // Quickly count the number of unsieved elements (in
  // the sieve array) up to a value that is close to
  // the stop number i.e. (stop - start) < counter_.dist.
  // We do this using the counter array, each element
  // of the counter array contains the number of
  // unsieved elements in the interval:
  // [i * counter_.dist, (i + 1) * counter_.dist[.
  while (counter_.stop <= stop)
  {
    start = counter_.stop;
    counter_.stop += counter_.dist;
    counter_.sum += counter_[counter_.i++];
    count_ = counter_.sum;
  }

  // Here the remaining distance is relatively small i.e.
  // (stop - start) < counter_.dist, hence we simply
  // count the remaining number of unsieved elements by
  // linearly iterating over the sieve array.
  SIEVE_COUNT_AVX512(start, stop);
  count_ += cnt;

  return count_;
}

#elif defined(ENABLE_ARM_SVE) || \
      defined(ENABLE_MULTIARCH_ARM_SVE)

/// Count 1 bits inside [0, stop]
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("arch=armv8-a+sve")))
#endif
ALWAYS_INLINE uint64_t Sieve::count_arm_sve(uint64_t stop)
{
  ASSERT(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

  if (start > stop)
    return count_;

  // Quickly count the number of unsieved elements (in
  // the sieve array) up to a value that is close to
  // the stop number i.e. (stop - start) < counter_.dist.
  // We do this using the counter array, each element
  // of the counter array contains the number of
  // unsieved elements in the interval:
  // [i * counter_.dist, (i + 1) * counter_.dist[.
  while (counter_.stop <= stop)
  {
    start = counter_.stop;
    counter_.stop += counter_.dist;
    counter_.sum += counter_[counter_.i++];
    count_ = counter_.sum;
  }

  // Here the remaining distance is relatively small i.e.
  // (stop - start) < counter_.dist, hence we simply
  // count the remaining number of unsieved elements by
  // linearly iterating over the sieve array.
  SIEVE_COUNT_ARM_SVE(start, stop);
  count_ += cnt;

  return count_;
}

#endif

} // namespace

#endif

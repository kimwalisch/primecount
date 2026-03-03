///
/// @file  Sieve_count_start_stop.hpp
/// @brief Highly optimized methods to count the number of 1 bits
///        in the sieve array in the special leaves algorithm
///        (used in the combinatorial prime counting algorithms
///        e.g. Lagarias-Miller-Odlyzko, Deleglise-Rivat, Gourdon).
///
///        The methods defined in this file are only called from
///        Sieve.cpp. Hence, "Sieve_count_start_stop.hpp" is only
///        included in Sieve.cpp. The Sieve::count(start, stop)
///        methods defined in this file are called much less
///        frequently than the Sieve::count(stop) methods. Hence,
///        the Sieve::count(start, stop) methods have not been
///        annotated using ALWAYS_INLINE.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_COUNT_START_STOP_HPP
#define SIEVE_COUNT_START_STOP_HPP

#include <Sieve.hpp>
#include <Sieve_count_simd.hpp>
#include <cpu_arch_macros.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <cpu_supports_avx512_vpopcnt.hpp>
#endif

namespace {

#if defined(ENABLE_MULTIARCH_ARM_SVE)

/// svcntd() requires -march=armv8-a+sve compiler flag
__attribute__ ((target ("arch=armv8-a+sve")))
uint64_t get_svcntd()
{
  return svcntd();
}

#endif

uint64_t bytes_per_count_instruction()
{
  #if defined(ENABLE_AVX512_VPOPCNT)
    // count_avx512() algorithm
    return sizeof(__m512i);
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    // count_avx512() algorithm
    if (cpu_supports_avx512_vpopcnt)
      return sizeof(__m512i);
  #elif defined(ENABLE_ARM_SVE)
    // count_arm_sve() algorithm
    return svcntd() * sizeof(uint64_t);
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    // count_arm_sve() algorithm
    if (cpu_supports_sve)
      return get_svcntd() * sizeof(uint64_t);
  #endif

  // Default count_popcnt64() algorithm
  return sizeof(uint64_t);
}

} // namespace

namespace primecount {

string_view_t Sieve::count_algo_name()
{
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    if (cpu_supports_avx512_vpopcnt)
      return "Algorithm: AVX512 bit counting";
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      return "Algorithm: ARM SVE bit counting";
  #endif

  return "Algorithm: " DEFAULT_SIEVE_COUNT_ALGO_NAME;
}

/// Count 1 bits inside [start, stop]
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  #if defined(ENABLE_ARM_SVE)
    return count_arm_sve(start, stop);
  #elif defined(ENABLE_AVX512_VPOPCNT)
    return count_avx512(start, stop);
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    return cpu_supports_sve ? count_arm_sve(start, stop) : count_popcnt64(start, stop);
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    return cpu_supports_avx512_vpopcnt ? count_avx512(start, stop) : count_popcnt64(start, stop);
  #else
    return count_popcnt64(start, stop);
  #endif
}

#if defined(ENABLE_PORTABLE_POPCNT64)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
uint64_t Sieve::count_popcnt64(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  SIEVE_COUNT_POPCNT64(start, stop);
  return cnt;
}

#endif

#if defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  __attribute__ ((target ("avx512f,avx512vpopcntdq")))
#endif
uint64_t Sieve::count_avx512(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  SIEVE_COUNT_AVX512(start, stop);
  return cnt;
}

#elif defined(ENABLE_ARM_SVE) || \
      defined(ENABLE_MULTIARCH_ARM_SVE)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
#if defined(ENABLE_MULTIARCH_ARM_SVE)
  __attribute__ ((target ("arch=armv8-a+sve")))
#endif
uint64_t Sieve::count_arm_sve(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  SIEVE_COUNT_ARM_SVE(start, stop);
  return cnt;
}

#endif

} // namespace

#endif

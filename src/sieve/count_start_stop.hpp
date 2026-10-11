///
/// @file  count_start_stop.hpp
/// @brief Highly optimized methods to count the number of 1 bits
///        in the sieve array in the special leaves algorithm
///        (used in the combinatorial prime counting algorithms
///        e.g. Lagarias-Miller-Odlyzko, Deleglise-Rivat, Gourdon).
///
///        The methods defined in this file are only called from
///        Sieve.cpp. Hence, "count_start_stop.hpp" is only
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

#ifndef COUNT_START_STOP_HPP
#define COUNT_START_STOP_HPP

#include "Sieve.hpp"
#include "count_simd.hpp"

#include <cpu_arch_macros.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <cpu_supports_avx512_vpopcnt.hpp>
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
#endif

namespace {

#if defined(ENABLE_MULTIARCH_ARM_SVE)

/// svcntd() requires ARM SVE compiler flag
PRIMECOUNT_TARGET_ENTRY("+sve")
uint64_t get_svcntd()
{
  return svcntd();
}

#endif

ALWAYS_INLINE uint64_t bytes_per_count_instruction()
{
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    // count_avx512() algorithm
    if (cpu_supports_avx512_vpopcnt)
      return sizeof(__m512i);
  #endif

  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    // count_arm_sve() algorithm
    if (cpu_supports_sve)
      return get_svcntd() * sizeof(uint64_t);
  #endif

  return DEFAULT_BYTES_PER_COUNT_INSTRUCTION;
}

} // namespace

namespace primecount {

/// Count 1 bits inside [start, stop]
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    if (cpu_supports_avx512_vpopcnt)
      return count_avx512(start, stop);
  #endif

  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      return count_arm_sve(start, stop);
  #endif

  return SIEVE_COUNT_DEFAULT(start, stop);
}

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
uint64_t Sieve::SIEVE_COUNT_DEFAULT(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  SIEVE_COUNT_KERNEL_DEFAULT(start, stop);
  return cnt;
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
PRIMECOUNT_TARGET_ENTRY("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
uint64_t Sieve::count_avx512(uint64_t start, uint64_t stop) const
{
  if (start > stop)
    return 0;

  SIEVE_COUNT_AVX512(start, stop);
  return cnt;
}

#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)

/// Count 1 bits inside [start, stop].
/// The distance [start, stop] is small here < sqrt(segment_size),
/// hence we simply count the number of unsieved elements
/// by linearly iterating over the sieve array.
///
PRIMECOUNT_TARGET_ENTRY("+sve")
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

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
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_COUNT_START_STOP_HPP
#define SIEVE_COUNT_START_STOP_HPP

#include <Sieve.hpp>
#include <cpu_arch_macros.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

#if defined(ENABLE_ARM_SVE)
  #include <arm_sve.h>
#elif defined(ENABLE_AVX512_VPOPCNT)
  #include <immintrin.h>
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <arm_sve.h>
  #include <cpu_supports_arm_sve.hpp>
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <immintrin.h>
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

  ASSERT(stop - start < segment_size());
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  // Branchfree bitmask calculation:
  // if (start_idx == stop_idx) m1 = m1 & m2;
  // if (start_idx == stop_idx) m2 = 0;
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2);
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0);

  const uint64_t* sieve64 = (const uint64_t*) sieve_.data();
  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  uint64_t cnt = popcnt64(start_bits);
  cnt += popcnt64(stop_bits);

  for (uint64_t i = start_idx + 1; i < stop_idx; i++)
    cnt += popcnt64(sieve64[i]);

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

  ASSERT(stop - start < segment_size());
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  // Branchfree bitmask calculation:
  // if (start_idx == stop_idx) m1 = m1 & m2;
  // if (start_idx == stop_idx) m2 = 0;
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2);
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0);

  const uint64_t* sieve64 = (const uint64_t*) sieve_.data();
  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  __m512i vec = _mm512_set_epi64(0, 0, 0, 0, 0, 0, stop_bits, start_bits);
  __m512i vcnt = _mm512_popcnt_epi64(vec);
  uint64_t i = start_idx + 1;

  // Compute this for loop using AVX512.
  // for (i = start_idx + 1; i < stop_idx; i++)
  //   cnt += popcnt64(sieve64[i]);

  for (; i + 8 < stop_idx; i += 8)
  {
    vec = _mm512_loadu_epi64(&sieve64[i]);
    vec = _mm512_popcnt_epi64(vec);
    vcnt = _mm512_add_epi64(vcnt, vec);
  }

  __mmask8 mask = (__mmask8) (0xff >> (i + 8 - stop_idx));
  vec = _mm512_maskz_loadu_epi64(mask, &sieve64[i]);
  vec = _mm512_popcnt_epi64(vec);
  vcnt = _mm512_add_epi64(vcnt, vec);
  return _mm512_reduce_add_epi64(vcnt);
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

  ASSERT(stop - start < segment_size());
  uint64_t start_idx = start / 240;
  uint64_t stop_idx = stop / 240;
  uint64_t m1 = unset_smaller[start % 240];
  uint64_t m2 = unset_larger[stop % 240];

  // Branchfree bitmask calculation:
  // if (start_idx == stop_idx) m1 = m1 & m2;
  // if (start_idx == stop_idx) m2 = 0;
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2);
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0);

  const uint64_t* sieve64 = (const uint64_t*) sieve_.data();
  uint64_t start_bits = sieve64[start_idx] & m1;
  uint64_t stop_bits = sieve64[stop_idx] & m2;
  ASSERT(svcntd() >= 2);
  svuint64_t vec = svinsr_n_u64(svdup_u64(start_bits), stop_bits);
  svuint64_t vcnt = svcnt_u64_z(svwhilelt_b64(0, 2), vec);
  uint64_t i = start_idx + 1;

  // Compute this for loop using ARM SVE.
  // for (i = start_idx + 1; i < stop_idx; i++)
  //   cnt += popcnt64(sieve64[i]);

  for (; i + svcntd() < stop_idx; i += svcntd())
  {
    vec = svld1_u64(svptrue_b64(), &sieve64[i]);
    vec = svcnt_u64_x(svptrue_b64(), vec);
    vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec);
  }

  svbool_t pg = svwhilelt_b64(i, stop_idx);
  vec = svld1_u64(pg, &sieve64[i]);
  vec = svcnt_u64_z(pg, vec);
  vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec);
  return svaddv_u64(svptrue_b64(), vcnt);
}

#endif

} // namespace

#endif

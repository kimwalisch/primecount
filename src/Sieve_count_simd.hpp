///
/// @file  Sieve_count_simd.hpp
/// @brief Highly optimized code to count the number of 1 bits in
///        the sieve array using SIMD instructions.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVE_COUNT_SIMD_HPP
#define SIEVE_COUNT_SIMD_HPP

#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

#if defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <arm_sve.h>
#elif defined(ENABLE_AVX512_VPOPCNT) || \
      defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <immintrin.h>
#endif

/// POPCNT64 /////////////////////////////////////////////////////////

/// Count 1 bits inside [start, stop] using POPCNT64
#define SIEVE_COUNT_POPCNT64(start, stop) \
  ASSERT(start <= stop); \
  ASSERT(stop - start < segment_size()); \
  uint64_t start_idx = start / 240; \
  uint64_t stop_idx = stop / 240; \
  uint64_t m1 = unset_smaller[start % 240]; \
  uint64_t m2 = unset_larger[stop % 240]; \
  \
  /* Branchfree bitmask calculation: */ \
  /* if (start_idx == stop_idx) m1 = m1 & m2; */ \
  /* if (start_idx == stop_idx) m2 = 0; */ \
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2); \
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0); \
  \
  const uint64_t* sieve64 = (const uint64_t*) sieve_.data(); \
  uint64_t start_bits = sieve64[start_idx] & m1; \
  uint64_t stop_bits = sieve64[stop_idx] & m2; \
  uint64_t cnt = popcnt64(start_bits); \
  cnt += popcnt64(stop_bits); \
  \
  for (uint64_t i = start_idx + 1; i < stop_idx; i++) \
    cnt += popcnt64(sieve64[i]);

/// AVX512 ///////////////////////////////////////////////////////////

/// Count 1 bits inside [start, stop] using AVX512
#define SIEVE_COUNT_AVX512(start, stop) \
  ASSERT(start <= stop); \
  ASSERT(stop - start < segment_size()); \
  uint64_t start_idx = start / 240; \
  uint64_t stop_idx = stop / 240; \
  uint64_t m1 = unset_smaller[start % 240]; \
  uint64_t m2 = unset_larger[stop % 240]; \
  \
  /* Branchfree bitmask calculation: */ \
  /* if (start_idx == stop_idx) m1 = m1 & m2; */ \
  /* if (start_idx == stop_idx) m2 = 0; */ \
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2); \
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0); \
  \
  const uint64_t* sieve64 = (const uint64_t*) sieve_.data(); \
  uint64_t start_bits = sieve64[start_idx] & m1; \
  uint64_t stop_bits = sieve64[stop_idx] & m2; \
  __m512i vec = _mm512_set_epi64(0, 0, 0, 0, 0, 0, stop_bits, start_bits); \
  __m512i vcnt = _mm512_popcnt_epi64(vec); \
  uint64_t i = start_idx + 1; \
  \
  /* Compute this for loop using AVX512. */ \
  /* for (i = start_idx + 1; i < stop_idx; i++) */ \
  /*   cnt += popcnt64(sieve64[i]); */ \
  for (; i + 8 < stop_idx; i += 8) \
  { \
    vec = _mm512_loadu_epi64(&sieve64[i]); \
    vec = _mm512_popcnt_epi64(vec); \
    vcnt = _mm512_add_epi64(vcnt, vec); \
  } \
  __mmask8 mask = (__mmask8) (0xff >> (i + 8 - stop_idx)); \
  vec = _mm512_maskz_loadu_epi64(mask, &sieve64[i]); \
  vec = _mm512_popcnt_epi64(vec); \
  vcnt = _mm512_add_epi64(vcnt, vec); \
  uint64_t cnt = _mm512_reduce_add_epi64(vcnt);

/// ARM SVE //////////////////////////////////////////////////////////

/// Count 1 bits inside [start, stop] using ARM SVE
#define SIEVE_COUNT_ARM_SVE(start, stop) \
  ASSERT(start <= stop); \
  ASSERT(stop - start < segment_size()); \
  uint64_t start_idx = start / 240; \
  uint64_t stop_idx = stop / 240; \
  uint64_t m1 = unset_smaller[start % 240]; \
  uint64_t m2 = unset_larger[stop % 240]; \
  \
  /* Branchfree bitmask calculation: */ \
  /* if (start_idx == stop_idx) m1 = m1 & m2; */ \
  /* if (start_idx == stop_idx) m2 = 0; */ \
  CONDITIONAL_MOVE(start_idx == stop_idx, m1, m1 & m2); \
  CONDITIONAL_MOVE(start_idx == stop_idx, m2, 0); \
  \
  const uint64_t* sieve64 = (const uint64_t*) sieve_.data(); \
  uint64_t start_bits = sieve64[start_idx] & m1; \
  uint64_t stop_bits = sieve64[stop_idx] & m2; \
  ASSERT(svcntd() >= 2); \
  svuint64_t vec = svinsr_n_u64(svdup_u64(start_bits), stop_bits); \
  svuint64_t vcnt = svcnt_u64_z(svwhilelt_b64(0, 2), vec); \
  uint64_t i = start_idx + 1; \
  \
  /* Compute this for loop using ARM SVE. */ \
  /* for (i = start_idx + 1; i < stop_idx; i++) */ \
  /*   cnt += popcnt64(sieve64[i]); */ \
  for (; i + svcntd() < stop_idx; i += svcntd()) \
  { \
    vec = svld1_u64(svptrue_b64(), &sieve64[i]); \
    vec = svcnt_u64_x(svptrue_b64(), vec); \
    vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec); \
  } \
  svbool_t pg = svwhilelt_b64(i, stop_idx); \
  vec = svld1_u64(pg, &sieve64[i]); \
  vec = svcnt_u64_z(pg, vec); \
  vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec); \
  uint64_t cnt = svaddv_u64(svptrue_b64(), vcnt);

#endif

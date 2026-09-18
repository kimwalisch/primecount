///
/// @file  count_simd.hpp
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

#ifndef COUNT_SIMD_HPP
#define COUNT_SIMD_HPP

#include <cpu_arch_macros.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

#if defined(ENABLE_ARM_NEON)
  #include <arm_neon.h>
#endif
#if defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <arm_sve.h>
#endif
#if defined(ENABLE_AVX512_VPOPCNT) || \
    defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <immintrin.h>
#endif

#if defined(ENABLE_ARM_NEON) && \
    defined(ENABLE_COUNT_DEFAULT)

/// ARM NEON /////////////////////////////////////////////////////////

/// Count 1 bits inside [start, stop] using ARM NEON
#define SIEVE_COUNT_DEFAULT(start, stop) \
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
  const uint64_t* sieve = sieve_.data(); \
  uint64_t start_bits = sieve[start_idx] & m1; \
  uint64_t stop_bits = sieve[stop_idx] & m2; \
  uint64x2_t bounds = vdupq_n_u64(start_bits); \
  bounds = vsetq_lane_u64(stop_bits, bounds, 1); \
  uint8x16_t bounds_cnt8 = vcntq_u8(vreinterpretq_u8_u64(bounds)); \
  uint16x8_t bounds_cnt16 = vpaddlq_u8(bounds_cnt8); \
  uint32x4_t bounds_cnt32 = vpaddlq_u16(bounds_cnt16); \
  uint64x2_t vcnt = vpaddlq_u32(bounds_cnt32); \
  uint64_t i = start_idx + 1; \
  \
  /* Compute this for loop using ARM NEON. */ \
  /* for (i = start_idx + 1; i < stop_idx; i++) */ \
  /*   cnt += popcnt64(sieve[i]); */ \
  NO_UNROLL_LOOP \
  for (; i + 2 <= stop_idx; i += 2) \
  { \
    uint64x2_t vec = vld1q_u64(&sieve[i]); \
    uint8x16_t cnt8 = vcntq_u8(vreinterpretq_u8_u64(vec)); \
    uint16x8_t cnt16 = vpaddlq_u8(cnt8); \
    uint32x4_t cnt32 = vpaddlq_u16(cnt16); \
    vcnt = vaddq_u64(vcnt, vpaddlq_u32(cnt32)); \
  } \
  /* Branchfree computation of: */ \
  /* if (i < stop_idx) */ \
  /*   cnt += popcnt64(sieve[i]); */ \
  uint64_t tail_bits = sieve[i & -(i < stop_idx)] & -(i < stop_idx); \
  uint64x2_t vec = vsetq_lane_u64(tail_bits, vdupq_n_u64(0), 0); \
  uint8x16_t cnt8 = vcntq_u8(vreinterpretq_u8_u64(vec)); \
  uint16x8_t cnt16 = vpaddlq_u8(cnt8); \
  uint32x4_t cnt32 = vpaddlq_u16(cnt16); \
  vcnt = vaddq_u64(vcnt, vpaddlq_u32(cnt32)); \
  uint64_t cnt = vaddvq_u64(vcnt);

#elif defined(ENABLE_COUNT_DEFAULT)

/// POPCNT64 /////////////////////////////////////////////////////////

/// Count 1 bits inside [start, stop] using POPCNT64
#define SIEVE_COUNT_DEFAULT(start, stop) \
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
  const uint64_t* sieve = sieve_.data(); \
  uint64_t start_bits = sieve[start_idx] & m1; \
  uint64_t stop_bits = sieve[stop_idx] & m2; \
  uint64_t cnt = popcnt64(start_bits); \
  cnt += popcnt64(stop_bits); \
  \
  NO_UNROLL_LOOP \
  for (uint64_t i = start_idx + 1; i < stop_idx; i++) \
    cnt += popcnt64(sieve[i]);

#endif

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
  const uint64_t* sieve = sieve_.data(); \
  uint64_t start_bits = sieve[start_idx] & m1; \
  uint64_t stop_bits = sieve[stop_idx] & m2; \
  uint64_t cnt = popcnt64_native(start_bits); \
  cnt += popcnt64_native(stop_bits); \
  __m512i vcnt = _mm512_setzero_si512(); \
  uint64_t i = start_idx + 1; \
  \
  /* Compute this for loop using AVX512. */ \
  /* for (i = start_idx + 1; i < stop_idx; i++) */ \
  /*   cnt += popcnt64(sieve[i]); */ \
  NO_UNROLL_LOOP \
  for (; i + 8 < stop_idx; i += 8) \
  { \
    __m512i vec = _mm512_loadu_epi64(&sieve[i]); \
    vec = _mm512_popcnt_epi64(vec); \
    vcnt = _mm512_add_epi64(vcnt, vec); \
  } \
  __mmask8 mask = (__mmask8) (0xff >> (i + 8 - stop_idx)); \
  __m512i vec = _mm512_maskz_loadu_epi64(mask, &sieve[i]); \
  vec = _mm512_popcnt_epi64(vec); \
  vcnt = _mm512_add_epi64(vcnt, vec); \
  cnt += _mm512_reduce_add_epi64(vcnt);

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
  const uint64_t* sieve = sieve_.data(); \
  uint64_t start_bits = sieve[start_idx] & m1; \
  uint64_t stop_bits = sieve[stop_idx] & m2; \
  ASSERT(svcntd() >= 2); \
  svuint64_t bounds = svinsr_n_u64(svdup_u64(start_bits), stop_bits); \
  bounds = svcnt_u64_z(svwhilelt_b64(0, 2), bounds); \
  svuint64_t vcnt = svdup_u64(0); \
  uint64_t i = start_idx + 1; \
  \
  /* Compute this for loop using ARM SVE. */ \
  /* for (i = start_idx + 1; i < stop_idx; i++) */ \
  /*   cnt += popcnt64(sieve[i]); */ \
  NO_UNROLL_LOOP \
  for (; i + svcntd() < stop_idx; i += svcntd()) \
  { \
    svuint64_t vec = svld1_u64(svptrue_b64(), &sieve[i]); \
    vec = svcnt_u64_x(svptrue_b64(), vec); \
    vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec); \
  } \
  svbool_t pg = svwhilelt_b64(i, stop_idx); \
  svuint64_t vec = svld1_u64(pg, &sieve[i]); \
  vec = svcnt_u64_z(pg, vec); \
  vcnt = svadd_u64_x(svptrue_b64(), vcnt, vec); \
  vcnt = svadd_u64_x(svptrue_b64(), vcnt, bounds); \
  uint64_t cnt = svaddv_u64(svptrue_b64(), vcnt);

#endif

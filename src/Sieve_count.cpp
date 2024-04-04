///
/// @file  Sieve_count.cpp
/// @brief Count the number of 1 bits in the sieve array. Since fast
///        bit counting is very important for primecount's performance
///        we use vector instructions to speed up the computation.
///
///        In-depth description of our counting algorithms:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.md
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <Sieve.hpp>
#include <SieveTables.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>

// x64 AVX512 vector popcount
#if __has_include(<immintrin.h>) && \
   (defined(__AVX512__) || (defined(__AVX512F__) && \
                            defined(__AVX512VPOPCNTDQ__)))
  #include <immintrin.h>
  #define HAS_AVX512_VPOPCNT

// GCC/Clang function multiversioning
#elif defined(MULTIARCH_TARGET_AVX512) && \
    __has_include(<immintrin.h>)
  #include <immintrin.h>

// ARM SVE vector popcount
#elif defined(__ARM_FEATURE_SVE) && \
      __has_include(<arm_sve.h>)
  #include <arm_sve.h>
  #define HAS_ARM_SVE

#else // Default portable popcount
  #define DEFAULT_CPU_ARCH
#endif

/// Count 1 bits inside [start, stop]
#define COUNT_1_BITS(start, stop, POPCNT_KERNEL)      \
  uint64_t res = 0;                                   \
                                                      \
  if (start <= stop)                                  \
  {                                                   \
    ASSERT(stop - start < segment_size());            \
    uint64_t start_idx = start / 240;                 \
    uint64_t stop_idx = stop / 240;                   \
    uint64_t m1 = unset_smaller[start % 240];         \
    uint64_t m2 = unset_larger[stop % 240];           \
    uint64_t* sieve64 = (uint64_t*) sieve_.data();    \
                                                      \
    if (start_idx == stop_idx)                        \
      res = popcnt64(sieve64[start_idx] & (m1 & m2)); \
    else                                              \
    {                                                 \
      res = popcnt64(sieve64[start_idx] & m1);        \
      POPCNT_KERNEL(sieve64, start_idx, stop_idx)     \
      res += popcnt64(sieve64[stop_idx] & m2);        \
    }                                                 \
  }

/// Default portable POPCNT kernel
#define DEFAULT_POPCNT_KERNEL(sieve64, start_idx, stop_idx)  \
  for (uint64_t i = start_idx + 1; i < stop_idx; i++)        \
    res += popcnt64(sieve64[i]);

/// Compute the loop below using ARM SVE.
/// for (i = start_idx + 1; i < stop_idx; i++)
///   res += popcnt64(sieve64[i]);
///
#define SVE_POPCNT_KERNEL(sieve64, start_idx, stop_idx)      \
  uint64_t i = start_idx + 1;                                \
  svuint64_t vcnt = svdup_u64(0);                            \
  for (; i < stop_idx; i += svcntd())                        \
  {                                                          \
    svbool_t pg = svwhilelt_b64(i, stop_idx);                \
    svuint64_t vec = svld1_u64(pg, &sieve64[i]);             \
    vec = svcnt_u64_z(pg, vec);                              \
    vcnt = svadd_u64_z(svptrue_b64(), vcnt, vec);            \
  }                                                          \
  res += svaddv_u64(svptrue_b64(), vcnt);

/// Compute the loop below using AVX512.
/// for (i = start_idx + 1; i < stop_idx; i++)
///   res += popcnt64(sieve64[i]);
///
#define AVX512_POPCNT_KERNEL(sieve64, start_idx, stop_idx)   \
  uint64_t i = start_idx + 1;                                \
  __m512i vcnt = _mm512_setzero_si512();                     \
  for (; i + 8 < stop_idx; i += 8)                           \
  {                                                          \
    __m512i vec = _mm512_loadu_epi64(&sieve64[i]);           \
    vec = _mm512_popcnt_epi64(vec);                          \
    vcnt = _mm512_add_epi64(vcnt, vec);                      \
  }                                                          \
  __mmask8 mask = 0xff >> (i + 8 - stop_idx);                \
  __m512i vec = _mm512_maskz_loadu_epi64(mask, &sieve64[i]); \
  vec = _mm512_popcnt_epi64(vec);                            \
  vcnt = _mm512_add_epi64(vcnt, vec);                        \
  res += _mm512_reduce_add_epi64(vcnt);

namespace primecount {

#if defined(DEFAULT_CPU_ARCH) || \
    defined(MULTIARCH_TARGET_DEFAULT)

#if defined(MULTIARCH_TARGET_DEFAULT)
  __attribute__ ((target ("default")))
#endif
uint64_t Sieve::count(uint64_t stop)
{
  // Count 1 bits inside [0, stop] in the sieve array.
  // The distance [0, stop] might be large > sqrt(segment_size).
  ASSERT(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

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
  COUNT_1_BITS(start, stop, DEFAULT_POPCNT_KERNEL)
  count_ += res;

  return count_;
}

#if defined(MULTIARCH_TARGET_DEFAULT)
  __attribute__ ((target ("default")))
#endif
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  // Count 1 bits inside [start, stop] in the sieve array.
  // The distance [start, stop] is small here < sqrt(segment_size).
  COUNT_1_BITS(start, stop, DEFAULT_POPCNT_KERNEL)
  return res;
}

#endif

#if defined(HAS_AVX512_VPOPCNT) || \
    defined(MULTIARCH_TARGET_AVX512)

#if defined(MULTIARCH_TARGET_AVX512)
  __attribute__ ((target ("avx512f,avx512vpopcntdq")))
#endif
uint64_t Sieve::count(uint64_t stop)
{
  // Count 1 bits inside [0, stop] in the sieve array.
  // The distance [0, stop] might be large > sqrt(segment_size).
  ASSERT(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

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
  COUNT_1_BITS(start, stop, AVX512_POPCNT_KERNEL)
  count_ += res;

  return count_;
}

#if defined(MULTIARCH_TARGET_AVX512)
  __attribute__ ((target ("avx512f,avx512vpopcntdq")))
#endif
uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  // Count 1 bits inside [start, stop] in the sieve array.
  // The distance [start, stop] is small here < sqrt(segment_size).
  COUNT_1_BITS(start, stop, AVX512_POPCNT_KERNEL)
  return res;
}

#elif defined(HAS_ARM_SVE)

uint64_t Sieve::count(uint64_t stop)
{
  // Count 1 bits inside [0, stop] in the sieve array.
  // The distance [0, stop] might be large > sqrt(segment_size).
  ASSERT(stop >= prev_stop_);
  uint64_t start = prev_stop_ + 1;
  prev_stop_ = stop;

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
  COUNT_1_BITS(start, stop, SVE_POPCNT_KERNEL)
  count_ += res;

  return count_;
}

uint64_t Sieve::count(uint64_t start, uint64_t stop) const
{
  // Count 1 bits inside [start, stop] in the sieve array.
  // The distance [start, stop] is small here < sqrt(segment_size).
  COUNT_1_BITS(start, stop, SVE_POPCNT_KERNEL)
  return res;
}

#endif

} // namespace

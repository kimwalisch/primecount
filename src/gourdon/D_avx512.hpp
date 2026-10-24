///
/// @file  D_avx512.hpp
/// @brief AVX512 implementation of the D formula (hard special
///        leaves) in Xavier Gourdon's prime counting algorithm. This
///        algorithm is identical to D_thread_default() in D_default.hpp
///        except that this algorithm has been partially vectorized
///        using AVX512.
///
///        For performance it is important that all AVX512 helper
///        functions are inlined by the compiler. We achieve this
///        by annotating all AVX512 helper functions using the
///        ALWAYS_INLINE macro and PRIMECOUNT_TARGET_ATTRIBUTE with
///        the same CPU instruction sets as PRIMECOUNT_TARGET_ENTRY.
///
///        primecount's D formula differs from Gourdon's original
///        D formula: it also computes part of the easy special
///        leaves that belong to the C formula in Gourdon's paper.
///        These are the leaves with (x/y)^(1/3) < p <= sqrt(y)
///        and x/p^3 < m <= y, where p is a prime and m is a prime
///        or a product of 2 distinct primes, all greater than p.
///        They are computed together with the hard special leaves
///        using the hyperbola method.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves-SIMD-Filtering.pdf
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_AVX512_HPP
#define D_AVX512_HPP

#include "DHyperbola.hpp"

#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <LoadBalancerS2.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <phi_vector.hpp>
#include <PiTable.hpp>
#include <sieve/Sieve.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <immintrin.h>

namespace {

using namespace primecount;

// Not a lambda, operator() requires the AVX512 target attribute
struct SieveCountAvx512
{
  Sieve& sieve;

  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
  #endif
  inline uint64_t operator()(uint64_t stop) const
  {
    return sieve.count_avx512(stop);
  }
};

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_epi32_avx512(const uint16_t* factor_table,
                                               __m512i reverse32)
{
  __m256i vec = _mm256_loadu_si256((const __m256i*) factor_table);
  return _mm512_permutexvar_epi32(reverse32, _mm512_cvtepu16_epi32(vec));
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_epi32_avx512(const uint32_t* factor_table,
                                               __m512i reverse32)
{
  __m512i vec = _mm512_loadu_si512((const void*) factor_table);
  return _mm512_permutexvar_epi32(reverse32, vec);
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_tail_epi32_avx512(const uint16_t* factor_table,
                                                    __mmask16 load_mask,
                                                    int count,
                                                    __m512i m_offsets32)
{
  __m256i vec = _mm256_maskz_loadu_epi16(load_mask, factor_table);
  __m512i reverse32 = _mm512_sub_epi32(_mm512_set1_epi32(count - 1), m_offsets32);
  return _mm512_maskz_permutexvar_epi32(load_mask,
                                        reverse32,
                                        _mm512_cvtepu16_epi32(vec));
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_tail_epi32_avx512(const uint32_t* factor_table,
                                                    __mmask16 load_mask,
                                                    int count,
                                                    __m512i m_offsets32)
{
  __m512i vec = _mm512_maskz_loadu_epi32(load_mask, factor_table);
  __m512i reverse32 = _mm512_sub_epi32(_mm512_set1_epi32(count - 1), m_offsets32);
  return _mm512_maskz_permutexvar_epi32(load_mask, reverse32, vec);
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_epi64_avx512(const uint16_t* factor_table,
                                               __m512i reverse64)
{
  __m128i vec = _mm_loadu_si128((const __m128i*) factor_table);
  return _mm512_permutexvar_epi64(reverse64, _mm512_cvtepu16_epi64(vec));
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_epi64_avx512(const uint32_t* factor_table,
                                               __m512i reverse64)
{
  __m256i vec = _mm256_loadu_si256((const __m256i*) factor_table);
  return _mm512_permutexvar_epi64(reverse64, _mm512_cvtepu32_epi64(vec));
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_tail_epi64_avx512(const uint16_t* factor_table,
                                                    __mmask8 load_mask,
                                                    int count,
                                                    __m512i m_offsets64)
{
  __m128i vec = _mm_maskz_loadu_epi16(load_mask, factor_table);
  __m512i reverse64 = _mm512_sub_epi64(_mm512_set1_epi64(count - 1), m_offsets64);
  return _mm512_maskz_permutexvar_epi64(load_mask,
                                        reverse64,
                                        _mm512_cvtepu16_epi64(vec));
}

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
ALWAYS_INLINE __m512i load_factor_tail_epi64_avx512(const uint32_t* factor_table,
                                                    __mmask8 load_mask,
                                                    int count,
                                                    __m512i m_offsets64)
{
  __m256i vec = _mm256_maskz_loadu_epi32(load_mask, factor_table);
  __m512i reverse64 = _mm512_sub_epi64(_mm512_set1_epi64(count - 1), m_offsets64);
  return _mm512_maskz_permutexvar_epi64(load_mask,
                                        reverse64,
                                        _mm512_cvtepu32_epi64(vec));
}

template <typename T, typename Primes, typename FactorTable>
#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  PRIMECOUNT_TARGET_ENTRY("avx512f,avx512bw,avx512vl,avx512vpopcntdq,popcnt")
#endif
T D_thread_avx512(T x,
                  int64_t y,
                  int64_t xy,
                  int64_t k,
                  int64_t x14,
                  const Primes& primes,
                  const PiTable& pi,
                  const FactorTable& factor,
                  const Vector<int64_t>& mertens_y,
                  ThreadData& thread)
{
  int64_t low = thread.low;
  int64_t low1 = max(low, 1);
  int64_t segments = thread.segments;
  int64_t segment_size = thread.segment_size;
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t limit = min(low + segment_size * segments, xy);
  int64_t max_b = pi[min3(isqrt(x / low1), isqrt(limit), x14)];
  int64_t min_b = pi[min(xy / limit, x14)];
  min_b = max(k, min_b) + 1;

  if (min_b > max_b)
    return 0;

  Vector<int64_t> phi = phi_vector(low, max_b, primes, pi);
  Sieve sieve(low, segment_size, max_b);
  Vector<SemiprimeEvent> events;
  SemiprimeEventSorter sorter;
  SieveCountAvx512 count_sieve{sieve};
  thread.init_time = get_time();

  INDETERMINATE Array<uint32_t, 128> m_indexes32;
  INDETERMINATE Array< int64_t, 128> m_indexes64;
  INDETERMINATE Array< int64_t, 128> xpm_low;
  const auto* factor_table = factor.data();

  __m512i reverse32 = _mm512_setr_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
  __m512i reverse64 = _mm512_setr_epi64(7, 6, 5, 4, 3, 2, 1, 0);
  __m512i m_offsets32 = _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
  __m512i m_offsets64 = _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 6, 7);

  T sum = 0;

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    low1 = max(low, 1);

    // For b < min_b there are no special leaves:
    // low <= x / (primes[b] * m) < high
    sieve.pre_sieve(primes, min_b - 1, low, high);
    sieve.init_counter(low, high);
    int64_t b = min_b;

    // For k + 1 <= b <= pi_sqrty
    // Find all special leaves in the current segment that are
    // composed of a prime and a square free number:
    // low <= x / (primes[b] * m) < high
    for (int64_t last = min(pi_sqrty, max_b); b <= last; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t divisor = max(low1, prime * prime);
      int64_t min_m = min(fast_div(xp, high), y);
      int64_t max_m = min(fast_div(xp, divisor), y);

      if (prime >= max_m)
        goto next_segment;

      // A composite m with all prime factors > prime exceeds
      // prime^2 and all hyperbola terms have positions < y.
      if (max_m <= prime * prime && low >= y)
        break;

      if (low < y)
        sum += hyperbola_sum(xp, prime, b, y, low, high, phi[b],
                             mertens_y[b], primes, pi, count_sieve,
                             events, sorter, max_m);

      min_m = FactorTable::to_index(min_m);
      max_m = FactorTable::to_index(max_m);
      int64_t encoded_prime = FactorTable::encode(b);
      int64_t m = max_m;
      std::size_t m_count = 0;

      // AVX512: 16-lane 32-bit
      if (max_m <= UINT32_MAX ||
          sizeof(T) <= sizeof(uint64_t))
      {
        __m512i encoded_prime_vec = _mm512_set1_epi32(uint32_t(encoded_prime));
        constexpr std::size_t max_m_count = m_indexes32.size() - 16;

        for (; m >= min_m + 16; m -= 16)
        {
          // Filter out square free m values using AVX512
          // that satisfy: factor_table[m] > encoded_prime
          __m512i m_vec = _mm512_sub_epi32(_mm512_set1_epi32(uint32_t(m)), m_offsets32);
          __m512i factor_vec = load_factor_epi32_avx512(&factor_table[m - 15], reverse32);
          __mmask16 mask = _mm512_cmpgt_epu32_mask(factor_vec, encoded_prime_vec);
          _mm512_mask_compressstoreu_epi32(&m_indexes32[m_count], mask, m_vec);
          m_count += _mm_popcnt_u64(mask);

          if (m_count > max_m_count)
          {
            // Batch calculate (xp/m - low) to improve CPU pipelining
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = factor.to_number(m_indexes32[i]);
              xpm_low[i] = fast_div64(xp, m) - low;
            }

            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              // sieve.count(xp/m - low)
              int64_t count = sieve.count_avx512(xpm_low[i]);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m_indexes32[i]) * phi_xpm;
            }

            m_count = 0;
          }
        }

        // Filter out last few square free m
        if (m > min_m)
        {
          int count = int(m - min_m);
          __mmask16 load_mask = __mmask16((1u << count) - 1);
          __m512i m_vec = _mm512_sub_epi32(_mm512_set1_epi32(uint32_t(m)), m_offsets32);
          __m512i factor_vec = load_factor_tail_epi32_avx512(&factor_table[min_m + 1], load_mask, count, m_offsets32);
          __mmask16 mask = _mm512_cmpgt_epu32_mask(factor_vec, encoded_prime_vec);
          _mm512_mask_compressstoreu_epi32(&m_indexes32[m_count], mask, m_vec);
          m_count += _mm_popcnt_u64(mask);
        }

        // Batch calculate (xp/m - low) to improve CPU pipelining
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = factor.to_number(m_indexes32[i]);
          xpm_low[i] = fast_div64(xp, m) - low;
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          // sieve.count(xp/m - low)
          int64_t count = sieve.count_avx512(xpm_low[i]);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m_indexes32[i]) * phi_xpm;
        }
      }
      else // AVX512: 8-lane 64-bit
      {
        __m512i encoded_prime_vec = _mm512_set1_epi64(encoded_prime);
        constexpr std::size_t max_m_count = m_indexes64.size() - 8;

        for (; m >= min_m + 8; m -= 8)
        {
          // Filter out square free m values using AVX512
          // that satisfy: factor_table[m] > encoded_prime
          __m512i m_vec = _mm512_sub_epi64(_mm512_set1_epi64(m), m_offsets64);
          __m512i factor_vec = load_factor_epi64_avx512(&factor_table[m - 7], reverse64);
          __mmask8 mask = _mm512_cmpgt_epi64_mask(factor_vec, encoded_prime_vec);
          _mm512_mask_compressstoreu_epi64(&m_indexes64[m_count], mask, m_vec);
          m_count += _mm_popcnt_u64(mask);

          if (m_count > max_m_count)
          {
            // Batch calculate (xp/m - low) to improve CPU pipelining
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = factor.to_number(m_indexes64[i]);
              xpm_low[i] = fast_div64(xp, m) - low;
            }

            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              // sieve.count(xp/m - low)
              int64_t count = sieve.count_avx512(xpm_low[i]);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m_indexes64[i]) * phi_xpm;
            }

            m_count = 0;
          }
        }

        // Filter out last few square free m
        if (m > min_m)
        {
          int count = int(m - min_m);
          __mmask8 load_mask = __mmask8((1u << count) - 1);
          __m512i m_vec = _mm512_sub_epi64(_mm512_set1_epi64(m), m_offsets64);
          __m512i factor_vec = load_factor_tail_epi64_avx512(&factor_table[min_m + 1], load_mask, count, m_offsets64);
          __mmask8 mask = _mm512_cmpgt_epi64_mask(factor_vec, encoded_prime_vec);
          _mm512_mask_compressstoreu_epi64(&m_indexes64[m_count], mask, m_vec);
          m_count += _mm_popcnt_u64(mask);
        }

        // Batch calculate (xp/m - low) to improve CPU pipelining
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = factor.to_number(m_indexes64[i]);
          xpm_low[i] = fast_div64(xp, m) - low;
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          // sieve.count(xp/m - low)
          int64_t count = sieve.count_avx512(xpm_low[i]);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m_indexes64[i]) * phi_xpm;
        }
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    // For the remaining b <= pi_x14
    // Find all special leaves in the current segment
    // that are composed of 2 primes:
    // low <= x / (primes[b] * primes[l]) < high
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t divisor = max(low1, prime * prime);
      int64_t xp_high = min(fast_div(xp, high), y);
      int64_t min_m = max(xp_high, prime);
      int64_t max_m = min(fast_div(xp, divisor), y);
      int64_t l = pi[max_m];

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xpq = fast_div64(xp, primes[l]);
        int64_t count = sieve.count_avx512(xpq - low);
        int64_t phi_xpq = phi[b] + count;
        sum += phi_xpq;
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    next_segment:;
  }

  return sum;
}

} // namespace

#endif

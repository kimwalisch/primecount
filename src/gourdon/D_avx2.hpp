///
/// @file  D_avx2.hpp
/// @brief AVX2 implementation of the D formula (hard special leaves)
///        in Xavier Gourdon's prime counting algorithm. This algorithm
///        is identical to D_thread_default() in D_default.hpp except for
///        AVX2 vectorization of index filtering.
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

#ifndef D_AVX2_HPP
#define D_AVX2_HPP

#include "DCentral.hpp"

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

/// AVX2 has no compress instruction, but since the filtered values
/// are the factor table indexes themselves, no permutation is
/// needed: the selected indexes are m - offsets[mask]. Bit i of
/// the mask corresponds to the offset 7 - i. The offsets are
/// stored in ascending order (descending indexes), unused slots
/// remain zero.
///
struct Avx2CompactOffsets
{
  alignas(64) uint8_t values[256][8];

  constexpr Avx2CompactOffsets() : values{}
  {
    for (int mask = 0; mask < 256; mask++)
    {
      int count = 0;
      for (int offset = 0; offset < 8; offset++)
        if (mask & (1 << (7 - offset)))
          values[mask][count++] = uint8_t(offset);
    }
  }
};

constexpr Avx2CompactOffsets avx2_compact_offsets;

/// Bit i of the returned mask is set if
/// factor_table[i] > encoded_prime.
///
#if defined(ENABLE_MULTIARCH_AVX2)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx2,popcnt")
#endif
ALWAYS_INLINE uint32_t mask8_avx2(const uint16_t* factor_table,
                                  uint32_t encoded_prime)
{
  // Zero extended 16-bit entries compare
  // correctly as signed 32-bit integers.
  ASSERT(encoded_prime <= UINT16_MAX);
  __m256i vec = _mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i*) factor_table));
  __m256i cmp = _mm256_cmpgt_epi32(vec, _mm256_set1_epi32(encoded_prime));
  return _mm256_movemask_ps(_mm256_castsi256_ps(cmp));
}

#if defined(ENABLE_MULTIARCH_AVX2)
  PRIMECOUNT_TARGET_ATTRIBUTE("avx2,popcnt")
#endif
ALWAYS_INLINE uint32_t mask8_avx2(const uint32_t* factor_table,
                                  uint32_t encoded_prime)
{
  // The factor table contains the unsigned sentinel values
  // UINT32_MAX - 1 and UINT32_MAX, hence we must use an
  // unsigned comparison: a > b <=> max(a, b + 1) == a.
  ASSERT(encoded_prime < UINT32_MAX);
  __m256i vec = _mm256_loadu_si256((const __m256i*) factor_table);
  __m256i vmax = _mm256_max_epu32(vec, _mm256_set1_epi32(encoded_prime + 1));
  __m256i cmp = _mm256_cmpeq_epi32(vmax, vec);
  return _mm256_movemask_ps(_mm256_castsi256_ps(cmp));
}

/// Compute the contribution of the hard special leaves using
/// a segmented sieve. Each thread processes the interval
/// [low, low + segment_size * segments[.
///
template <typename T, typename Primes, typename FactorTable>
#if defined(ENABLE_MULTIARCH_AVX2)
  PRIMECOUNT_TARGET_ENTRY("avx2,popcnt")
#endif
T D_thread_avx2(T x,
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
  T sum = 0;

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
  Vector<CentralEvent> events;
  SemiprimeEventBitmap sorter;
  auto count_sieve = [&sieve](uint64_t stop) { return sieve.count(stop); };
  thread.init_time = get_time();

  INDETERMINATE Array<uint32_t, 128> m_indexes32;
  INDETERMINATE Array< int64_t, 128> m_indexes64;
  INDETERMINATE Array< int64_t, 128> xpm_low;
  const auto* factor_table = factor.data();

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
      int64_t xp_low = min(fast_div(xp, low1), y);
      int64_t xp_high = min(fast_div(xp, high), y);
      int64_t min_m = max(xp_high, y / prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);
      int64_t original_max_m = max_m;

      if (prime >= max_m)
        goto next_segment;

      if (central_row(xp, prime, b, y, low, high, phi[b], primes, pi,
                      mertens_y, count_sieve, events, sorter, sum, max_m))
      {
        phi[b] += sieve.get_total_count();
        sieve.cross_off_count(prime, b);
        continue;
      }

      // A composite m with all prime factors > prime exceeds prime^2.
      // Test the uncapped bound: later rows may still have central
      // leaves even when the capped max_m is at most prime^2.
      if (original_max_m <= prime * prime && low >= y)
        break;

      min_m = FactorTable::to_index(min_m);
      max_m = FactorTable::to_index(max_m);
      int64_t encoded_prime = FactorTable::encode(b);
      int64_t m = max_m;
      std::size_t m_count = 0;

      // AVX2 32-bit
      if (max_m <= UINT32_MAX ||
          sizeof(T) <= sizeof(uint64_t))
      {
        ASSERT(encoded_prime < UINT32_MAX);
        constexpr std::size_t max_m_count = m_indexes32.size() - 8;
        __m256i m_vec = _mm256_set1_epi32(uint32_t(m));
        __m256i lane_step = _mm256_set1_epi32(8);

        // Filter out square free m values using AVX2
        // that satisfy: factor_table[m] > encoded_prime
        for (; m >= min_m + 8; m -= 8)
        {
          uint32_t mask = mask8_avx2(&factor_table[m - 7], uint32_t(encoded_prime));
          __m128i offsets = _mm_loadl_epi64((const __m128i*) avx2_compact_offsets.values[mask]);
          __m256i indexes = _mm256_sub_epi32(m_vec, _mm256_cvtepu8_epi32(offsets));
          _mm256_storeu_si256((__m256i*) &m_indexes32[m_count], indexes);
          m_count += _mm_popcnt_u32(mask);
          m_vec = _mm256_sub_epi32(m_vec, lane_step);

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
              int64_t count = sieve.count(xpm_low[i]);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m_indexes32[i]) * phi_xpm;
            }

            m_count = 0;
          }
        }

        // Filter out the last few square free m
        for (; m > min_m; m--)
        {
          m_indexes32[m_count] = uint32_t(m);
          m_count += (factor_table[m] > encoded_prime);
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
          int64_t count = sieve.count(xpm_low[i]);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m_indexes32[i]) * phi_xpm;
        }
      }
      else // AVX2 64-bit
      {
        ASSERT(encoded_prime < UINT32_MAX);
        constexpr std::size_t max_m_count = m_indexes64.size() - 8;
        __m256i m_vec = _mm256_set1_epi64x(m);
        __m256i lane_step = _mm256_set1_epi64x(8);

        // Filter out square free m values using AVX2
        // that satisfy: factor_table[m] > encoded_prime
        for (; m >= min_m + 8; m -= 8)
        {
          uint32_t mask = mask8_avx2(&factor_table[m - 7], uint32_t(encoded_prime));
          __m128i offsets = _mm_loadl_epi64((const __m128i*) avx2_compact_offsets.values[mask]);
          __m256i indexes = _mm256_sub_epi64(m_vec, _mm256_cvtepu8_epi64(offsets));
          _mm256_storeu_si256((__m256i*) &m_indexes64[m_count], indexes);
          offsets = _mm_srli_si128(offsets, 4);
          indexes = _mm256_sub_epi64(m_vec, _mm256_cvtepu8_epi64(offsets));
          _mm256_storeu_si256((__m256i*) &m_indexes64[m_count + 4], indexes);
          m_count += _mm_popcnt_u32(mask);
          m_vec = _mm256_sub_epi64(m_vec, lane_step);

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
              int64_t count = sieve.count(xpm_low[i]);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m_indexes64[i]) * phi_xpm;
            }

            m_count = 0;
          }
        }

        // Filter out the last few square free m
        for (; m > min_m; m--)
        {
          m_indexes64[m_count] = m;
          m_count += (factor_table[m] > encoded_prime);
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
          int64_t count = sieve.count(xpm_low[i]);
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
      int64_t xp_low = min(fast_div(xp, low1), y);
      int64_t xp_high = min(fast_div(xp, high), y);
      int64_t min_m = max3(xp_high, y / prime, prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);
      int64_t l = pi[max_m];

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xpq = fast_div64(xp, primes[l]);
        int64_t count = sieve.count(xpq - low);
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

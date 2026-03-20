///
/// @file  D_avx512.hpp
/// @brief AVX-512 implementation of the D(x, y) thread routine used
///        in Xavier Gourdon's prime counting algorithm.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_AVX512_HPP
#define D_AVX512_HPP

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_epi32_avx512(const uint16_t* factor_data,
                                 __m512i reverse32)
{
  __m256i vec = _mm256_loadu_si256((const __m256i*) factor_data);
  return _mm512_permutexvar_epi32(reverse32, _mm512_cvtepu16_epi32(vec));
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_epi32_avx512(const uint32_t* factor_data,
                                 __m512i reverse32)
{
  __m512i vec = _mm512_loadu_si512((const void*) factor_data);
  return _mm512_permutexvar_epi32(reverse32, vec);
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_tail_epi32_avx512(const uint16_t* factor_data,
                                      __mmask16 load_mask,
                                      int count,
                                      __m512i m_offsets32)
{
  __m256i vec = _mm256_maskz_loadu_epi16(load_mask, factor_data);
  __m512i reverse32 = _mm512_sub_epi32(_mm512_set1_epi32(count - 1), m_offsets32);
  return _mm512_maskz_permutexvar_epi32(load_mask,
                                        reverse32,
                                        _mm512_cvtepu16_epi32(vec));
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_tail_epi32_avx512(const uint32_t* factor_data,
                                      __mmask16 load_mask,
                                      int count,
                                      __m512i m_offsets32)
{
  __m512i vec = _mm512_maskz_loadu_epi32(load_mask, factor_data);
  __m512i reverse32 = _mm512_sub_epi32(_mm512_set1_epi32(count - 1), m_offsets32);
  return _mm512_maskz_permutexvar_epi32(load_mask, reverse32, vec);
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_epi64_avx512(const uint16_t* factor_data,
                                 __m512i reverse64)
{
  __m128i vec = _mm_loadu_si128((const __m128i*) factor_data);
  return _mm512_permutexvar_epi64(reverse64, _mm512_cvtepu16_epi64(vec));
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_epi64_avx512(const uint32_t* factor_data,
                                 __m512i reverse64)
{
  __m256i vec = _mm256_loadu_si256((const __m256i*) factor_data);
  return _mm512_permutexvar_epi64(reverse64, _mm512_cvtepu32_epi64(vec));
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_tail_epi64_avx512(const uint16_t* factor_data,
                                      __mmask8 load_mask,
                                      int count,
                                      __m512i m_offsets64)
{
  __m128i vec = _mm_maskz_loadu_epi16(load_mask, factor_data);
  __m512i reverse64 = _mm512_sub_epi64(_mm512_set1_epi64(count - 1), m_offsets64);
  return _mm512_maskz_permutexvar_epi64(load_mask,
                                        reverse64,
                                        _mm512_cvtepu16_epi64(vec));
}

ALWAYS_INLINE __attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
__m512i load_factor_tail_epi64_avx512(const uint32_t* factor_data,
                                      __mmask8 load_mask,
                                      int count,
                                      __m512i m_offsets64)
{
  __m256i vec = _mm256_maskz_loadu_epi32(load_mask, factor_data);
  __m512i reverse64 = _mm512_sub_epi64(_mm512_set1_epi64(count - 1), m_offsets64);
  return _mm512_maskz_permutexvar_epi64(load_mask,
                                        reverse64,
                                        _mm512_cvtepu32_epi64(vec));
}

/// This algorithm computes the hard special leaves in the Gourdon
/// prime counting algorithm. This algorithm is identical to
/// D_thread_default() except that parts of the algorithm have
/// been vectorized using AVX512.
///
/// For performance it is important that all AVX512 helper functions
/// are inlined by the compiler. We achieve this by annotating all
/// AVX512 helper functions using the same AVX512 __attribute__ and
/// the ALWAYS_INLINE macro.
///
template <typename T, typename Primes, typename FactorTableD>
__attribute__ ((target ("avx512f,avx512bw,avx512vl,avx512vpopcntdq")))
T D_thread_avx512(T x,
                  int64_t x_star,
                  int64_t xz,
                  int64_t y,
                  int64_t z,
                  int64_t k,
                  const Primes& primes,
                  const PiTable& pi,
                  const FactorTableD& factor,
                  ThreadData& thread)
{
  int64_t low = thread.low;
  int64_t low1 = max(low, 1);
  int64_t segments = thread.segments;
  int64_t segment_size = thread.segment_size;
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t limit = min(low + segment_size * segments, xz);
  int64_t max_b = pi[min3(isqrt(x / low1), isqrt(limit), x_star)];
  int64_t min_b = pi[min(xz / limit, x_star)];
  min_b = max(k, min_b) + 1;

  if (min_b > max_b)
    return 0;

  Vector<int64_t> phi = phi_vector(low, max_b, primes, pi);
  Sieve sieve(low, segment_size, max_b);
  thread.init_finished();

  Array<uint32_t, 256> m_indexes32;
  Array< int64_t, 128> m_indexes64;
  const auto* factor_data = factor.factor_data();

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

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves in the current segment that are
    // composed of a prime and a square free number:
    // low <= x / (primes[b] * m) < high
    for (int64_t last = min(pi_sqrtz, max_b); b <= last; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_low = min(fast_div(xp, low1), z);
      int64_t xp_high = min(fast_div(xp, high), z);
      int64_t min_m = max(xp_high, z / prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);

      if (prime >= max_m)
        goto next_segment;

      min_m = factor.to_index(min_m);
      max_m = factor.to_index(max_m);
      std::size_t m_count = 0;
      int64_t m = max_m;

      if (FactorTableD::max() <= UINT32_MAX ||
          max_m <= UINT32_MAX)
      {
        __m512i prime_vec = _mm512_set1_epi32(uint32_t(prime));
        constexpr std::size_t max_m_count = m_indexes32.size() - 31;

        // AVX512: 16-lane 32-bit
        for (; m > min_m + 15; m -= 16)
        {
          if (m_count > max_m_count)
          {
            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = m_indexes32[i];
              int64_t xpm = fast_div64(xp, factor.to_number(m));
              int64_t count = sieve.count_avx512(xpm - low);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m) * phi_xpm;
            }
            m_count = 0;
          }

          // Filter out square free m values using AVX512
          // that satisfy: prime < factor.is_leaf(m)
          __m512i m_vec = _mm512_sub_epi32(_mm512_set1_epi32(uint32_t(m)), m_offsets32);
          __m512i factor_vec = load_factor_epi32_avx512(&factor_data[m - 15], reverse32);
          __mmask16 mask = _mm512_cmpgt_epu32_mask(factor_vec, prime_vec);
          _mm512_mask_compressstoreu_epi32(&m_indexes32[m_count], mask, m_vec);
          m_count += popcnt64_native(mask);
        }

        if (m > min_m)
        {
          // Filter out the last few square free m
          int count = int(m - min_m);
          __mmask16 load_mask = __mmask16((1u << count) - 1);
          __m512i m_vec = _mm512_sub_epi32(_mm512_set1_epi32(uint32_t(m)), m_offsets32);
          __m512i factor_vec = load_factor_tail_epi32_avx512(&factor_data[min_m + 1], load_mask, count, m_offsets32);
          __mmask16 mask = _mm512_cmpgt_epu32_mask(factor_vec, prime_vec);
          _mm512_mask_compressstoreu_epi32(&m_indexes32[m_count], mask, m_vec);
          m_count += popcnt64_native(mask);
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = m_indexes32[i];
          int64_t xpm = fast_div64(xp, factor.to_number(m));
          int64_t count = sieve.count_avx512(xpm - low);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m) * phi_xpm;
        }
      }
      else
      {
        __m512i prime_vec = _mm512_set1_epi64(prime);
        constexpr std::size_t max_m_count = m_indexes64.size() - 15;

        // AVX512: 8-lane 64-bit
        for (; m > min_m + 7; m -= 8)
        {
          if (m_count > max_m_count)
          {
            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = m_indexes64[i];
              int64_t xpm = fast_div64(xp, factor.to_number(m));
              int64_t count = sieve.count_avx512(xpm - low);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m) * phi_xpm;
            }
            m_count = 0;
          }

          // Filter out square free m values using AVX512
          // that satisfy: prime < factor.is_leaf(m)
          __m512i m_vec = _mm512_sub_epi64(_mm512_set1_epi64(m), m_offsets64);
          __m512i factor_vec = load_factor_epi64_avx512(&factor_data[m - 7], reverse64);
          __mmask8 mask = _mm512_cmpgt_epi64_mask(factor_vec, prime_vec);
          _mm512_mask_compressstoreu_epi64(&m_indexes64[m_count], mask, m_vec);
          m_count += popcnt64_native(mask);
        }

        if (m > min_m)
        {
          // Filter out the last few square free m
          int count = int(m - min_m);
          __mmask8 load_mask = __mmask8((1u << count) - 1);
          __m512i m_vec = _mm512_sub_epi64(_mm512_set1_epi64(m), m_offsets64);
          __m512i factor_vec = load_factor_tail_epi64_avx512(&factor_data[min_m + 1], load_mask, count, m_offsets64);
          __mmask8 mask = _mm512_cmpgt_epi64_mask(factor_vec, prime_vec);
          _mm512_mask_compressstoreu_epi64(&m_indexes64[m_count], mask, m_vec);
          m_count += popcnt64_native(mask);
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = m_indexes64[i];
          int64_t xpm = fast_div64(xp, factor.to_number(m));
          int64_t count = sieve.count_avx512(xpm - low);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m) * phi_xpm;
        }
      }

      phi[b] += sieve.get_total_count();
      sieve.cross_off_count(prime, b);
    }

    // For pi_sqrtz < b <= pi_x_star
    // Find all special leaves in the current segment
    // that are composed of 2 primes:
    // low <= x / (primes[b] * primes[l]) < high
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_low = min(fast_div(xp, low1), y);
      int64_t xp_high = min(fast_div(xp, high), y);
      int64_t min_m = max(xp_high, prime);
      int64_t max_m = min(fast_div(xp, prime * prime), xp_low);
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

#endif

#endif

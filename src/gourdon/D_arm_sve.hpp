///
/// @file  D_arm_sve.hpp
/// @brief ARM SVE implementation of the D formula (hard special
///        leaves) in Xavier Gourdon's prime counting algorithm. This
///        algorithm is identical to D_thread_default() in D.cpp
///        except that this algorithm has been partially vectorized
///        using ARM SVE.
///
///        For performance it is important that all ARM SVE helper
///        functions are inlined by the compiler. We achieve this
///        by annotating all ARM SVE helper functions using the same
///        ARM SVE __attribute__ and the ALWAYS_INLINE macro.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_ARM_SVE_HPP
#define D_ARM_SVE_HPP

#include <arm_sve.h>

namespace {

using namespace primecount;

#if defined(ENABLE_MULTIARCH_ARM_SVE)
__attribute__ ((target ("arch=armv8-a+sve")))
#endif
ALWAYS_INLINE svuint32_t load_factor_u32_arm_sve(svbool_t pg,
                                                 const uint16_t* factor_data)
{
  return svrev_u32(svld1uh_u32(pg, factor_data));
}

#if defined(ENABLE_MULTIARCH_ARM_SVE)
__attribute__ ((target ("arch=armv8-a+sve")))
#endif
ALWAYS_INLINE svuint32_t load_factor_u32_arm_sve(svbool_t pg,
                                                 const uint32_t* factor_data)
{
  return svrev_u32(svld1_u32(pg, factor_data));
}

#if defined(ENABLE_MULTIARCH_ARM_SVE)
__attribute__ ((target ("arch=armv8-a+sve")))
#endif
ALWAYS_INLINE svuint64_t load_factor_u64_arm_sve(svbool_t pg,
                                                 const uint16_t* factor_data)
{
  return svrev_u64(svld1uh_u64(pg, factor_data));
}

#if defined(ENABLE_MULTIARCH_ARM_SVE)
__attribute__ ((target ("arch=armv8-a+sve")))
#endif
ALWAYS_INLINE svuint64_t load_factor_u64_arm_sve(svbool_t pg,
                                                 const uint32_t* factor_data)
{
  return svrev_u64(svld1uw_u64(pg, factor_data));
}

template <typename T, typename Primes, typename FactorTableD>
#if defined(ENABLE_MULTIARCH_ARM_SVE)
__attribute__ ((target ("arch=armv8-a+sve")))
#endif
T D_thread_arm_sve(T x,
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
  T sum = 0;

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
        int64_t lanes32 = svcntw();
        ASSERT(svcntw() <= m_indexes32.size());
        std::size_t max_m_count = m_indexes32.size() - lanes32;
        svbool_t all32 = svptrue_b32();
        svuint32_t m_offsets32 = svindex_u32(0, 1);

        // ARM SVE 32-bit loop
        while (m > min_m)
        {
          // Filter out square free m values using ARM SVE
          // that satisfy: prime < factor.is_leaf(m)
          int64_t lane_count = min(m - min_m, lanes32);
          svbool_t load_pg = svwhilelt_b32(int64_t(0), lane_count);
          svbool_t store_pg = svrev_b32(load_pg);
          uint32_t base = uint32_t(m + lanes32 - lane_count);
          svuint32_t m_vec = svsub_u32_x(all32, svdup_n_u32(base), m_offsets32);
          svuint32_t factor_vec = load_factor_u32_arm_sve(load_pg, &factor_data[m + 1 - lane_count]);
          svbool_t mask = svcmpgt_n_u32(store_pg, factor_vec, uint32_t(prime));
          int64_t matches = svcntp_b32(store_pg, mask);
          svuint32_t compact = svcompact_u32(mask, m_vec);
          svbool_t compact_pg = svwhilelt_b32(int64_t(0), matches);
          svst1_u32(compact_pg, &m_indexes32[m_count], compact);

          m_count += matches;
          m -= lane_count;

          if (m_count > max_m_count)
          {
            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = m_indexes32[i];
              int64_t xpm = fast_div64(xp, factor.to_number(m));
              int64_t count = sieve.count_arm_sve(xpm - low);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m) * phi_xpm;
            }
            m_count = 0;
          }
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = m_indexes32[i];
          int64_t xpm = fast_div64(xp, factor.to_number(m));
          int64_t count = sieve.count_arm_sve(xpm - low);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m) * phi_xpm;
        }
      }
      else
      {
        int64_t lanes64 = svcntd();
        ASSERT(svcntd() <= m_indexes64.size());
        std::size_t max_m_count = m_indexes64.size() - lanes64;
        svbool_t all64 = svptrue_b64();
        svuint64_t m_offsets64 = svindex_u64(0, 1);

        // ARM SVE 64-bit loop
        while (m > min_m)
        {
          // Filter out square free m values using ARM SVE
          // that satisfy: prime < factor.is_leaf(m)
          int64_t lane_count = min(m - min_m, lanes64);
          svbool_t load_pg = svwhilelt_b64(int64_t(0), lane_count);
          svbool_t store_pg = svrev_b64(load_pg);
          uint64_t base = uint64_t(m + lanes64 - lane_count);
          svuint64_t m_vec = svsub_u64_x(all64, svdup_n_u64(base), m_offsets64);
          svuint64_t factor_vec = load_factor_u64_arm_sve(load_pg, &factor_data[m + 1 - lane_count]);
          svbool_t mask = svcmpgt_n_u64(store_pg, factor_vec, uint64_t(prime));
          int64_t matches = svcntp_b64(store_pg, mask);
          svint64_t compact = svreinterpret_s64_u64(svcompact_u64(mask, m_vec));
          svbool_t compact_pg = svwhilelt_b64(int64_t(0), matches);
          svst1_s64(compact_pg, &m_indexes64[m_count], compact);

          m_count += matches;
          m -= lane_count;

          if (m_count > max_m_count)
          {
            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = m_indexes64[i];
              int64_t xpm = fast_div64(xp, factor.to_number(m));
              int64_t count = sieve.count_arm_sve(xpm - low);
              int64_t phi_xpm = phi[b] + count;
              sum -= factor.mu(m) * phi_xpm;
            }
            m_count = 0;
          }
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = m_indexes64[i];
          int64_t xpm = fast_div64(xp, factor.to_number(m));
          int64_t count = sieve.count_arm_sve(xpm - low);
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
        int64_t count = sieve.count_arm_sve(xpq - low);
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

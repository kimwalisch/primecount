///
/// @file  D_default.hpp
/// @brief This is a highly optimized implementation of the D(x, y)
///        formula in Xavier Gourdon's prime counting algorithm. The D
///        formula is very similar to the formula of the hard special
///        leaves in the Deleglise-Rivat algorithm. Hence this
///        algorithm is very similar to S2_hard.cpp, except that in
///        this implementation the square free leaves have been more
///        heavily optimized (branchfree + CPU pipelining).
///
///        This implementation uses multi-threading with advanced load
///        balancing, it scales well up to a large number of CPU cores
///        because the compute threads are completely independent from
///        each other. This implementation also uses the highly
///        optimized Sieve class and the FactorTableD class which is a
///        compressed lookup table of moebius function values,
///        least prime factors and max prime factors.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves.pdf
///        https://github.com/kimwalisch/primecount/blob/master/doc/Hard-Special-Leaves-SIMD-Filtering.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef D_DEFAULT_HPP
#define D_DEFAULT_HPP

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

namespace {

using namespace primecount;

/// Compute the contribution of the hard special leaves using
/// a segmented sieve. Each thread processes the interval
/// [low, low + segment_size * segments[.
///
template <typename T, typename Primes, typename FactorTable>
T D_thread_default(T x,
                   int64_t x_star,
                   int64_t xz,
                   int64_t y,
                   int64_t z,
                   int64_t k,
                   const Primes& primes,
                   const PiTable& pi,
                   const FactorTable& factor,
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

      min_m = FactorTable::to_index(min_m);
      max_m = FactorTable::to_index(max_m);
      int64_t encoded_prime = FactorTable::encode(b);
      int64_t m = max_m;
      std::size_t m_count = 0;

      // 32-bit code path
      if (max_m <= UINT32_MAX ||
          sizeof(T) <= sizeof(uint64_t))
      {
        constexpr std::size_t max_m_count = m_indexes32.size() - 4;

        // Filter out square free m values branchlessly
        // that satisfy: factor_table[m] > encoded_prime
        for (; m >= min_m + 4; m -= 4)
        {
          m_indexes32[m_count] = uint32_t(m);
          m_count += (factor_table[m] > encoded_prime);
          m_indexes32[m_count] = uint32_t(m - 1);
          m_count += (factor_table[m - 1] > encoded_prime);
          m_indexes32[m_count] = uint32_t(m - 2);
          m_count += (factor_table[m - 2] > encoded_prime);
          m_indexes32[m_count] = uint32_t(m - 3);
          m_count += (factor_table[m - 3] > encoded_prime);

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
      else // 64-bit code path
      {
        constexpr std::size_t max_m_count = m_indexes64.size() - 4;

        // Filter out square free m values branchlessly
        // that satisfy: factor_table[m] > encoded_prime
        for (; m >= min_m + 4; m -= 4)
        {
          m_indexes64[m_count] = m;
          m_count += (factor_table[m] > encoded_prime);
          m_indexes64[m_count] = m - 1;
          m_count += (factor_table[m - 1] > encoded_prime);
          m_indexes64[m_count] = m - 2;
          m_count += (factor_table[m - 2] > encoded_prime);
          m_indexes64[m_count] = m - 3;
          m_count += (factor_table[m - 3] > encoded_prime);

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

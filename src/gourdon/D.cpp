///
/// @file  D.cpp
/// @brief This is a highly optimized implementation of the D(x, y)
///        formula in Xavier Gourdon's prime counting algorithm. The D
///        formula is very similar to the formula of the hard special
///        leaves in the Deleglise-Rivat algorithm. Hence this
///        algorithm is very similar to S2_hard.cpp, expect that in
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
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "FactorTableD.hpp"

#include <primecount-internal.hpp>
#include <PiTable.hpp>
#include <Sieve.hpp>
#include <LoadBalancerS2.hpp>
#include <fast_div.hpp>
#include <generate_primes.hpp>
#include <phi_vector.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>

#include <stdint.h>
#include <utility>

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <cpu_supports_avx512_vpopcnt.hpp>
#endif

#if defined(ENABLE_ARM_SVE) || \
    defined(ENABLE_MULTIARCH_ARM_SVE)
  #include "D_arm_sve.hpp"
#elif defined(ENABLE_AVX512_VPOPCNT) || \
      defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include "D_avx512.hpp"
#endif

namespace {

using namespace primecount;

#if !defined(ENABLE_AVX512_VPOPCNT) && \
    !defined(ENABLE_ARM_SVE)

/// Compute the contribution of the hard special leaves using
/// a segmented sieve. Each thread processes the interval
/// [low, low + segment_size * segments[.
///
template <typename T, typename Primes, typename FactorTableD>
T D_thread_default(T x,
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

  Array<uint32_t, 128> m_indexes32;
  Array< int64_t, 128> m_indexes64;
  Array< int64_t, 128> xpm_cache;
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

      // 32-bit code path
      if (FactorTableD::max() <= UINT32_MAX ||
          max_m <= UINT32_MAX)
      {
        constexpr std::size_t max_m_count = m_indexes32.size() - 4;

        // Filter out square free m values branchlessly
        // that satisfy: prime < factor.is_leaf(m)
        for (; m > min_m + 3; m -= 4)
        {
          m_indexes32[m_count] = uint32_t(m);
          m_count += (prime < factor_data[m]);
          m_indexes32[m_count] = uint32_t(m - 1);
          m_count += (prime < factor_data[m - 1]);
          m_indexes32[m_count] = uint32_t(m - 2);
          m_count += (prime < factor_data[m - 2]);
          m_indexes32[m_count] = uint32_t(m - 3);
          m_count += (prime < factor_data[m - 3]);

          if (m_count > max_m_count)
          {
            // Batch calculate xp/m to improve CPU pipelining
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = factor.to_number(m_indexes32[i]);
              xpm_cache[i] = fast_div64(xp, m);
            }

            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t xpm = xpm_cache[i];
              int64_t count = sieve.count(xpm - low);
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
          m_count += (prime < factor_data[m]);
        }

        // Batch calculate xp/m to improve CPU pipelining
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = factor.to_number(m_indexes32[i]);
          xpm_cache[i] = fast_div64(xp, m);
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t xpm = xpm_cache[i];
          int64_t count = sieve.count(xpm - low);
          int64_t phi_xpm = phi[b] + count;
          sum -= factor.mu(m_indexes32[i]) * phi_xpm;
        }
      }
      else // 64-bit code path
      {
        constexpr std::size_t max_m_count = m_indexes64.size() - 4;

        // Filter out square free m values branchlessly
        // that satisfy: prime < factor.is_leaf(m)
        for (; m > min_m + 3; m -= 4)
        {
          m_indexes64[m_count] = m;
          m_count += (prime < factor_data[m]);
          m_indexes64[m_count] = m - 1;
          m_count += (prime < factor_data[m - 1]);
          m_indexes64[m_count] = m - 2;
          m_count += (prime < factor_data[m - 2]);
          m_indexes64[m_count] = m - 3;
          m_count += (prime < factor_data[m - 3]);

          if (m_count > max_m_count)
          {
            // Batch calculate xp/m to improve CPU pipelining
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t m = factor.to_number(m_indexes64[i]);
              xpm_cache[i] = fast_div64(xp, m);
            }

            // Process the next few special leaves that are
            // composed of a prime and a square free number:
            // low <= x / (primes[b] * m) < high
            for (std::size_t i = 0; i < m_count; i++)
            {
              int64_t xpm = xpm_cache[i];
              int64_t count = sieve.count(xpm - low);
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
          m_count += (prime < factor_data[m]);
        }

        // Batch calculate xp/m to improve CPU pipelining
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t m = factor.to_number(m_indexes64[i]);
          xpm_cache[i] = fast_div64(xp, m);
        }

        // Process the last few m values
        for (std::size_t i = 0; i < m_count; i++)
        {
          int64_t xpm = xpm_cache[i];
          int64_t count = sieve.count(xpm - low);
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

#endif

/// Runtime dispatch to highly optimized SIMD algorithm if the CPU
/// supports the required instruction set.
///
template <typename T, typename... Args>
T D_thread(Args&&... args)
{
  #if defined(ENABLE_AVX512_VPOPCNT)
    return D_thread_avx512<T>(std::forward<Args>(args)...);
  #elif defined(ENABLE_ARM_SVE)
    return D_thread_arm_sve<T>(std::forward<Args>(args)...);
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    return cpu_supports_avx512_vpopcnt
      ? D_thread_avx512<T>(std::forward<Args>(args)...)
      : D_thread_default<T>(std::forward<Args>(args)...);
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    return cpu_supports_sve
      ? D_thread_arm_sve<T>(std::forward<Args>(args)...)
      : D_thread_default<T>(std::forward<Args>(args)...);
  #else
    return D_thread_default<T>(std::forward<Args>(args)...);
  #endif
}

string_view_t D_algo_name()
{
  #if defined(ENABLE_AVX512_VPOPCNT)
    return "Algorithm: AVX512";
  #elif defined(ENABLE_ARM_SVE)
    return "Algorithm: ARM SVE";
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    return cpu_supports_avx512_vpopcnt
      ? "Algorithm: AVX512"
      : "Algorithm: POPCNT64";
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    return cpu_supports_sve
      ? "Algorithm: ARM SVE"
      : "Algorithm: POPCNT64";
  #else
    return "Algorithm: POPCNT64";
  #endif
}

/// Calculate the contribution of the hard special leaves.
///
/// This is a parallel D(x, y) implementation with advanced load
/// balancing. As most special leaves tend to be in the first segments
/// we start off with a tiny segment size and one segment per thread.
/// After each iteration we dynamically increase the segment size (until
/// it reaches some limit) or the number of segments.
///
/// D(x, y) has been parallelized using an idea devised by Xavier
/// Gourdon. The idea is to make the individual threads completely
/// independent from each other so that no thread depends on values
/// calculated by another thread. The benefit of this approach is that
/// the algorithm will scale well up to a very large number of CPU
/// cores. In order to make the threads independent from each other
/// each thread needs to precompute a lookup table of phi(x, a) values
/// (this is done in D_thread(x, y)) every time the thread starts
/// a new computation.
///
template <typename T, typename Primes, typename FactorTableD>
T D_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t k,
           T d_approx,
           const Primes& primes,
           const FactorTableD& factor,
           int threads,
           bool is_print)
{
  int64_t xz = x / z;
  int64_t x_star = get_x_star_gourdon(x, y);

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  int64_t thread_threshold = 1 << 20;
  int max_threads = (int) std::pow(xz, 1 / 3.7);
  threads = std::min(threads, max_threads);
  threads = ideal_num_threads(xz, threads, thread_threshold);
  LoadBalancerS2 loadBalancer(x, xz, d_approx, threads, is_print);
  PiTable pi(y, threads);

  #pragma omp parallel num_threads(threads)
  {
    ThreadData thread;

    while (loadBalancer.get_work(thread))
    {
      // Unsigned integer division is usually slightly
      // faster than signed integer division
      using UT = typename pstd::make_unsigned<T>::type;

      thread.start_time();
      UT sum = D_thread<UT>(x, x_star, xz, y, z, k, primes, pi, factor, thread);
      thread.sum = (T) sum;
      thread.stop_time();
    }
  }

  T sum = (T) loadBalancer.get_sum();

  return sum;
}

} // namespace

namespace primecount {

int64_t D(int64_t x,
          int64_t y,
          int64_t z,
          int64_t k,
          int64_t d_approx,
          int threads,
          bool is_print)
{
  double time;

  if (is_print)
  {
    print("");
    print("=== D(x, y) ===");
    print(D_algo_name());
    print_gourdon_vars(x, y, z, k, threads);
    time = get_time();
  }

  FactorTableD<uint16_t> factor(y, z, threads);
  auto primes = generate_primes<uint32_t>(y);
  int64_t sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, is_print);

  if (is_print)
    print("D", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t D(int128_t x,
                   int64_t y,
                   int64_t z,
                   int64_t k,
                   int128_t d_approx,
                   int threads,
                   bool is_print)
{
  double time;

  if (is_print)
  {
    print("");
    print("=== D(x, y) ===");
    print(D_algo_name());
    print_gourdon_vars(x, y, z, k, threads);
    time = get_time();
  }

  int128_t sum;

  // uses less memory
  if (z <= FactorTableD<uint16_t>::max())
  {
    FactorTableD<uint16_t> factor(y, z, threads);
    auto primes = generate_primes<uint32_t>(y);
    sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, is_print);
  }
  else
  {
    FactorTableD<uint32_t> factor(y, z, threads);
    auto primes = generate_primes<int64_t>(y);
    sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, is_print);
  }

  if (is_print)
    print("D", sum, time);

  return sum;
}

#endif

} // namespace

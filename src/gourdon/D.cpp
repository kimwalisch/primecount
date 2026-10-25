///
/// @file  D.cpp
/// @brief Implementation of the D formula (hard special leaves)
///        in Xavier Gourdon's prime counting algorithm.
///
///        This file handles runtime dispatch to optimized SIMD
///        implementations, thread scheduling, and load balancing.
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
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "DHyperbola.hpp"

#include <primecount-internal.hpp>
#include <cpu_arch_macros.hpp>
#include <FactorTable.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <LoadBalancerS2.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <print.hpp>
#include <PiTable.hpp>
#include <RelaxedAtomic.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <utility>

#if defined(ENABLE_AVX512_VPOPCNT)
  #define D_DEFAULT_HEADER "D_avx512.hpp"
  #define D_DEFAULT_ALGO_NAME "Algorithm: AVX512"
  #define D_THREAD_DEFAULT D_thread_avx512
#elif defined(ENABLE_ARM_SVE)
  #define D_DEFAULT_HEADER "D_arm_sve.hpp"
  #define D_DEFAULT_ALGO_NAME "Algorithm: ARM SVE"
  #define D_THREAD_DEFAULT D_thread_arm_sve
#elif defined(ENABLE_AVX2)
  #define D_DEFAULT_HEADER "D_avx2.hpp"
  #define D_DEFAULT_ALGO_NAME "Algorithm: AVX2"
  #define D_THREAD_DEFAULT D_thread_avx2
#elif defined(ENABLE_ARM_NEON)
  #define D_DEFAULT_HEADER "D_arm_neon.hpp"
  #define D_DEFAULT_ALGO_NAME "Algorithm: ARM NEON"
  #define D_THREAD_DEFAULT D_thread_arm_neon
#else
  #define D_DEFAULT_HEADER "D_default.hpp"
  #define D_DEFAULT_ALGO_NAME "Algorithm: POPCNT64"
  #define D_THREAD_DEFAULT D_thread_default
#endif

#if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include <cpu_supports_avx512_vpopcnt.hpp>
  #include "D_avx512.hpp"
#endif

#if defined(ENABLE_MULTIARCH_AVX2)
  #include <cpu_supports_avx2.hpp>
  #include "D_avx2.hpp"
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
  #include "D_arm_sve.hpp"
#endif

// Fastest algorithm supported by the compiler flags
#include D_DEFAULT_HEADER

namespace {

using namespace primecount;

/// Runtime dispatch to highly optimized SIMD algorithm if the CPU
/// supports the required instruction set.
///
template <typename T, typename... Args>
T D_thread(Args&&... args)
{
  // Unsigned integer division is usually
  // faster than signed integer division.
  using UT = typename pstd::make_unsigned<T>::type;

  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    if (cpu_supports_avx512_vpopcnt)
      return D_thread_avx512<UT>(std::forward<Args>(args)...);
  #endif

  #if defined(ENABLE_MULTIARCH_AVX2)
    if (cpu_supports_avx2)
      return D_thread_avx2<UT>(std::forward<Args>(args)...);
  #endif

  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      return D_thread_arm_sve<UT>(std::forward<Args>(args)...);
  #endif

  return D_THREAD_DEFAULT<UT>(std::forward<Args>(args)...);
}

string_view_t D_algo_name()
{
  #if defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    if (cpu_supports_avx512_vpopcnt)
      return "Algorithm: AVX512";
  #endif

  #if defined(ENABLE_MULTIARCH_AVX2)
    if (cpu_supports_avx2)
      return "Algorithm: AVX2";
  #endif

  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      return "Algorithm: ARM SVE";
  #endif

  return D_DEFAULT_ALGO_NAME;
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
template <typename T, typename Primes, typename FactorTable>
T D_OpenMP(T x,
           int64_t y,
           int64_t k,
           const PiTable& pi,
           const Primes& primes,
           const FactorTable& factor,
           int threads,
           bool is_print)
{
  using UT = typename pstd::make_unsigned<T>::type;

  // Hard special leaves satisfy m > p and
  // p^3 * m <= x, hence p^4 < x.
  int64_t x14 = iroot<4>(x);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t xy = x / y;
  int64_t pi_root3_xy = pi[iroot<3>(xy)];

  // Two-factor Moebius prefixes, adjusted for triple products below.
  // Only entries with primes[b]^2 <= y < primes[b]^4 are used.
  Vector<int64_t> mertens_y(pi_sqrty + 1);
  int64_t mertens = 1 - pi[y] + pi_sqrty;
  for (int64_t b = pi_sqrty; b > k; b--)
  {
    mertens_y[b] = mertens;
    mertens += pi[y / primes[b]] - b - 1;
  }

  // For three-prime leaves the boundary prefix is Phi(y) - Q(y),
  // which adds all triple products to the two-factor Moebius sum.
  int64_t max_triple_prime = iroot<3>(y);
  int64_t min_triple_prime = max(iroot<4>(y), xy / y);
  int64_t min_triple_b = max(k, pi[min(min_triple_prime, max_triple_prime)]);
  int64_t triples = 0;

  for (int64_t b = pi[max_triple_prime]; b > min_triple_b; b--)
  {
    int64_t yp = y / primes[b];
    for (int64_t i = b, last = pi[isqrt(yp)]; i <= last; i++)
      triples += pi[yp / primes[i]] - i + 1;
    mertens_y[b] += triples;
  }

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  int64_t thread_threshold = 1 << 20;
  int max_threads = (int) std::pow(xy, 1 / 3.7);
  threads = std::min(threads, max_threads);
  threads = ideal_num_threads(xy, threads, thread_threshold);
  INDETERMINATE LoadBalancerS2 loadBalancer(x, y, xy, threads, is_print);
  INDETERMINATE RelaxedAtomic<int64_t> min_b(max(k, pi_root3_xy) + 1);
  UT sum = 0;

  #pragma omp parallel num_threads(threads) reduction(+: sum)
  {
    ThreadData thread;

    while (loadBalancer.get_work(thread))
    {
      thread.start_time = get_time();
      thread.sum = D_thread<T>(x, y, xy, k, x14, primes, pi, factor, mertens_y, thread);
      thread.stop_time = get_time();
      sum += thread.sum;
    }

    // Compute hyperbola terms with positions < prime^2 once per prime.
    // These terms are skipped by the segmented kernels.
    for (int64_t b = min_b++; b <= pi_sqrty; b = min_b++)
    {
      UT xp = UT(x / primes[b]);
      sum += easy_hyperbola_sum(xp, primes[b], b, y, mertens_y[b], primes, pi);
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t D(int64_t x,
          int64_t y,
          int64_t k,
          const PiTable& pi,
          int threads,
          bool is_print)
{
  ASSERT(uint64_t(y) < pi.size());
  double time;

  if (is_print)
  {
    print("");
    print("=== D(x, y) ===");
    print(D_algo_name());
    print_gourdon_vars(x, y, k, threads);
    time = get_time();
  }

  FactorTable<uint16_t> factor(y, threads);
  auto primes = pi.get_primes<uint32_t>(y, threads);
  int64_t sum = D_OpenMP(x, y, k, pi, primes, factor, threads, is_print);

  if (is_print)
    print("D", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t D(int128_t x,
           int64_t y,
           int64_t k,
           const PiTable& pi,
           int threads,
           bool is_print)
{
  ASSERT(uint64_t(y) < pi.size());
  double time;

  if (is_print)
  {
    print("");
    print("=== D(x, y) ===");
    print(D_algo_name());
    print_gourdon_vars(x, y, k, threads);
    time = get_time();
  }

  int128_t sum;

  // Use 16-bit factor table entries whenever possible.
  if (y <= FactorTable<uint16_t>::max())
  {
    FactorTable<uint16_t> factor(y, threads);

    if (y <= UINT32_MAX)
    {
      auto primes = pi.get_primes<uint32_t>(y, threads);
      sum = D_OpenMP(x, y, k, pi, primes, factor, threads, is_print);
    }
    else
    {
      auto primes = pi.get_primes<int64_t>(y, threads);
      sum = D_OpenMP(x, y, k, pi, primes, factor, threads, is_print);
    }
  }
  else
  {
    FactorTable<uint32_t> factor(y, threads);
    auto primes = pi.get_primes<int64_t>(y, threads);
    sum = D_OpenMP(x, y, k, pi, primes, factor, threads, is_print);
  }

  if (is_print)
    print("D", sum, time);

  return sum;
}

#endif

} // namespace

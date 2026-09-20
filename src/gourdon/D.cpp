///
/// @file  D.cpp
/// @brief Implementation of the D formula (hard special leaves)
///        in Xavier Gourdon's prime counting algorithm.
///
///        This file handles runtime dispatch to optimized SIMD
///        implementations, thread scheduling, and load balancing.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "FactorTableD.hpp"

#include <primecount-internal.hpp>
#include <cpu_arch_macros.hpp>
#include <macros.hpp>
#include <PiTable.hpp>
#include <sieve/Sieve.hpp>
#include <LoadBalancerS2.hpp>
#include <fast_div.hpp>
#include <phi_vector.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>

#include <stdint.h>
#include <utility>

#if defined(ENABLE_ARM_SVE)
  #include "D_arm_sve.hpp"
#elif defined(ENABLE_AVX512_VPOPCNT)
  #include "D_avx512.hpp"
#elif defined(ENABLE_MULTIARCH_ARM_SVE)
  #include "D_arm_sve.hpp"
  #include <cpu_supports_arm_sve.hpp>
#elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
  #include "D_avx512.hpp"
  #include <cpu_supports_avx512_vpopcnt.hpp>
#endif

#if !defined(ENABLE_ARM_SVE) && \
    !defined(ENABLE_AVX512_VPOPCNT)
  #include "D_default.hpp"
#endif

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

  #if defined(ENABLE_AVX512_VPOPCNT)
    return D_thread_avx512<UT>(std::forward<Args>(args)...);
  #elif defined(ENABLE_ARM_SVE)
    return D_thread_arm_sve<UT>(std::forward<Args>(args)...);
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    return cpu_supports_avx512_vpopcnt
      ? D_thread_avx512<UT>(std::forward<Args>(args)...)
      : D_thread_default<UT>(std::forward<Args>(args)...);
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    return cpu_supports_sve
      ? D_thread_arm_sve<UT>(std::forward<Args>(args)...)
      : D_thread_default<UT>(std::forward<Args>(args)...);
  #else
    return D_thread_default<UT>(std::forward<Args>(args)...);
  #endif
}

string_view_t D_algo_name()
{
  #if defined(ENABLE_ARM_NEON)
    #define DEFAULT_ALGO_NAME "Algorithm: ARM NEON"
  #else
    #define DEFAULT_ALGO_NAME "Algorithm: POPCNT64"
  #endif

  #if defined(ENABLE_AVX512_VPOPCNT)
    return "Algorithm: AVX512";
  #elif defined(ENABLE_ARM_SVE)
    return "Algorithm: ARM SVE";
  #elif defined(ENABLE_MULTIARCH_AVX512_VPOPCNT)
    return cpu_supports_avx512_vpopcnt
      ? "Algorithm: AVX512"
      : DEFAULT_ALGO_NAME;
  #elif defined(ENABLE_MULTIARCH_ARM_SVE)
    return cpu_supports_sve
      ? "Algorithm: ARM SVE"
      : DEFAULT_ALGO_NAME;
  #else
    return DEFAULT_ALGO_NAME;
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
template <typename T, typename Primes, typename FactorTable>
T D_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t k,
           const PiTable& pi,
           const Primes& primes,
           const FactorTable& factor,
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
  INDETERMINATE LoadBalancerS2 loadBalancer(x, y, xz, threads, is_print);
  T sum = 0;

  #pragma omp parallel num_threads(threads) reduction(+: sum)
  {
    ThreadData thread;

    while (loadBalancer.get_work(thread))
    {
      thread.start_time = get_time();
      thread.sum = D_thread<T>(x, x_star, xz, y, z, k, primes, pi, factor, thread);
      thread.stop_time = get_time();
      sum += thread.sum;
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t D(int64_t x,
          int64_t y,
          int64_t z,
          int64_t k,
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
  PiTable pi(y, threads);
  auto primes = pi.get_primes<uint32_t>(y, threads);
  int64_t sum = D_OpenMP(x, y, z, k, pi, primes, factor, threads, is_print);

  if (is_print)
    print("D", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t D(int128_t x,
           int64_t y,
           int64_t z,
           int64_t k,
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

  // Use 16-bit factor table entries whenever possible.
  if (z <= FactorTableD<uint16_t>::max())
  {
    FactorTableD<uint16_t> factor(y, z, threads);
    PiTable pi(y, threads);

    if (y <= UINT32_MAX)
    {
      auto primes = pi.get_primes<uint32_t>(y, threads);
      sum = D_OpenMP(x, y, z, k, pi, primes, factor, threads, is_print);
    }
    else
    {
      auto primes = pi.get_primes<int64_t>(y, threads);
      sum = D_OpenMP(x, y, z, k, pi, primes, factor, threads, is_print);
    }
  }
  else
  {
    FactorTableD<uint32_t> factor(y, z, threads);
    PiTable pi(y, threads);
    auto primes = pi.get_primes<int64_t>(y, threads);
    sum = D_OpenMP(x, y, z, k, pi, primes, factor, threads, is_print);
  }

  if (is_print)
    print("D", sum, time);

  return sum;
}

#endif

} // namespace

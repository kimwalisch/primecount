///
/// @file  AC.cpp
/// @brief Setup and dispatch for the A + C formulas (easy special
///        leaves) in Xavier Gourdon's prime counting algorithm.
///
///        This file prepares the primes, dispatches to the scalar or
///        SIMD implementation, and prints algorithm details and the
///        result when requested.
///
///        primecount's C formula differs from Gourdon's original
///        C formula: the easy special leaves with
///        (x/y)^(1/3) < p <= sqrt(y) and x/p^3 < m <= y
///        are computed as part of the D formula using the hyperbola
///        method. Here p is a prime and m is a prime or a product
///        of 2 distinct primes, all greater than p.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <cpu_arch_macros.hpp>
#include <gourdon.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>
#include <PiTable.hpp>

#include <stdint.h>
#include <utility>

#if defined(ENABLE_ARM_SVE)
  #define AC_DEFAULT_HEADER "AC_arm_sve.hpp"
  #define AC_DEFAULT_ALGO_NAME "Algorithm: ARM SVE"
  #define AC_OPENMP_DEFAULT AC_OpenMP_arm_sve
#elif defined(ENABLE_LIBDIVIDE)
  #define AC_DEFAULT_HEADER "AC_libdivide.hpp"
  #define AC_DEFAULT_ALGO_NAME "Algorithm: libdivide"
  #define AC_OPENMP_DEFAULT AC_OpenMP_libdivide
#else
  #define AC_DEFAULT_HEADER "AC_default.hpp"
  #define AC_DEFAULT_ALGO_NAME "Algorithm: CPU div"
  #define AC_OPENMP_DEFAULT AC_OpenMP_default
#endif

#if defined(ENABLE_MULTIARCH_ARM_SVE)
  #include <cpu_supports_arm_sve.hpp>
  #include "AC_arm_sve.hpp"
#endif

// Fastest algorithm supported by the compiler flags
#include AC_DEFAULT_HEADER

namespace {

/// Runtime dispatch to highly optimized SIMD algorithm if
/// the CPU supports the required instruction set.
///
template <typename T, typename... Args>
T AC_OpenMP(T x, Args&&... args)
{
  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      return AC_OpenMP_arm_sve(x, std::forward<Args>(args)...);
  #endif

  return AC_OPENMP_DEFAULT(x, std::forward<Args>(args)...);
}

string_view_t AC_algo_name()
{
  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    if (cpu_supports_sve)
      return "Algorithm: ARM SVE";
  #endif

  return AC_DEFAULT_ALGO_NAME;
}

} // namespace

namespace primecount {

int64_t AC(int64_t x,
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
    print("=== AC(x, y) ===");
    print(AC_algo_name());
    print_gourdon_vars(x, y, k, threads);
    time = get_time();
  }

  // The A and C2 formulas use the primes p > sqrt(y) and
  // q <= sqrt(x / p), hence all primes are <= sqrt(x / sqrt(y)).
  int64_t max_prime = min(isqrt(x / isqrt(y)), y);
  int64_t x_star = get_x_star_gourdon(x, y);
  auto primes = pi.get_primes<uint32_t>(max_prime, threads);
  int64_t sum = AC_OpenMP((uint64_t) x, y, k, x_star, pi, primes, threads, is_print);

  if (is_print)
    print("A + C", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t AC(int128_t x,
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
    print("=== AC(x, y) ===");
    print(AC_algo_name());
    print_gourdon_vars(x, y, k, threads);
    time = get_time();
  }

  // The A and C2 formulas use the primes p > sqrt(y) and
  // q <= sqrt(x / p), hence all primes are <= sqrt(x / sqrt(y)).
  int64_t max_prime = min(isqrt(x / isqrt(y)), y);
  int64_t x_star = get_x_star_gourdon(x, y);
  int128_t sum;

  // uses less memory
  if (max_prime <= pstd::numeric_limits<uint32_t>::max())
  {
    auto primes = pi.get_primes<uint32_t>(max_prime, threads);
    sum = AC_OpenMP((uint128_t) x, y, k, x_star, pi, primes, threads, is_print);
  }
  else
  {
    auto primes = pi.get_primes<int64_t>(max_prime, threads);
    sum = AC_OpenMP((uint128_t) x, y, k, x_star, pi, primes, threads, is_print);
  }

  if (is_print)
    print("A + C", sum, time);

  return sum;
}

#endif

} // namespace

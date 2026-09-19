///
/// @file  AC.cpp
/// @brief Implementation of the A + C formulas in Xavier Gourdon's
///        prime counting algorithm. In this implementation the memory
///        usage of the pi[x] lookup table has been reduced from
///        O(x^(1/2)) to O(x^(1/4)) by using a segmented pi[x] lookup
///        table. In each segment we process the leaves that satisfy:
///        low <= x / (prime * m) < high.
///
///        The A & C formulas roughly correspond to the easy special
///        leaves in the Deleglise-Rivat algorithm. Since both
///        formulas use a very similar segmented algorithm that goes
///        up to x^(1/2) it makes sense to merge the A & C formulas
///        hence reducing the runtime complexity by a factor of
///        O(x^(1/2) * ln ln x^(1/2)) and avoiding initializing some
///        data structures twice. Merging the A & C formulas also
///        improves scaling on systems with many CPU cores.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "LoadBalancerAC.hpp"
#include "SegmentedPiTable.hpp"

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <macros.hpp>
#include <cpu_arch_macros.hpp>
#include <fast_div.hpp>
#include <gourdon.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <utility>

#if defined(ENABLE_ARM_SVE)
  #include "AC_arm_sve.hpp"
#else
  #if defined(ENABLE_MULTIARCH_ARM_SVE)
    #include "AC_arm_sve.hpp"
    #include <cpu_supports_arm_sve.hpp>
  #endif
  // Portable fallback AC algorithms
  #if defined(ENABLE_LIBDIVIDE)
    #include "AC_libdivide.hpp"
  #else
    #include "AC_default.hpp"
  #endif
#endif

namespace {

/// Runtime dispatch to highly optimized SIMD algorithm if
/// the CPU supports the required instruction set.
///
template <typename T, typename... Args>
T AC_OpenMP(T x, Args&&... args)
{
  #if defined(ENABLE_ARM_SVE)
    return AC_OpenMP_arm_sve(x, std::forward<Args>(args)...);
  #else
    #if defined(ENABLE_MULTIARCH_ARM_SVE)
      if (cpu_supports_sve)
        return AC_OpenMP_arm_sve(x, std::forward<Args>(args)...);
    #endif
    #if defined(ENABLE_LIBDIVIDE)
      return AC_OpenMP_libdivide(x, std::forward<Args>(args)...);
    #else
      return AC_OpenMP_default(x, std::forward<Args>(args)...);
    #endif
  #endif
}

string_view_t AC_algo_name()
{
  #if defined(ENABLE_ARM_SVE)
    return "Algorithm: ARM SVE";
  #else
    #if defined(ENABLE_MULTIARCH_ARM_SVE)
      if (cpu_supports_sve)
        return "Algorithm: ARM SVE";
    #endif
    #if defined(ENABLE_LIBDIVIDE)
      return "Algorithm: libdivide";
    #else
      return "Algorithm: CPU div";
    #endif
  #endif
}

} // namespace

namespace primecount {

int64_t AC(int64_t x,
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
    print("=== AC(x, y) ===");
    print(AC_algo_name());
    print_gourdon_vars(x, y, z, k, threads);
    time = get_time();
  }

  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);

  // The A and C algorithms use the large PiTable only
  // for initialization. The inner-most loops of those
  // algorithms use the small SegmentedPiTable instead
  // which fits into the CPU's cache.
  PiTable pi(max_prime, threads);

  auto primes = pi.get_primes<uint32_t>(max_prime, threads);
  int64_t sum = AC_OpenMP((uint64_t) x, y, z, k, x_star, pi, primes, threads, is_print);

  if (is_print)
    print("A + C", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t AC(int128_t x,
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
    print("=== AC(x, y) ===");
    print(AC_algo_name());
    print_gourdon_vars(x, y, z, k, threads);
    time = get_time();
  }

  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);

  // The A and C algorithms use the large PiTable only
  // for initialization. The inner-most loops of those
  // algorithms use the small SegmentedPiTable instead
  // which fits into the CPU's cache.
  PiTable pi(max_prime, threads);
  int128_t sum;

  // uses less memory
  if (max_prime <= pstd::numeric_limits<uint32_t>::max())
  {
    auto primes = pi.get_primes<uint32_t>(max_prime, threads);
    sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, pi, primes, threads, is_print);
  }
  else
  {
    auto primes = pi.get_primes<int64_t>(max_prime, threads);
    sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, pi, primes, threads, is_print);
  }

  if (is_print)
    print("A + C", sum, time);

  return sum;
}

#endif

} // namespace

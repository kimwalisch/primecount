///
/// @file  pi_gourdon.cpp
/// @brief Implementation of Xavier Gourdon's prime counting
///        algorithm. Xavier Gourdon's algorithm is an improved
///        version of the Deleglise-Rivat algorithm, according to my
///        benchmarks Gourdon's algorithm runs up to 2x faster than
///        the Deleglise-Rivat algorithm.
///
///        Xavier Gourdon formula:
///        pi(x) = A - B + C + D + Phi0 + Sigma
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <PhiTiny.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <string>

namespace primecount {

/// Calculate the number of primes below x using
/// Xavier Gourdon's algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_gourdon_64(int64_t x,
                      int threads,
                      bool is_print)
{
  if (x < 2)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print)
  {
    print("");
    print("=== pi_gourdon_64(x) ===");
    print("pi(x) = A - B + C + D + Phi0 + Sigma");
    print_gourdon(x, y, z, k, threads);
  }

  // For very short computations (< 1 second) we achieve the best
  // performance by executing the different algorithms in increasing
  // order of their memory and power usage. This effect is mainly
  // measurable on server CPUs with a large number of CPU cores. On
  // my AMD EPYC 7642 server with 192 CPU cores I measured up to 2x
  // less context switches and cpu migrations using this method. If
  // we would start with the algorithm that puts the highest load on
  // the CPU and memory (i.e. the B algorithm) we would overload
  // both the CPU and operating system.

  int64_t lix = li(x);
  int64_t sigma = Sigma(x, y, threads, is_print);
  int64_t phi0 = Phi0(x, y, z, k, threads, is_print);
  int64_t ac = AC(x, y, z, k, threads, is_print);
  int64_t b = B(x, y, threads, is_print);
  int64_t d_approx = D_approx(x, lix, sigma, phi0, ac, b);
  int64_t d = D(x, y, z, k, d_approx, threads, is_print);
  int64_t pix = ac - b + d + phi0 + sigma;

  verify_pix("pi_gourdon_64", x, pix, lix);

  return pix;
}

#if defined(HAVE_INT128_T)

/// Calculate the number of primes below x using
/// Xavier Gourdon's algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int128_t pi_gourdon_128(int128_t x,
                        int threads,
                        bool is_print)
{
  if (x < 2)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if_unlikely(x > limit)
    throw primecount_error("pi(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t) 1);

  if (is_print)
  {
    print("");
    print("=== pi_gourdon_128(x) ===");
    print("pi(x) = A - B + C + D + Phi0 + Sigma");
    print_gourdon(x, y, z, k, threads);
  }

  // For very short computations (< 1 second) we achieve the best
  // performance by executing the different algorithms in increasing
  // order of their memory and power usage. This effect is mainly
  // measurable on server CPUs with a large number of CPU cores. On
  // my AMD EPYC 7642 server with 192 CPU cores I measured up to 2x
  // less context switches and cpu migrations using this method. If
  // we would start with the algorithm that puts the highest load on
  // the CPU and memory (i.e. the B algorithm) we would overload
  // both the CPU and operating system.

  int128_t lix = li(x);
  int128_t sigma = Sigma(x, y, threads, is_print);
  int128_t phi0 = Phi0(x, y, z, k, threads, is_print);
  int128_t ac = AC(x, y, z, k, threads, is_print);
  int128_t b = B(x, y, threads, is_print);
  int128_t d_approx = D_approx(x, lix, sigma, phi0, ac, b);
  int128_t d = D(x, y, z, k, d_approx, threads, is_print);
  int128_t pix = ac - b + d + phi0 + sigma;

  verify_pix("pi_gourdon_128", x, pix, lix);

  return pix;
}

#endif

} // namespace

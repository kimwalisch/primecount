///
/// @file  pi_gourdon.cpp
/// @brief Simple demonstration implementation of Xavier Gourdon's
///        prime counting algorithm. Xavier Gourdon's algorithm is
///        an improved version of the Deleglise-Rivat algorithm.
///
///        Xavier Gourdon formula:
///        pi(x) = A - B + C + D + phi0 + Sigma
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {

/// Calculate the number of primes below x using
/// Xavier Gourdon's algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_gourdon(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha_y = get_alpha_y_gourdon(x);
  double alpha_z = get_alpha_z_gourdon(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, isqrt(x) - 1);
  y = std::max(y, (int64_t) 1);

  int64_t k = PhiTiny::get_k(y);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, isqrt(x) - 1);
  z = std::max(z, (int64_t) 1);

  print("");
  print("=== pi_gourdon(x) ===");
  print("pi(x) = A - B + C + D + phi0 + Sigma");
  print(x, y, z, k, alpha_y, alpha_z, threads);

  int64_t sigma = Sigma(x, y, threads);
  int64_t phi0 = Phi0(x, y, z, k, threads);
  int64_t a = A(x, y, threads);
  int64_t b = B(x, y, threads);
  int64_t c = C(x, y, z, k, threads);
  int64_t d_approx = D_approx(x, sigma, phi0, a, b, c);
  int64_t d = D(x, y, z, k, d_approx, threads);
  int64_t sum = a - b + c + d + phi0 + sigma;

  return sum;
}

} // namespace

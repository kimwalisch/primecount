///
/// @file  Sigma.cpp
///        The 7 sigma formulas are the least computationally
///        expensive formulas in Gourdon's algorithm. Sigma0 has a
///        runtime complexity of O(x(1/2)), all other formulas
///        have a runtime complexity of O(y) and hence it does not
///        make much sense to use multi-threading.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <PiTable.hpp>
#include <print.hpp>

#include <stdint.h>

using namespace primecount;

namespace {

template <typename T>
T Sigma0(T x, T a, int threads)
{
  T pi_sqrtx = pi_noprint(isqrt(x), threads);
  return a - 1 + (pi_sqrtx * (pi_sqrtx - 1)) / 2 - (a * (a - 1)) / 2;
}

template <typename T>
T Sigma1(T a, T b)
{
  return (a - b) * (a - b - 1) / 2;
}

template <typename T>
T Sigma2(T a, T b, T c, T d)
{
  return a * (b - c - (c * (c - 3)) / 2 + (d * (d - 3)) / 2);
}

template <typename T>
T Sigma3(T b, T d)
{
  return (b * (b - 1) * (2 * b - 1)) / 6 - b - (d * (d - 1) * (2 * d - 1)) / 6 + d;
}

/// Memory usage: O(x^(3/8))
template <typename T>
T Sigma456(T x,
           int64_t y,
           int64_t a,
           int64_t x_star,
           const PiTable& pi)
{
  T sigma4 = 0;
  T sigma5 = 0;
  T sigma6 = 0;

  int64_t x13 = iroot<3>(x);
  int64_t sqrt_xy = isqrt(x / y);
  primesieve::iterator it(x_star, x13);
  int64_t prime = it.next_prime();

  // Sigma4: x_star < prime <= sqrt(x / y)
  // Sigma5: sqrt(x / y) < prime <= x^(1/3)
  // Sigma6: x_star < prime <= x^(1/3)
  for (; prime <= x13; prime = it.next_prime())
  {
    if (prime <= sqrt_xy)
      sigma4 += pi[x / (prime * (T) y)];
    else
      sigma5 += pi[x / (prime * (T) prime)];

    // Note that in Xavier Gourdon's paper the actual
    // formula for Σ6 is: sum += pi(x^(1/2) / prime^(1/2))^2.
    // However when implemented this way using integers
    // the formula returns incorrect results.
    // Hence the formula must be implemented as below:
    int64_t sqrt_xp = isqrt(x / prime);
    int64_t pi_sqrt_xp = pi[sqrt_xp];
    sigma6 += pi_sqrt_xp * (T) pi_sqrt_xp;
  }

  sigma4 *= a;
  sigma6 = -sigma6;

  return sigma4 + sigma5 + sigma6;
}

} // namespace

namespace primecount {

int64_t Sigma(int64_t x,
              int64_t y,
              int threads,
              bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== Sigma(x, y) ===");
    print_gourdon_vars(x, y, threads);
  }

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_pix_sigma4 = x / (x_star * y);
  int64_t max_pix_sigma5 = y;
  int64_t max_pix_sigma6 = isqrt(x / x_star);
  int64_t max_pix = max3(max_pix_sigma4, max_pix_sigma5, max_pix_sigma6);
  PiTable pi(max_pix, threads);

  int64_t a = pi[y];
  int64_t b = pi[iroot<3>(x)];
  int64_t c = pi[isqrt(x / y)];
  int64_t d = pi[x_star];

  int64_t sum = Sigma0(x, a, threads) +
                Sigma1(a, b) +
                Sigma2(a, b, c, d) +
                Sigma3(b, d) +
                Sigma456(x, y, a, x_star, pi);

  if (is_print)
    print("Sigma", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t Sigma(int128_t x,
               int64_t y,
               int threads,
               bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== Sigma(x, y) ===");
    print_gourdon_vars(x, y, threads);
  }

  int128_t x_star = get_x_star_gourdon(x, y);
  int64_t max_pix_sigma4 = x / (x_star * y);
  int64_t max_pix_sigma5 = y;
  int64_t max_pix_sigma6 = isqrt(x / x_star);
  int64_t max_pix = max3(max_pix_sigma4, max_pix_sigma5, max_pix_sigma6);
  PiTable pi(max_pix, threads);

  int128_t a = pi[y];
  int128_t b = pi[iroot<3>(x)];
  int128_t c = pi[isqrt(x / y)];
  int128_t d = pi[x_star];

  double time = get_time();
  int128_t sum = Sigma0(x, a, threads) +
                 Sigma1(a, b) +
                 Sigma2(a, b, c, d) +
                 Sigma3(b, d) +
                 Sigma456(x, y, a, x_star, pi);

  if (is_print)
    print("Sigma", sum, time);

  return sum;
}

#endif

} // namespace

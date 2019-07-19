///
/// @file  Sigma.cpp
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
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
#include <algorithm>

using namespace std;
using namespace primecount;

namespace {

template <typename T>
T Sigma0(T x, int64_t y, int threads)
{
  T a = pi_legendre(y, threads);
  T pi_sqrtx = pi_legendre(isqrt(x), threads);

  return a - 1 + (pi_sqrtx * (pi_sqrtx - 1)) / 2 - (a * (a - 1)) / 2;
}

template <typename T>
T Sigma1(T x, int64_t y, int threads)
{
  T a = pi_legendre(y, threads);
  T b = pi_legendre(iroot<3>(x), threads);

  return (a - b) * (a - b - 1) / 2;
}

template <typename T>
T Sigma2(T x, int64_t y, int threads)
{
  T a = pi_legendre(y, threads);
  T b = pi_legendre(iroot<3>(x), threads);
  T c = pi_legendre(isqrt(x / y), threads);
  T x_star = max(iroot<4>(x), x / ((T) y * (T) y));
  T d = pi_legendre(x_star, threads);

  return a * (b - c - (c * (c - 3)) / 2 + (d * (d - 3)) / 2);
}

template <typename T>
T Sigma3(T x, int64_t y, int threads)
{
  T b = pi_legendre(iroot<3>(x), threads);
  T x_star = max(iroot<4>(x), x / ((T) y * (T) y));
  T d = pi_legendre(x_star, threads);

  return (b * (b - 1) * (2 * b - 1)) / 6 - b - (d * (d - 1) * (2 * d - 1)) / 6 + d;
}

/// Memory usage: O(x^(1/3)) or less
template <typename T>
T Sigma4(T x, int64_t y, int threads)
{
  T sum = 0;
  int64_t pi_y = pi_legendre(y, threads);
  int64_t x_star = max(iroot<4>(x), x / ((T) y * (T) y));
  int64_t sqrt_xy = isqrt(x / y);

  PiTable pi(x / (x_star * (T) y));
  primesieve::iterator it(x_star, sqrt_xy);
  int64_t prime = it.next_prime();

  for (; prime <= sqrt_xy; prime = it.next_prime())
    sum += pi[x / (prime * (T) y)];

  return pi_y * sum;
}

/// Memory usage: O(y)
template <typename T>
T Sigma5(T x, int64_t y)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t sqrt_xy = isqrt(x / y);

  PiTable pi(y);
  primesieve::iterator it(sqrt_xy, x13);
  int64_t prime = it.next_prime();

  for (; prime <= x13; prime = it.next_prime())
    sum += pi[x / (prime * (T) prime)];

  return sum;
}

/// Memory usage: O(x^(3/8))
template <typename T>
T Sigma6(T x, int64_t y)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t x_star = max(iroot<4>(x), x / ((T) y * (T) y));

  PiTable pi(isqrt(x / x_star));
  primesieve::iterator it(x_star, x13);
  int64_t prime = it.next_prime();

  for (; prime <= x13; prime = it.next_prime())
  {
    // Note that in Xavier Gourdon's paper the actual
    // formula for Σ6 is: sum += pi(x^(1/2) / prime^(1/2))^2.
    // However when implemented this way using integers
    // the formula returns erroneous results.
    // Hence the formula must be implemented as below:
    T pix = pi[isqrt(x / prime)]; 
    sum += pix * pix;
  }
 
  return -sum;
}

} // namespace

namespace primecount {

int64_t Sigma(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== Sigma(x, y) ===");
  print(x, y, threads);

  double time = get_time();
  int64_t sum = Sigma0(x, y, threads) +
                Sigma1(x, y, threads) +
                Sigma2(x, y, threads) +
                Sigma3(x, y, threads) +
                Sigma4(x, y, threads) +
                Sigma5(x, y) +
                Sigma6(x, y);
  
  print("Sigma", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t Sigma(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== Sigma(x, y) ===");
  print(x, y, threads);

  double time = get_time();
  int128_t sum = Sigma0(x, y, threads) +
                 Sigma1(x, y, threads) +
                 Sigma2(x, y, threads) +
                 Sigma3(x, y, threads) +
                 Sigma4(x, y, threads) +
                 Sigma5(x, y) +
                 Sigma6(x, y);

  print("Sigma", sum, time);
  return sum;
}

#endif

} // namespace

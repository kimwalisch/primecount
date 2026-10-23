///
/// @file   D.cpp
/// @brief  Test the D function used in Gourdon's algorithm.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <gourdon.hpp>
#include <generate_primes.hpp>
#include <imath.hpp>
#include <PiTable.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <array>

using namespace primecount;

struct D_formula_params
{
  int64_t x;
  int64_t y;
  int64_t k;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon1.sh
/// For each input x=10^n we test using:
/// 1) The default alpha
/// 2) The minimum alpha=1
/// 3) The maximum alpha
std::array<D_formula_params, 69> test_cases =
{{
  { 10, 3, 0, 0 },
  { 10, 3, 0, 0 },
  { 10, 3, 0, 0 },
  { 100, 5, 2, 0 },
  { 100, 5, 2, 0 },
  { 100, 8, 2, 0 },
  { 1000, 15, 3, 0 },
  { 1000, 11, 3, 0 },
  { 1000, 30, 3, 0 },
  { 10000, 36, 4, 0 },
  { 10000, 22, 4, 0 },
  { 10000, 84, 4, 0 },
  { 100000, 87, 7, 0 },
  { 100000, 47, 7, 0 },
  { 100000, 276, 7, 0 },
  { 1000000, 207, 8, 2465 },
  { 1000000, 101, 8, 2465 },
  { 1000000, 1000, 8, 11215 },
  { 10000000, 485, 8, 132692 },
  { 10000000, 216, 8, 112667 },
  { 10000000, 3010, 8, 199723 },
  { 100000000, 1131, 8, 2413042 },
  { 100000000, 465, 8, 2141021 },
  { 100000000, 9744, 8, 2542718 },
  { 1000000000, 2619, 8, 30871820 },
  { 1000000000, 1001, 8, 30228636 },
  { 1000000000, 31000, 8, 27700392 },
  { 10000000000, 6029, 8, 351726346 },
  { 10000000000, 2155, 8, 365911138 },
  { 10000000000, 99084, 8, 279536758 },
  { 100000000000, 13825, 8, 3738964518 },
  { 100000000000, 4642, 8, 4018018477 },
  { 100000000000, 315588, 8, 2648288971 },
  { 1000000000000, 76500, 8, 32914684646 },
  { 1000000000000, 10001, 8, 42262337684 },
  { 1000000000000, 1000000, 8, 25185974064 },
  { 10000000000000, 241615, 8, 307179531887 },
  { 100000000000000, 660067, 8, 2902060663462 },
  { 1000000000000000, 1724700, 8, 27502746342885 },
  // Small inputs, including y = 1
  { 2, 1, 0, 0 },
  { 3, 1, 0, 0 },
  { 4, 2, 0, 0 },
  { 8, 2, 0, 0 },
  { 49, 7, 1, 0 },
  // The first hard leaf at x = 23^3 * 29 crosses the
  // prime-only/squarefree split at y = 23^2.
  { 352842, 528, 8, 0 },
  { 352843, 528, 8, 93 },
  { 352844, 528, 8, 93 },
  { 352842, 529, 8, 1617 },
  { 352843, 529, 8, 1618 },
  { 352844, 529, 8, 1618 },
  { 352842, 530, 8, 1617 },
  { 352843, 530, 8, 1618 },
  { 352844, 530, 8, 1618 },
  // y = sqrt(x) and its neighbors
  { 1000000, 999, 8, 11215 },
  { 100000000, 9999, 8, 2528919 },
  { 100000000, 10000, 8, 2528919 },
  { 148035889, 12167, 8, 4040900 },
  { 148035890, 12167, 8, 4040900 },
  { 148048055, 12167, 8, 4041242 },
  { 148048056, 12167, 8, 4041242 },
  // Prime fourth powers with the minimum y
  { 707280, 90, 8, 953 },
  { 707281, 90, 8, 953 },
  { 707282, 90, 8, 953 },
  { 923520, 98, 8, 2070 },
  { 923521, 98, 8, 2070 },
  { 923522, 98, 8, 2070 },
  { 88529280, 446, 8, 1836639 },
  { 88529281, 446, 8, 1836639 },
  { 88529282, 446, 8, 1836639 }
}};

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int64_t D_reference(int64_t x, int64_t y, int64_t k)
{
  auto primes = generate_primes<uint32_t>(y);
  auto lpf = generate_lpf(y, primes);
  auto mu = generate_moebius(y, primes);
  int64_t max_prime = iroot<4>(x);
  int64_t sum = 0;

  for (std::size_t b = k + 1; b < primes.size(); b++)
  {
    int64_t prime = primes[b];
    if (prime > max_prime) break;

    // Primes with x/p^3 < y and sqrt(x/p) < y < p^3 also
    // include their easy leaves (AC does not compute them).
    int64_t xp = x / prime;
    bool easy_leaves = prime * prime <= y &&
                        xp / (prime * prime) < y &&
                        prime * prime > y / prime &&
                        isqrt(xp) < y;

    for (int64_t m = prime + 1; m <= y; m++)
      if (mu[m] != 0 &&
          lpf[m] > prime &&
          prime * m > y &&
          (easy_leaves || x / (prime * m) >= prime * prime))
        sum -= mu[m] * phi(x / (prime * m), b - 1);
  }

  return sum;
}

int main()
{
  int threads = get_num_threads();

  for (const D_formula_params& params : test_cases)
  {
    const PiTable pi(params.y, threads);
    int64_t res = D(params.x, params.y, params.k, pi, threads);
    std::cout << "D_64bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = D((int128_t) params.x, params.y, params.k, pi, threads);
      std::cout << "D_128bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  // Check the transition from square-free to prime leaves, including
  // the first square-free composite m beyond each prime squared.
  const Array<int64_t, 5> primes = { 23, 29, 31, 37, 41 };

  for (std::size_t i = 0; i + 2 < primes.size(); i++)
  {
    int64_t p = primes[i];
    int64_t q = primes[i + 1];
    int64_t r = primes[i + 2];

    for (int64_t center : { ipow<5>(p), p * p * p * q * r })
      for (int64_t delta : { -1, 0, 1 })
        for (int64_t y : { p * p - 1, p * p, p * p + 1,
                           q * r - 1, q * r, q * r + 1 })
        {
          int64_t x = center + delta;
          int64_t expected = D_reference(x, y, 8);
          const PiTable pi(y, threads);
          int64_t res = D(x, y, 8, pi, threads);
          std::cout << "D boundary(" << x << ", " << y << ") = " << res;
          check(res == expected);

          #ifdef HAVE_INT128_T
            int128_t res2 = D((int128_t) x, y, 8, pi, threads);
            check(res2 == expected);
          #endif
        }
  }

  struct HyperbolaParams
  {
    int64_t xp;
    int64_t y;
    int64_t prime;
    int hyperbola;
  };

  // Hyperbola method: 0 = not used, 1 = also easy leaves (x/p^3 < y),
  // 2 = hard leaves only (x/p^3 >= y).
  // Include x/(p*y) = p^2, x/(p*y) = sqrt(x/p), prime/semiprime
  // endpoints, and boundaries immediately below, at, and above
  // multiples of the sieve alignment.
  const Array<HyperbolaParams, 30> hyperbola_cases = {{
    { 279840, 530, 23, 1 },
    { 279841, 530, 23, 1 },
    { 280370, 530, 23, 2 },
    { 280898, 530, 23, 2 },
    { 280899, 530, 23, 2 },
    { 280900, 530, 23, 0 },
    { 280901, 530, 23, 0 },
    { 559681, 1058, 23, 1 },
    { 559682, 1058, 23, 2 },
    { 560739, 1058, 23, 2 },
    { 560740, 1058, 23, 2 },
    { 7653060, 12167, 23, 0 },
    { 292140, 541, 23, 2 },
    { 807302, 899, 23, 2 },
    { 864017, 1200, 23, 2 },
    { 396017, 720, 23, 2 },
    { 444888, 767, 23, 2 },
    { 444889, 767, 23, 2 },
    { 444890, 767, 23, 2 },
    { 707280, 941, 23, 2 },
    { 707281, 941, 23, 2 },
    { 707282, 941, 23, 2 },
    { 808200, 999, 23, 2 },
    { 808201, 999, 23, 2 },
    { 808202, 999, 23, 2 },
    { 921599, 1200, 23, 2 },
    { 921600, 1200, 23, 2 },
    { 921601, 1200, 23, 2 },
    { 1414562, 1682, 29, 2 },
    { 1416243, 1682, 29, 2 }
  }};

  for (const HyperbolaParams& params : hyperbola_cases)
  {
    int64_t p = params.prime;
    int64_t x = params.xp * p;
    int64_t sqrt_xp = isqrt(params.xp);
    int hyperbola = 0;

    if (p * p <= params.y &&
        sqrt_xp < params.y &&
        params.y < p * p * p)
      hyperbola = (params.xp / (p * p) < params.y) ? 1 : 2;

    std::cout << "D hyperbola " << params.hyperbola << " ("
              << x << ", " << params.y << ")";
    check(hyperbola == params.hyperbola &&
          iroot<3>(x) < params.y &&
          params.y <= isqrt(x));
    int64_t expected = D_reference(x, params.y, 8);
    const PiTable pi(params.y, threads);
    check(D(x, params.y, 8, pi, threads) == expected);

    #ifdef HAVE_INT128_T
      check(D((int128_t) x, params.y, 8, pi, threads) == expected);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

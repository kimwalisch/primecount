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
std::array<D_formula_params, 65> test_cases =
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
  { 1000000, 1000, 8, 1246 },
  { 10000000, 485, 8, 132692 },
  { 10000000, 216, 8, 112667 },
  { 10000000, 3010, 8, 67108 },
  { 100000000, 1131, 8, 2413042 },
  { 100000000, 465, 8, 2141021 },
  { 100000000, 9744, 8, 1204711 },
  { 1000000000, 2619, 8, 30871820 },
  { 1000000000, 1001, 8, 30228636 },
  { 1000000000, 31000, 8, 15033924 },
  { 10000000000, 6029, 8, 351726346 },
  { 10000000000, 2155, 8, 365911138 },
  { 10000000000, 99084, 8, 158874158 },
  { 100000000000, 13825, 8, 3738964518 },
  { 100000000000, 4642, 8, 4018018477 },
  { 100000000000, 315588, 8, 1556900264 },
  { 1000000000000, 76500, 8, 32856197918 },
  { 1000000000000, 10001, 8, 42262337684 },
  { 1000000000000, 1000000, 8, 14815465134 },
  { 10000000000000, 241615, 8, 304739017794 },
  { 100000000000000, 660067, 8, 2876480095175 },
  { 1000000000000000, 1724700, 8, 27248002319287 },
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
  { 352842, 529, 8, 0 },
  { 352843, 529, 8, 93 },
  { 352844, 529, 8, 93 },
  { 352842, 530, 8, 0 },
  { 352843, 530, 8, 93 },
  { 352844, 530, 8, 93 },
  // y = sqrt(x) and its neighbors
  { 1000000, 999, 8, 1246 },
  { 100000000, 9999, 8, 1181490 },
  { 100000000, 10000, 8, 1181490 },
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

int64_t D_reference(int64_t x, int64_t y, int64_t z, int64_t k)
{
  auto primes = generate_primes<uint32_t>(z);
  auto lpf = generate_lpf(z, primes);
  auto mpf = generate_mpf(z, primes);
  auto mu = generate_moebius(z, primes);
  int64_t max_prime = iroot<4>(x);
  int64_t sum = 0;

  for (std::size_t b = k + 1; b < primes.size(); b++)
  {
    int64_t prime = primes[b];
    if (prime > max_prime) break;

    for (int64_t m = prime + 1; m <= z; m++)
      if (mu[m] != 0 &&
          lpf[m] > prime &&
          mpf[m] <= y &&
          prime * m > z &&
          x / (prime * m) >= prime * prime)
        sum -= mu[m] * phi(x / (prime * m), b - 1);
  }

  return sum;
}

int main()
{
  int threads = get_num_threads();

  for (const D_formula_params& params : test_cases)
  {
    int64_t res = D(params.x, params.y, params.k, threads);
    std::cout << "D_64bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = D((int128_t) params.x, params.y, params.k, threads);
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
        for (int64_t z : { p * p - 1, p * p, p * p + 1,
                           q * r - 1, q * r, q * r + 1 })
        {
          int64_t x = center + delta;
          int64_t y = iroot<3>(x) + 1;
          int64_t expected = D_reference(x, y, z, 8);
          int64_t res = D(x, y, z, 8, threads);
          std::cout << "D boundary(" << x << ", " << y << ", " << z << ") = " << res;
          check(res == expected);

          #ifdef HAVE_INT128_T
            int128_t res2 = D((int128_t) x, y, z, 8, threads);
            check(res2 == expected);
          #endif
        }
  }

  struct CentralParams
  {
    int64_t xp;
    int64_t y;
    int64_t z;
    int64_t prime;
    int region;
  };

  // Regions: 0 = fallback, 1 = natural upper bound, 2 = upper bound y.
  // Include A = p^2, A = H, prime/semiprime endpoints, and boundaries
  // immediately below, at, and above multiples of the sieve alignment.
  const Array<CentralParams, 33> central_cases = {{
    { 279840, 530, 530, 23, 0 },
    { 279841, 530, 530, 23, 0 },
    { 280370, 530, 530, 23, 1 },
    { 280898, 530, 530, 23, 1 },
    { 280899, 530, 530, 23, 2 },
    { 280900, 530, 530, 23, 0 },
    { 280901, 530, 530, 23, 0 },
    { 559681, 1058, 1058, 23, 1 },
    { 559682, 1058, 1058, 23, 1 },
    { 560739, 1058, 1058, 23, 2 },
    { 560740, 1058, 1058, 23, 2 },
    { 280370, 600, 700, 23, 1 },
    { 665499, 1058, 2116, 23, 0 },
    { 7653060, 12167, 12167, 23, 0 },
    { 292140, 541, 541, 23, 2 },
    { 807302, 899, 899, 23, 2 },
    { 864017, 1200, 1200, 23, 2 },
    { 396017, 720, 720, 23, 2 },
    { 444888, 767, 767, 23, 2 },
    { 444889, 767, 767, 23, 2 },
    { 444890, 767, 767, 23, 2 },
    { 707280, 941, 941, 23, 2 },
    { 707281, 941, 941, 23, 2 },
    { 707282, 941, 941, 23, 2 },
    { 808200, 999, 999, 23, 2 },
    { 808201, 999, 999, 23, 2 },
    { 808202, 999, 999, 23, 2 },
    { 921599, 1200, 1200, 23, 2 },
    { 921600, 1200, 1200, 23, 2 },
    { 921601, 1200, 1200, 23, 2 },
    { 1414562, 1682, 1682, 29, 1 },
    { 1416243, 1682, 1682, 29, 2 },
    { 1582779, 1682, 3364, 29, 0 }
  }};

  for (const CentralParams& params : central_cases)
  {
    int64_t p = params.prime;
    int64_t x = params.xp * p;
    int64_t upper = params.xp / (p * p);
    int64_t lower = params.xp / params.y;
    int64_t middle = isqrt(params.xp);
    int region = 0;

    if (p * p <= upper && middle < upper && upper <= params.y)
      region = 1;
    else if (params.z == params.y &&
             p * p <= lower && lower <= middle &&
             middle < params.y && params.y < p * p * p)
      region = 2;

    std::cout << "D central region " << params.region << " ("
              << x << ", " << params.y << ", " << params.z << ")";
    check(region == params.region &&
          iroot<3>(x) < params.y &&
          params.y <= params.z && params.z < isqrt(x));
    int64_t expected = D_reference(x, params.y, params.z, 8);
    check(D(x, params.y, params.z, 8, threads) == expected);

    #ifdef HAVE_INT128_T
      check(D((int128_t) x, params.y, params.z, 8, threads) == expected);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

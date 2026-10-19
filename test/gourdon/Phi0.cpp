///
/// @file   Phi0.cpp
/// @brief  Test the Phi0 function used in Gourdon's algorithm.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <gourdon.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <array>

using namespace primecount;

struct Phi0_formula_params
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
std::array<Phi0_formula_params, 39> test_cases =
{{
  { 10, 2, 0, 5 },
  { 10, 2, 0, 5 },
  { 10, 2, 0, 5 },
  { 100, 5, 2, 26 },
  { 100, 5, 2, 26 },
  { 100, 8, 2, 21 },
  { 1000, 15, 3, 184 },
  { 1000, 11, 3, 204 },
  { 1000, 30, 3, 134 },
  { 10000, 36, 4, 1396 },
  { 10000, 22, 4, 1647 },
  { 10000, 84, 4, 906 },
  { 100000, 87, 7, 11248 },
  { 100000, 47, 7, 13391 },
  { 100000, 276, 7, 7329 },
  { 1000000, 207, 8, 86595 },
  { 1000000, 101, 8, 109894 },
  { 1000000, 999, 8, 45589 },
  { 10000000, 485, 8, 623863 },
  { 10000000, 216, 8, 858009 },
  { 10000000, 3010, 8, 266214 },
  { 100000000, 1131, 8, 4221122 },
  { 100000000, 465, 8, 6312352 },
  { 100000000, 9744, 8, 1433229 },
  { 1000000000, 2619, 8, 28775469 },
  { 1000000000, 1001, 8, 45262927 },
  { 1000000000, 31000, 8, 7564137 },
  { 10000000000, 6029, 8, 186957171 },
  { 10000000000, 2155, 8, 311485874 },
  { 10000000000, 99084, 8, 37685290 },
  { 100000000000, 13825, 8, 1185193538 },
  { 100000000000, 4642, 8, 2145153126 },
  { 100000000000, 315588, 8, 186320545 },
  { 1000000000000, 76500, 8, 4423248789 },
  { 1000000000000, 10001, 8, 14168456261 },
  { 1000000000000, 999999, 8, 891878665 },
  { 10000000000000, 241615, 8, 21938749396 },
  { 100000000000000, 660067, 8, 116101684200 },
  { 1000000000000000, 1724700, 8, 618276610001 }
}};

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  int threads = get_num_threads();

  for (const Phi0_formula_params& params : test_cases)
  {
    int64_t res = Phi0(params.x, params.y, params.k, threads);
    std::cout << "Phi0_64bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = Phi0((int128_t) params.x, params.y, params.k, threads);
      std::cout << "Phi0_128bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

#ifdef HAVE_INT128_T
  {
    // Test Phi0(1e20) and compare with known correct value
    int128_t x = ((int128_t) 10000000000) * ((int128_t) 10000000000);
    int64_t y = 100615703;
    int64_t k = 8;
    int128_t res1 = Phi0(x, y, k, threads);
    int128_t res2 = 3742831815767347;

    std::cout << "Phi0_128bit(" << x << ", " << y << ", " << k << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

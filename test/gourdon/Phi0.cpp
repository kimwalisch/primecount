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
  int64_t z;
  int64_t k;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon1.sh
/// For each input x=10^n we test using:
/// 1) The default alpha_y & alpha_z
/// 2) The minimum alpha_y=1 & alpha_z=1
/// 3) The maximum alpha_y
/// 4) The maximum alpha_z
std::array<Phi0_formula_params, 51> test_cases =
{{
  { 10, 2, 2, 0, 5 },
  { 10, 2, 2, 0, 5 },
  { 10, 2, 2, 0, 5 },
  { 10, 2, 2, 0, 5 },
  { 100, 5, 5, 2, 26 },
  { 100, 5, 5, 2, 26 },
  { 100, 8, 8, 2, 21 },
  { 100, 5, 9, 2, 26 },
  { 1000, 15, 15, 3, 184 },
  { 1000, 11, 11, 3, 204 },
  { 1000, 30, 30, 3, 134 },
  { 1000, 11, 30, 3, 204 },
  { 10000, 36, 36, 4, 1396 },
  { 10000, 22, 22, 4, 1647 },
  { 10000, 84, 84, 4, 906 },
  { 10000, 22, 88, 4, 1647 },
  { 100000, 87, 87, 7, 11248 },
  { 100000, 47, 47, 7, 13391 },
  { 100000, 276, 276, 7, 7329 },
  { 100000, 47, 282, 7, 13391 },
  { 1000000, 207, 207, 8, 86595 },
  { 1000000, 101, 101, 8, 109894 },
  { 1000000, 999, 999, 8, 45589 },
  { 1000000, 101, 999, 8, 111121 },
  { 10000000, 485, 485, 8, 623863 },
  { 10000000, 216, 216, 8, 858009 },
  { 10000000, 3010, 3010, 8, 266214 },
  { 10000000, 216, 3024, 8, 934624 },
  { 100000000, 1131, 1131, 8, 4221122 },
  { 100000000, 465, 465, 8, 6312352 },
  { 100000000, 9744, 9744, 8, 1433229 },
  { 100000000, 465, 9765, 8, 8174210 },
  { 1000000000, 2619, 2619, 8, 28775469 },
  { 1000000000, 1001, 1001, 8, 45262927 },
  { 1000000000, 31000, 31000, 8, 7564137 },
  { 1000000000, 1001, 31031, 8, 75599746 },
  { 10000000000, 6029, 6029, 8, 186957171 },
  { 10000000000, 2155, 2155, 8, 311485874 },
  { 10000000000, 99084, 99084, 8, 37685290 },
  { 10000000000, 2155, 99130, 8, 709831674 },
  { 100000000000, 13825, 13825, 8, 1185193538 },
  { 100000000000, 4642, 4642, 8, 2145153126 },
  { 100000000000, 315588, 315588, 8, 186320545 },
  { 100000000000, 4642, 315656, 8, 6788948344 },
  { 1000000000000, 50000, 70850, 8, 10073346812 },
  { 1000000000000, 10001, 10001, 8, 14168456261 },
  { 1000000000000, 999999, 999999, 8, 891878665 },
  { 1000000000000, 10001, 999999, 8, 65140722670 },
  { 10000000000000, 107720, 209946, 8, 119423210693 },
  { 100000000000000, 282435, 564870, 8, 1045985238238 },
  { 1000000000000000, 737200, 1474400, 8, 9230903137263 }
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
    int64_t res = Phi0(params.x, params.y, params.z, params.k, threads);
    std::cout << "Phi0_64bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = Phi0((int128_t) params.x, params.y, params.z, params.k, threads);
      std::cout << "Phi0_128bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.k << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

#ifdef HAVE_INT128_T
  {
    // Test Phi0(1e20) and compare with known correct value
    int128_t x = ((int128_t) 10000000000) * ((int128_t) 10000000000);
    int64_t y = 100615703;
    int64_t z = 201231406;
    int64_t k = 8;
    int128_t res1 = Phi0(x, y, z, k, threads);
    int128_t res2 = 633772346752344505;

    std::cout << "Phi0_128bit(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

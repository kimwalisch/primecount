///
/// @file   Phi0.cpp
/// @brief  Test the Phi0 function used in Gourdon's algorithm.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
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
  { 10LL, 2, 2, 0, 5LL },
  { 10LL, 2, 2, 0, 5LL },
  { 10LL, 2, 2, 0, 5LL },
  { 10LL, 2, 2, 0, 5LL },
  { 100LL, 5, 5, 2, 26LL },
  { 100LL, 5, 5, 2, 26LL },
  { 100LL, 8, 8, 2, 21LL },
  { 100LL, 5, 9, 2, 26LL },
  { 1000LL, 15, 15, 3, 184LL },
  { 1000LL, 11, 11, 3, 204LL },
  { 1000LL, 30, 30, 3, 134LL },
  { 1000LL, 11, 30, 3, 204LL },
  { 10000LL, 36, 36, 4, 1396LL },
  { 10000LL, 22, 22, 4, 1647LL },
  { 10000LL, 84, 84, 4, 906LL },
  { 10000LL, 22, 88, 4, 1647LL },
  { 100000LL, 87, 87, 7, 11248LL },
  { 100000LL, 47, 47, 7, 13391LL },
  { 100000LL, 276, 276, 7, 7329LL },
  { 100000LL, 47, 282, 7, 13391LL },
  { 1000000LL, 207, 207, 8, 86595LL },
  { 1000000LL, 101, 101, 8, 109894LL },
  { 1000000LL, 999, 999, 8, 45589LL },
  { 1000000LL, 101, 999, 8, 111121LL },
  { 10000000LL, 485, 485, 8, 623863LL },
  { 10000000LL, 216, 216, 8, 858009LL },
  { 10000000LL, 3010, 3010, 8, 266214LL },
  { 10000000LL, 216, 3024, 8, 934624LL },
  { 100000000LL, 1131, 1131, 8, 4221122LL },
  { 100000000LL, 465, 465, 8, 6312352LL },
  { 100000000LL, 9744, 9744, 8, 1433229LL },
  { 100000000LL, 465, 9765, 8, 8174210LL },
  { 1000000000LL, 2619, 2619, 8, 28775469LL },
  { 1000000000LL, 1001, 1001, 8, 45262927LL },
  { 1000000000LL, 31000, 31000, 8, 7564137LL },
  { 1000000000LL, 1001, 31031, 8, 75599746LL },
  { 10000000000LL, 6029, 6029, 8, 186957171LL },
  { 10000000000LL, 2155, 2155, 8, 311485874LL },
  { 10000000000LL, 99084, 99084, 8, 37685290LL },
  { 10000000000LL, 2155, 99130, 8, 709831674LL },
  { 100000000000LL, 13825, 13825, 8, 1185193538LL },
  { 100000000000LL, 4642, 4642, 8, 2145153126LL },
  { 100000000000LL, 315588, 315588, 8, 186320545LL },
  { 100000000000LL, 4642, 315656, 8, 6788948344LL },
  { 1000000000000LL, 50000, 70850, 8, 10073346812LL },
  { 1000000000000LL, 10001, 10001, 8, 14168456261LL },
  { 1000000000000LL, 999999, 999999, 8, 891878665LL },
  { 1000000000000LL, 10001, 999999, 8, 65140722670LL },
  { 10000000000000LL, 107720, 209946, 8, 119423210693LL },
  { 100000000000000LL, 282435, 564870, 8, 1045985238238LL },
  { 1000000000000000LL, 737200, 1474400, 8, 9230903137263LL }
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
    int128_t res2 = 633772346752344505ll;

    std::cout << "Phi0_128bit(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

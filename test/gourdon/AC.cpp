///
/// @file   AC.cpp
/// @brief  Test the AC function used in Gourdon's algorithm.
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

struct AC_formula_params
{
  int64_t x;
  int64_t y;
  int64_t z;
  int64_t k;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon1.sh
/// For each input x we test using:
/// 1) The default alpha_y & alpha_z
/// 2) The maximum alpha_y
/// 3) The maximum alpha_z
std::array<AC_formula_params, 36> test_cases =
{{
  { 10LL, 2, 2, 0, 0LL },
  { 10LL, 2, 2, 0, 0LL },
  { 10LL, 2, 2, 0, 0LL },
  { 100LL, 5, 5, 2, 0LL },
  { 100LL, 8, 8, 2, 0LL },
  { 100LL, 5, 9, 2, 0LL },
  { 1000LL, 15, 15, 3, 10LL },
  { 1000LL, 30, 30, 3, 10LL },
  { 1000LL, 11, 30, 3, 3LL },
  { 10000LL, 36, 36, 4, 170LL },
  { 10000LL, 84, 84, 4, 258LL },
  { 10000LL, 22, 88, 4, 64LL },
  { 100000LL, 87, 87, 7, 1331LL },
  { 100000LL, 276, 276, 7, 1886LL },
  { 100000LL, 47, 282, 7, 507LL },
  { 1000000LL, 207, 207, 8, 18065LL },
  { 1000000LL, 999, 999, 8, 27607LL },
  { 1000000LL, 101, 999, 8, 7197LL },
  { 10000000LL, 485, 485, 8, 175136LL },
  { 10000000LL, 3010, 3010, 8, 322447LL },
  { 10000000LL, 216, 3024, 8, 64127LL },
  { 100000000LL, 1131, 1131, 8, 1563000LL },
  { 100000000LL, 9744, 9744, 8, 3077951LL },
  { 100000000LL, 465, 9765, 8, 548239LL },
  { 1000000000LL, 2619, 2619, 8, 13875464LL },
  { 1000000000LL, 31000, 31000, 8, 28111648LL },
  { 1000000000LL, 1001, 31031, 8, 4838942LL },
  { 10000000000LL, 6029, 6029, 8, 124680594LL },
  { 10000000000LL, 99084, 99084, 8, 257272174LL },
  { 10000000000LL, 2155, 99130, 8, 45096604LL },
  { 100000000000LL, 13825, 13825, 8, 1122177179LL },
  { 100000000000LL, 315588, 315588, 8, 2378181717LL },
  { 100000000000LL, 4642, 315656, 8, 412323998LL },
  { 10000000000000LL, 107720, 209946, 8, 106430408717LL },
  { 100000000000000LL, 282435, 564870, 8, 1008985328656LL },
  { 1000000000000000LL, 737200, 1474400, 8, 9561261537251LL }
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

  for (const AC_formula_params& params : test_cases)
  {
    int64_t res = AC(params.x, params.y, params.z, params.k, threads);
    std::cout << "AC_64bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = AC((int128_t) params.x, params.y, params.z, params.k, threads);
      std::cout << "AC_128bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.k << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

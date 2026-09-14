///
/// @file   AC.cpp
/// @brief  Test the AC function used in Gourdon's algorithm.
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

struct AC_formula_params
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
std::array<AC_formula_params, 51> test_cases =
{{
  { 10, 2, 2, 0, 0 },
  { 10, 2, 2, 0, 0 },
  { 10, 2, 2, 0, 0 },
  { 10, 2, 2, 0, 0 },
  { 100, 5, 5, 2, 0 },
  { 100, 5, 5, 2, 0 },
  { 100, 8, 8, 2, 0 },
  { 100, 5, 9, 2, 0 },
  { 1000, 15, 15, 3, 10 },
  { 1000, 11, 11, 3, 3 },
  { 1000, 30, 30, 3, 10 },
  { 1000, 11, 30, 3, 3 },
  { 10000, 36, 36, 4, 170 },
  { 10000, 22, 22, 4, 64 },
  { 10000, 84, 84, 4, 258 },
  { 10000, 22, 88, 4, 64 },
  { 100000, 87, 87, 7, 1331 },
  { 100000, 47, 47, 7, 507 },
  { 100000, 276, 276, 7, 1886 },
  { 100000, 47, 282, 7, 507 },
  { 1000000, 207, 207, 8, 18065 },
  { 1000000, 101, 101, 8, 7205 },
  { 1000000, 999, 999, 8, 27607 },
  { 1000000, 101, 999, 8, 7197 },
  { 10000000, 485, 485, 8, 175136 },
  { 10000000, 216, 216, 8, 68724 },
  { 10000000, 3010, 3010, 8, 322447 },
  { 10000000, 216, 3024, 8, 64127 },
  { 100000000, 1131, 1131, 8, 1563000 },
  { 100000000, 465, 465, 8, 657446 },
  { 100000000, 9744, 9744, 8, 3077951 },
  { 100000000, 465, 9765, 8, 548239 },
  { 1000000000, 2619, 2619, 8, 13875464 },
  { 1000000000, 1001, 1001, 8, 6023539 },
  { 1000000000, 31000, 31000, 8, 28111648 },
  { 1000000000, 1001, 31031, 8, 4838942 },
  { 10000000000, 6029, 6029, 8, 124680594 },
  { 10000000000, 2155, 2155, 8, 56822547 },
  { 10000000000, 99084, 99084, 8, 257272174 },
  { 10000000000, 2155, 99130, 8, 45096604 },
  { 100000000000, 13825, 13825, 8, 1122177179 },
  { 100000000000, 4642, 4642, 8, 526077035 },
  { 100000000000, 315588, 315588, 8, 2378181717 },
  { 100000000000, 4642, 315656, 8, 412323998 },
  { 1000000000000, 50000, 70850, 8, 12040548047 },
  { 1000000000000, 10001, 10001, 8, 4906487209 },
  { 1000000000000, 999999, 999999, 8, 21948644076 },
  { 1000000000000, 10001, 999999, 8, 3808940423 },
  { 10000000000000, 107720, 209946, 8, 106430408717 },
  { 100000000000000, 282435, 564870, 8, 1008985328656 },
  { 1000000000000000, 737200, 1474400, 8, 9561261537251 }
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

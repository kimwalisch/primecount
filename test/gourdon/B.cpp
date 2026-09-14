///
/// @file   B.cpp
/// @brief  Test the B function used in Gourdon's algorithm.
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

struct B_formula_params
{
  int64_t x;
  int64_t y;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon2.sh
/// For each input x=10^n we test using:
/// 1) The default alpha_y & alpha_z
/// 2) The maximum alpha_y
/// 3) The maximum alpha_z
std::array<B_formula_params, 41> test_cases =
{{
  { 10, 2, 2 },
  { 10, 2, 2 },
  { 10, 2, 2 },
  { 100, 5, 6 },
  { 100, 8, 0 },
  { 100, 5, 6 },
  { 1000, 15, 67 },
  { 1000, 30, 11 },
  { 1000, 11, 88 },
  { 10000, 36, 543 },
  { 10000, 84, 56 },
  { 10000, 22, 761 },
  { 100000, 87, 4403 },
  { 100000, 276, 480 },
  { 100000, 47, 6295 },
  { 1000000, 207, 37293 },
  { 1000000, 999, 0 },
  { 1000000, 101, 54794 },
  { 10000000, 485, 325348 },
  { 10000000, 3010, 6887 },
  { 10000000, 216, 473021 },
  { 100000000, 1131, 2876542 },
  { 100000000, 9744, 33602 },
  { 100000000, 465, 4100054 },
  { 1000000000, 2619, 25991893 },
  { 1000000000, 31000, 209274 },
  { 1000000000, 1001, 36435407 },
  { 10000000000, 6029, 235385820 },
  { 10000000000, 99084, 770317 },
  { 10000000000, 2155, 325113158 },
  { 100000000000, 13825, 2151216255 },
  { 100000000000, 315588, 1420565 },
  { 100000000000, 4642, 2943439103 },
  { 1000000000000, 50000, 17133805730 },
  { 1000000000000, 999999, 0 },
  { 1000000000000, 10001, 26809544511 },
  { 10000000000000, 107720, 163974930685 },
  { 10000000000000, 3145424, 255862065 },
  { 10000000000000, 21545, 246427408287 },
  { 100000000000000, 282435, 1483796135572 },
  { 1000000000000000, 737200, 13558621700511 }
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

  for (const B_formula_params& params : test_cases)
  {
    int64_t res = B(params.x, params.y, threads);
    std::cout << "B_64bit(" << params.x << ", " << params.y << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = B((int128_t) params.x, params.y, threads);
      std::cout << "B_128bit(" << params.x << ", " << params.y << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

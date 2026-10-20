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
std::array<D_formula_params, 39> test_cases =
{{
  { 10, 2, 0, 0 },
  { 10, 2, 0, 0 },
  { 10, 2, 0, 0 },
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
  { 1000000000000000, 1724700, 8, 27248002319287 }
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

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

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
  int64_t k;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon1.sh
/// For each input x=10^n we test using:
/// 1) The default alpha
/// 2) The minimum alpha=1
/// 3) The maximum alpha
std::array<AC_formula_params, 43> test_cases =
{{
  { 10, 3, 0, 0 },
  { 10, 3, 0, 0 },
  { 10, 3, 0, 0 },
  { 100, 5, 2, 0 },
  { 100, 5, 2, 0 },
  { 100, 8, 2, 0 },
  { 1000, 15, 3, 10 },
  { 1000, 11, 3, 3 },
  { 1000, 30, 3, 10 },
  { 10000, 36, 4, 170 },
  { 10000, 22, 4, 64 },
  { 10000, 84, 4, 258 },
  { 100000, 87, 7, 1331 },
  { 100000, 47, 7, 507 },
  { 100000, 276, 7, 1886 },
  { 1000000, 207, 8, 18065 },
  { 1000000, 101, 8, 7205 },
  { 1000000, 1000, 8, 27607 },
  { 10000000, 485, 8, 175136 },
  { 10000000, 216, 8, 68724 },
  { 10000000, 3010, 8, 322447 },
  { 100000000, 1131, 8, 1563000 },
  { 100000000, 465, 8, 657446 },
  { 100000000, 9744, 8, 3077951 },
  { 1000000000, 2619, 8, 13875464 },
  { 1000000000, 1001, 8, 6023539 },
  { 1000000000, 31000, 8, 28111648 },
  { 10000000000, 6029, 8, 124680594 },
  { 10000000000, 2155, 8, 56822547 },
  { 10000000000, 99084, 8, 257272174 },
  { 100000000000, 13825, 8, 1122177179 },
  { 100000000000, 4642, 8, 526077035 },
  { 100000000000, 315588, 8, 2378181717 },
  { 1000000000000, 76500, 8, 13752431086 },
  { 1000000000000, 10001, 8, 4906487209 },
  { 1000000000000, 1000000, 8, 21948644076 },
  { 10000000000000, 241615, 8, 134251530576 },
  { 100000000000000, 660067, 8, 1260415861796 },
  { 1000000000000000, 1724700, 8, 11750232830271 },
  // C2 reflected endpoint: x / 23 crosses 131^2 and 131 * 132.
  { 394702, 149, 8, 6475 },
  { 394703, 149, 8, 6476 },
  { 397715, 149, 8, 6570 },
  { 397716, 149, 8, 6570 }
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
    int64_t res = AC(params.x, params.y, params.k, threads);
    std::cout << "AC_64bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = AC((int128_t) params.x, params.y, params.k, threads);
      std::cout << "AC_128bit(" << params.x << ", " << params.y << ", " << params.k << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

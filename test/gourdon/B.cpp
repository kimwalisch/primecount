///
/// @file   B.cpp
/// @brief  Test the B function used in Gourdon's algorithm.
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
  { 10LL, 2, 2LL },
  { 10LL, 2, 2LL },
  { 10LL, 2, 2LL },
  { 100LL, 5, 6LL },
  { 100LL, 8, 0LL },
  { 100LL, 5, 6LL },
  { 1000LL, 15, 67LL },
  { 1000LL, 30, 11LL },
  { 1000LL, 11, 88LL },
  { 10000LL, 36, 543LL },
  { 10000LL, 84, 56LL },
  { 10000LL, 22, 761LL },
  { 100000LL, 87, 4403LL },
  { 100000LL, 276, 480LL },
  { 100000LL, 47, 6295LL },
  { 1000000LL, 207, 37293LL },
  { 1000000LL, 999, 0LL },
  { 1000000LL, 101, 54794LL },
  { 10000000LL, 485, 325348LL },
  { 10000000LL, 3010, 6887LL },
  { 10000000LL, 216, 473021LL },
  { 100000000LL, 1131, 2876542LL },
  { 100000000LL, 9744, 33602LL },
  { 100000000LL, 465, 4100054LL },
  { 1000000000LL, 2619, 25991893LL },
  { 1000000000LL, 31000, 209274LL },
  { 1000000000LL, 1001, 36435407LL },
  { 10000000000LL, 6029, 235385820LL },
  { 10000000000LL, 99084, 770317LL },
  { 10000000000LL, 2155, 325113158LL },
  { 100000000000LL, 13825, 2151216255LL },
  { 100000000000LL, 315588, 1420565LL },
  { 100000000000LL, 4642, 2943439103LL },
  { 1000000000000LL, 50000, 17133805730LL },
  { 1000000000000LL, 999999, 0LL },
  { 1000000000000LL, 10001, 26809544511LL },
  { 10000000000000LL, 107720, 163974930685LL },
  { 10000000000000LL, 3145424, 255862065LL },
  { 10000000000000LL, 21545, 246427408287LL },
  { 100000000000000LL, 282435, 1483796135572LL },
  { 1000000000000000LL, 737200, 13558621700511LL }
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

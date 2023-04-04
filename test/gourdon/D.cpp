///
/// @file   D.cpp
/// @brief  Test the D function used in Gourdon's algorithm.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
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
  int64_t z;
  int64_t k;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon1.sh
/// For each input x we test using:
/// 1) The default alpha_y & alpha_z
/// 2) The maximum alpha_y
/// 3) The maximum alpha_z
std::array<D_formula_params, 33> test_cases =
{{
  { 10LL, 2, 2, 0, 0LL },
  { 10LL, 2, 2, 0, 0LL },
  { 10LL, 2, 2, 0, 0LL },
  { 100LL, 5, 5, 2, 0LL },
  { 100LL, 8, 8, 2, 0LL },
  { 100LL, 5, 9, 2, 0LL },
  { 1000LL, 15, 15, 3, 0LL },
  { 1000LL, 30, 30, 3, 0LL },
  { 1000LL, 11, 30, 3, 0LL },
  { 10000LL, 36, 36, 4, 0LL },
  { 10000LL, 84, 84, 4, 0LL },
  { 10000LL, 22, 88, 4, 0LL },
  { 100000LL, 87, 87, 7, 0LL },
  { 100000LL, 276, 276, 7, 0LL },
  { 100000LL, 47, 282, 7, 0LL },
  { 1000000LL, 207, 207, 8, 2465LL },
  { 1000000LL, 999, 999, 8, 1246LL },
  { 1000000LL, 101, 999, 8, 1246LL },
  { 10000000LL, 485, 485, 8, 132692LL },
  { 10000000LL, 3010, 3010, 8, 67108LL },
  { 10000000LL, 216, 3024, 8, 40649LL },
  { 100000000LL, 1131, 1131, 8, 2413042LL },
  { 100000000LL, 9744, 9744, 8, 1204711LL },
  { 100000000LL, 465, 9765, 8, 388370LL },
  { 1000000000LL, 2619, 2619, 8, 30871820LL },
  { 1000000000LL, 31000, 31000, 8, 15033924LL },
  { 1000000000LL, 1001, 31031, 8, 1076414LL },
  { 10000000000LL, 6029, 6029, 8, 351726346LL },
  { 10000000000LL, 99084, 99084, 8, 158874158LL },
  { 10000000000LL, 2155, 99130, 8, -20708719LL },
  { 10000000000000LL, 107720, 209946, 8, 270354670695LL },
  { 100000000000000LL, 282435, 564870, 8, 2518169986968LL },
  { 1000000000000000LL, 737200, 1474400, 8, 23628309295271LL }
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
    int64_t res = D(params.x, params.y, params.z, params.k, Li(params.x), threads);
    std::cout << "D_64bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.k << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = D((int128_t) params.x, params.y, params.z, params.k, (int128_t) Li(params.x), threads);
      std::cout << "D_128bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.k << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

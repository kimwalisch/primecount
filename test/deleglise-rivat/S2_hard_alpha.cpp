///
/// @file   S2_hard_alpha.cpp
/// @brief  Test the S2_hard function used in the
///         Deleglise-Rivat algorithm.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <S.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <array>

using namespace primecount;

struct formula_params
{
  int64_t x;
  int64_t y;
  int64_t z;
  int64_t c;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_dr.sh
/// For each input x we test using:
/// 1) The default alpha
/// 2) The minimum alpha=1
/// 3) The maximum alpha
std::array<formula_params, 35> test_cases =
{{
  { 10LL, 2, 5, 1, 0LL },
  { 10LL, 2, 5, 1, 0LL },
  { 10LL, 2, 5, 1, 0LL },
  { 100LL, 5, 20, 3, 0LL },
  { 100LL, 4, 25, 2, 0LL },
  { 100LL, 8, 12, 4, 0LL },
  { 1000LL, 15, 66, 6, 0LL },
  { 1000LL, 10, 100, 4, 0LL },
  { 1000LL, 30, 33, 8, 0LL },
  { 10000LL, 36, 277, 8, 0LL },
  { 10000LL, 21, 476, 8, 0LL },
  { 10000LL, 84, 119, 8, 0LL },
  { 100000LL, 87, 1149, 8, 185LL },
  { 100000LL, 46, 2173, 8, 242LL },
  { 100000LL, 276, 362, 8, 0LL },
  { 1000000LL, 207, 4830, 8, 11557LL },
  { 1000000LL, 100, 10000, 8, 9171LL },
  { 1000000LL, 1000, 1000, 8, 11215LL },
  { 10000000LL, 485, 20618, 8, 233493LL },
  { 10000000LL, 215, 46511, 8, 181391LL },
  { 10000000LL, 3010, 3322, 8, 199723LL },
  { 100000000LL, 1131, 88417, 8, 3353160LL },
  { 100000000LL, 464, 215517, 8, 2798467LL },
  { 100000000LL, 9744, 10262, 8, 2542718LL },
  { 1000000000LL, 2619, 381825, 8, 39599180LL },
  { 1000000000LL, 1000, 1000000, 8, 36252175LL },
  { 1000000000LL, 31000, 32258, 8, 27700392LL },
  { 10000000000LL, 10621, 941530, 8, 389393048LL },
  { 10000000000LL, 2154, 4642525, 8, 422733685LL },
  { 10000000000LL, 99084, 100924, 8, 279536758LL },
  { 100000000000LL, 25766, 3881083, 8, 4000882894LL },
  { 100000000000LL, 4641, 21547080, 8, 4544095512LL },
  { 100000000000LL, 315588, 316868, 8, 2648288971LL },
  { 10000000000000LL, 178815, 55923720, 8, 371595005834LL },
  { 100000000000000LL, 494134, 202374254, 8, 3474606376629LL }
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

  for (const formula_params& params : test_cases)
  {
    int64_t res = S2_hard(params.x, params.y, params.z, params.c, Li(params.x), threads);
    std::cout << "S2_hard_64bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.c << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = S2_hard((int128_t) params.x, params.y, params.z, params.c, (int128_t) Li(params.x), threads);
      std::cout << "S2_hard_128bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.c << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

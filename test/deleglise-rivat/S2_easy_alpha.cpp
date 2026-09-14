///
/// @file   S2_easy_alpha.cpp
/// @brief  Test the S2_easy function used in the
///         Deleglise-Rivat algorithm.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
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
/// For each input x=10^n we test using:
/// 1) The default alpha
/// 2) The minimum alpha=1
/// 3) The maximum alpha
std::array<formula_params, 35> test_cases =
{{
  { 10, 2, 5, 1, 0 },
  { 10, 2, 5, 1, 0 },
  { 10, 2, 5, 1, 0 },
  { 100, 5, 20, 3, 0 },
  { 100, 4, 25, 2, 0 },
  { 100, 8, 12, 4, 0 },
  { 1000, 15, 66, 6, 0 },
  { 1000, 10, 100, 4, 0 },
  { 1000, 30, 33, 8, 0 },
  { 10000, 36, 277, 8, 0 },
  { 10000, 21, 476, 8, 0 },
  { 10000, 84, 119, 8, 0 },
  { 100000, 87, 1149, 8, 328 },
  { 100000, 46, 2173, 8, 0 },
  { 100000, 276, 362, 8, 618 },
  { 1000000, 207, 4830, 8, 4330 },
  { 1000000, 100, 10000, 8, 0 },
  { 1000000, 1000, 1000, 8, 9895 },
  { 10000000, 485, 20618, 8, 36104 },
  { 10000000, 215, 46511, 8, 0 },
  { 10000000, 3010, 3322, 8, 116253 },
  { 100000000, 1131, 88417, 8, 320032 },
  { 100000000, 464, 215517, 8, 0 },
  { 100000000, 9744, 10262, 8, 1110900 },
  { 1000000000, 2619, 381825, 8, 2725380 },
  { 1000000000, 1000, 1000000, 8, 0 },
  { 1000000000, 31000, 32258, 8, 10224606 },
  { 10000000000, 10621, 941530, 8, 69354279 },
  { 10000000000, 2154, 4642525, 8, 0 },
  { 10000000000, 99084, 100924, 8, 93607845 },
  { 100000000000, 25766, 3881083, 8, 622734970 },
  { 100000000000, 4641, 21547080, 8, 0 },
  { 100000000000, 315588, 316868, 8, 917197198 },
  { 10000000000000, 178815, 55923720, 8, 60888055472 },
  { 100000000000000, 494134, 202374254, 8, 617442826127 }
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
    int64_t res = S2_easy(params.x, params.y, params.z, params.c, threads);
    std::cout << "S2_easy_64bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.c << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = S2_easy((int128_t) params.x, params.y, params.z, params.c, threads);
      std::cout << "S2_easy_128bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.c << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

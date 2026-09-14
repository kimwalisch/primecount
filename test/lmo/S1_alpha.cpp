///
/// @file   S1_alpha.cpp
/// @brief  Test the S1 function used in the
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
  int64_t c;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_s1.sh
/// For each input x=10^n we test using:
/// 1) The default alpha
/// 2) The minimum alpha=1
/// 3) The maximum alpha
std::array<formula_params, 45> test_cases =
{{
  { 10, 2, 1, 5 },
  { 10, 2, 1, 5 },
  { 10, 2, 1, 5 },
  { 100, 5, 3, 26 },
  { 100, 4, 2, 33 },
  { 100, 8, 4, 22 },
  { 1000, 15, 6, 190 },
  { 1000, 10, 4, 228 },
  { 1000, 30, 8, 159 },
  { 10000, 36, 8, 1514 },
  { 10000, 21, 8, 1711 },
  { 10000, 84, 8, 1111 },
  { 100000, 87, 8, 11582 },
  { 100000, 46, 8, 13962 },
  { 100000, 276, 8, 7829 },
  { 1000000, 207, 8, 86595 },
  { 1000000, 100, 8, 111589 },
  { 1000000, 1000, 8, 45589 },
  { 10000000, 485, 8, 623863 },
  { 10000000, 215, 8, 858009 },
  { 10000000, 3010, 8, 266214 },
  { 100000000, 1131, 8, 4221122 },
  { 100000000, 464, 8, 6312352 },
  { 100000000, 9744, 8, 1433229 },
  { 1000000000, 2619, 8, 28775469 },
  { 1000000000, 1000, 8, 45262927 },
  { 1000000000, 31000, 8, 7564137 },
  { 10000000000, 10621, 8, 137353002 },
  { 10000000000, 2154, 8, 311485874 },
  { 10000000000, 99084, 8, 37685290 },
  { 100000000000, 25766, 8, 837262286 },
  { 100000000000, 4641, 8, 2145153126 },
  { 100000000000, 315588, 8, 186320545 },
  { 1000000000000, 66380, 8, 4843382848 },
  { 1000000000000, 10000, 8, 14151355571 },
  { 1000000000000, 1000000, 8, 891878665 },
  { 10000000000000, 178815, 8, 26584164412 },
  { 10000000000000, 21544, 8, 92959518290 },
  { 10000000000000, 3145424, 8, 4150683115 },
  { 100000000000000, 494134, 8, 138733088111 },
  { 100000000000000, 46415, 8, 598951445224 },
  { 100000000000000, 9979225, 8, 18992123716 },
  { 1000000000000000, 1378500, 8, 714283960231 },
  { 1000000000000000, 100000, 8, 3740876972904 },
  { 1000000000000000, 31600000, 8, 83148668863 }
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
    int64_t res = S1(params.x, params.y, params.c, threads);
    std::cout << "S1_64bit(" << params.x << ", " << params.y << ", " << params.c << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = S1((int128_t) params.x, params.y, params.c, threads);
      std::cout << "S1_128bit(" << params.x << ", " << params.y << ", " << params.c << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

#ifdef HAVE_INT128_T
  {
    // Test S1(1e20) and compare with known correct value
    int128_t x = ((int128_t) 10000000000) * ((int128_t) 10000000000);
    int64_t y = 209809060;
    int64_t c = 8;
    int128_t res1 = S1(x, y, c, threads);
    int128_t res2 = 2141872489903326;

    std::cout << "S1_128bit(" << x << ", " << y << ", " << c << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

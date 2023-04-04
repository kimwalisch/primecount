///
/// @file   S1_alpha.cpp
/// @brief  Test the S1 function used in the
///         Deleglise-Rivat algorithm.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
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
  { 10LL, 2, 1, 5LL },
  { 10LL, 2, 1, 5LL },
  { 10LL, 2, 1, 5LL },
  { 100LL, 5, 3, 26LL },
  { 100LL, 4, 2, 33LL },
  { 100LL, 8, 4, 22LL },
  { 1000LL, 15, 6, 190LL },
  { 1000LL, 10, 4, 228LL },
  { 1000LL, 30, 8, 159LL },
  { 10000LL, 36, 8, 1514LL },
  { 10000LL, 21, 8, 1711LL },
  { 10000LL, 84, 8, 1111LL },
  { 100000LL, 87, 8, 11582LL },
  { 100000LL, 46, 8, 13962LL },
  { 100000LL, 276, 8, 7829LL },
  { 1000000LL, 207, 8, 86595LL },
  { 1000000LL, 100, 8, 111589LL },
  { 1000000LL, 1000, 8, 45589LL },
  { 10000000LL, 485, 8, 623863LL },
  { 10000000LL, 215, 8, 858009LL },
  { 10000000LL, 3010, 8, 266214LL },
  { 100000000LL, 1131, 8, 4221122LL },
  { 100000000LL, 464, 8, 6312352LL },
  { 100000000LL, 9744, 8, 1433229LL },
  { 1000000000LL, 2619, 8, 28775469LL },
  { 1000000000LL, 1000, 8, 45262927LL },
  { 1000000000LL, 31000, 8, 7564137LL },
  { 10000000000LL, 10621, 8, 137353002LL },
  { 10000000000LL, 2154, 8, 311485874LL },
  { 10000000000LL, 99084, 8, 37685290LL },
  { 100000000000LL, 25766, 8, 837262286LL },
  { 100000000000LL, 4641, 8, 2145153126LL },
  { 100000000000LL, 315588, 8, 186320545LL },
  { 1000000000000LL, 66380, 8, 4843382848LL },
  { 1000000000000LL, 10000, 8, 14151355571LL },
  { 1000000000000LL, 1000000, 8, 891878665LL },
  { 10000000000000LL, 178815, 8, 26584164412LL },
  { 10000000000000LL, 21544, 8, 92959518290LL },
  { 10000000000000LL, 3145424, 8, 4150683115LL },
  { 100000000000000LL, 494134, 8, 138733088111LL },
  { 100000000000000LL, 46415, 8, 598951445224LL },
  { 100000000000000LL, 9979225, 8, 18992123716LL },
  { 1000000000000000LL, 1378500, 8, 714283960231LL },
  { 1000000000000000LL, 100000, 8, 3740876972904LL },
  { 1000000000000000LL, 31600000, 8, 83148668863LL }
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
    int128_t res2 = 2141872489903326ll;

    std::cout << "S1_128bit(" << x << ", " << y << ", " << c << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

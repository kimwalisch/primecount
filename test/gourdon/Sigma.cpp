///
/// @file   Sigma.cpp
/// @brief  Test the Sigma function used in Gourdon's algorithm.
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

struct Sigma_formula_params
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
std::array<Sigma_formula_params, 41> test_cases =
{{
  { 10, 2, 1 },
  { 10, 2, 1 },
  { 10, 2, 1 },
  { 100, 5, 5 },
  { 100, 8, 4 },
  { 100, 5, 5 },
  { 1000, 15, 41 },
  { 1000, 30, 35 },
  { 1000, 11, 49 },
  { 10000, 36, 206 },
  { 10000, 84, 121 },
  { 10000, 22, 279 },
  { 100000, 87, 1416 },
  { 100000, 276, 857 },
  { 100000, 47, 1989 },
  { 1000000, 207, 8666 },
  { 1000000, 999, 4056 },
  { 1000000, 101, 13728 },
  { 10000000, 485, 58236 },
  { 10000000, 3010, 15697 },
  { 10000000, 216, 98200 },
  { 100000000, 1131, 440833 },
  { 100000000, 9744, 79166 },
  { 100000000, 465, 750690 },
  { 1000000000, 2619, 3316674 },
  { 1000000000, 31000, 347099 },
  { 1000000000, 1001, 5767839 },
  { 10000000000, 6029, 27074220 },
  { 10000000000, 99084, 1991206 },
  { 10000000000, 2155, 45946110 },
  { 100000000000, 13825, 222935833 },
  { 100000000000, 315588, -1927148 },
  { 100000000000, 4642, 372245278 },
  { 1000000000000, 50000, 1541740088 },
  { 1000000000000, 999999, -48075857 },
  { 1000000000000, 10001, 3080175375 },
  { 10000000000000, 107720, 13832177419 },
  { 10000000000000, 3145424, -749805339 },
  { 10000000000000, 21545, 25908547161 },
  { 100000000000000, 282435, 115597332512 },
  { 1000000000000000, 737200, 982718153395 }
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

  for (const Sigma_formula_params& params : test_cases)
  {
    int64_t res = Sigma(params.x, params.y, threads);
    std::cout << "Sigma_64bit(" << params.x << ", " << params.y << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = Sigma((int128_t) params.x, params.y, threads);
      std::cout << "Sigma_128bit(" << params.x << ", " << params.y << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

#ifdef HAVE_INT128_T
  {
    // Test Sigma(1e20) and compare with known correct value
    int128_t x = ((int128_t) 10000000000) * ((int128_t) 10000000000);
    int64_t y = 100615703;
    int128_t res1 = Sigma(x, y, threads);
    int128_t res2 = 49384621237095387;

    std::cout << "Sigma(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

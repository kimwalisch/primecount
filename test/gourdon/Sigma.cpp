///
/// @file   Sigma.cpp
/// @brief  Test the Sigma function used in Gourdon's algorithm.
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

struct Sigma_formula_params
{
  int64_t x;
  int64_t y;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_gourdon2.sh
/// For each input x we test using:
/// 1) The default alpha_y & alpha_z
/// 2) The maximum alpha_y
/// 3) The maximum alpha_z
std::array<Sigma_formula_params, 41> test_cases =
{{
  { 10LL, 2, 1LL },
  { 10LL, 2, 1LL },
  { 10LL, 2, 1LL },
  { 100LL, 5, 5LL },
  { 100LL, 8, 4LL },
  { 100LL, 5, 5LL },
  { 1000LL, 15, 41LL },
  { 1000LL, 30, 35LL },
  { 1000LL, 11, 49LL },
  { 10000LL, 36, 206LL },
  { 10000LL, 84, 121LL },
  { 10000LL, 22, 279LL },
  { 100000LL, 87, 1416LL },
  { 100000LL, 276, 857LL },
  { 100000LL, 47, 1989LL },
  { 1000000LL, 207, 8666LL },
  { 1000000LL, 999, 4056LL },
  { 1000000LL, 101, 13728LL },
  { 10000000LL, 485, 58236LL },
  { 10000000LL, 3010, 15697LL },
  { 10000000LL, 216, 98200LL },
  { 100000000LL, 1131, 440833LL },
  { 100000000LL, 9744, 79166LL },
  { 100000000LL, 465, 750690LL },
  { 1000000000LL, 2619, 3316674LL },
  { 1000000000LL, 31000, 347099LL },
  { 1000000000LL, 1001, 5767839LL },
  { 10000000000LL, 6029, 27074220LL },
  { 10000000000LL, 99084, 1991206LL },
  { 10000000000LL, 2155, 45946110LL },
  { 100000000000LL, 13825, 222935833LL },
  { 100000000000LL, 315588, -1927148LL },
  { 100000000000LL, 4642, 372245278LL },
  { 1000000000000LL, 50000, 1541740088LL },
  { 1000000000000LL, 999999, -48075857LL },
  { 1000000000000LL, 10001, 3080175375LL },
  { 10000000000000LL, 107720, 13832177419LL },
  { 10000000000000LL, 3145424, -749805339LL },
  { 10000000000000LL, 21545, 25908547161LL },
  { 100000000000000LL, 282435, 115597332512LL },
  { 1000000000000000LL, 737200, 982718153395LL }
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
    int128_t res2 = 49384621237095387ll;

    std::cout << "Sigma(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

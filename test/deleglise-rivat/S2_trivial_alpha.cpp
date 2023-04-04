///
/// @file   S2_trivial_alpha.cpp
/// @brief  Test the S2_trivial function used in the
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
  int64_t z;
  int64_t c;
  int64_t res;
};

/// Known correct results generated using: scripts/gen_tests_dr.sh
/// For each input x we test using:
/// 1) The default alpha
/// 2) The maximum alpha
std::array<formula_params, 24> test_cases =
{{
  { 10LL, 2, 5, 1, 0LL },
  { 10LL, 2, 5, 1, 0LL },
  { 100LL, 5, 20, 3, 0LL },
  { 100LL, 8, 12, 4, 0LL },
  { 1000LL, 15, 66, 6, 0LL },
  { 1000LL, 30, 33, 8, 1LL },
  { 10000LL, 36, 277, 8, 3LL },
  { 10000LL, 84, 119, 8, 105LL },
  { 100000LL, 87, 1149, 8, 51LL },
  { 100000LL, 276, 362, 8, 1141LL },
  { 1000000LL, 207, 4830, 8, 271LL },
  { 1000000LL, 1000, 1000, 8, 11632LL },
  { 10000000LL, 485, 20618, 8, 1327LL },
  { 10000000LL, 3010, 3322, 8, 82276LL },
  { 100000000LL, 1131, 88417, 8, 6655LL },
  { 100000000LL, 9744, 10262, 8, 674204LL },
  { 1000000000LL, 2619, 381825, 8, 29329LL },
  { 1000000000LL, 31000, 32258, 8, 5358764LL },
  { 10000000000LL, 10621, 941530, 8, 574931LL },
  { 10000000000LL, 99084, 100924, 8, 44219304LL },
  { 100000000000LL, 25766, 3881083, 8, 2935021LL },
  { 100000000000LL, 315588, 316868, 8, 366223566LL },
  { 10000000000000LL, 178815, 55923720, 8, 110007115LL },
  { 100000000000000LL, 494134, 202374254, 8, 742709619LL }
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
    int64_t res = S2_trivial(params.x, params.y, params.z, params.c, threads);
    std::cout << "S2_trivial_64bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.c << ") = " << res;
    check(res == params.res);

    #ifdef HAVE_INT128_T
      int128_t res2 = S2_trivial((int128_t) params.x, params.y, params.z, params.c, threads);
      std::cout << "S2_trivial_128bit(" << params.x << ", " << params.y << ", " << params.z << ", " << params.c << ") = " << res2;
      check(res2 == params.res);
    #endif
  }

#ifdef HAVE_INT128_T
  {
    // Test S2_trivial(1e20) and compare with known correct value
    int128_t x = ((int128_t) 10000000000) * ((int128_t) 10000000000);
    int64_t y = 209809060;
    int64_t z = 476623840743;
    int64_t c = 8;
    int128_t res1 = S2_trivial(x, y, z, c, threads);
    int128_t res2 = 66066585011132ll;

    std::cout << "S2_trivial_128bit(" << x << ", " << y << ", " << z << ", " << c << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

///
/// @file   AC.cpp
/// @brief  Test the AC function used in Gourdon's algorithm.
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
#include <random>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  int threads = get_num_threads();

  {
    // Test AC(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 107720;
    int64_t z = 209946;
    int64_t k = 8;
    int64_t res1 = AC(x, y, z, k, threads);
    int64_t res2 = 106430408717ll;

    std::cout << "AC(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test AC(1e14) and compare with known correct value
    int128_t x = 100000000000000ll;
    int64_t y = 282435;
    int64_t z = 564870;
    int64_t k = 8;
    int128_t res1 = AC(x, y, z, k, threads);
    int128_t res2 = 1008985328656ll;

    std::cout << "AC(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

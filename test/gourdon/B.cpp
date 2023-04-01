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
    // Test B(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 107720;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 163974930685ll;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test B(1e14) and compare with known correct value
    int128_t x = 100000000000000ll;
    int64_t y = 282435;
    int128_t res1 = B(x, y, threads);
    int128_t res2 = 1483796135572ll;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

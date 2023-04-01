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
    // Test B(24) and compare with known correct value
    int64_t x = 24;
    int64_t y = 3;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 0;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(25) and compare with known correct value
    int64_t x = 25;
    int64_t y = 3;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 3;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(1e2) and compare with known correct value
    int64_t x = 100;
    int64_t y = 5;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 6;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(1e3) and compare with known correct value
    int64_t x = 1000;
    int64_t y = 15;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 67;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(1e5) and compare with known correct value
    int64_t x = 100000;
    int64_t y = 87;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 4403;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(1e7) and compare with known correct value
    int64_t x = 10000000;
    int64_t y = 323;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 397078;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 107720;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 163974930685ll;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test B(1e14) and compare with known correct value
    int64_t x = 100000000000000ll;
    int64_t y = 282435;
    int64_t res1 = B(x, y, threads);
    int64_t res2 = 1483796135572ll;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test B(1e15) and compare with known correct value
    int128_t x = 1000000000000000ll;
    int64_t y = 737200;
    int128_t res1 = B(x, y, threads);
    int128_t res2 = 13558621700511ll;

    std::cout << "B(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

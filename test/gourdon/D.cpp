///
/// @file   D.cpp
/// @brief  Test the D function used in Gourdon's algorithm.
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
    // Test D(1e2) and compare with known correct value
    int64_t x = 100;
    int64_t y = 5;
    int64_t z = 5;
    int64_t k = 2;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 0;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test D(1e5) and compare with known correct value
    int64_t x = 100000;
    int64_t y = 87;
    int64_t z = 87;
    int64_t k = 7;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 0;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test D(352843) and compare with known correct value
    int64_t x = 352843;
    int64_t y = 93;
    int64_t z = 139;
    int64_t k = 8;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 93;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test D(1e6) and compare with known correct value
    int64_t x = 1000000;
    int64_t y = 207;
    int64_t z = 207;
    int64_t k = 8;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 2465;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test D(1e7) and compare with known correct value
    int64_t x = 10000000;
    int64_t y = 323;
    int64_t z = 484;
    int64_t k = 8;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 125610;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test D(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 107720;
    int64_t z = 209946;
    int64_t k = 8;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 270354670695ll;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test D(1e14) and compare with known correct value
    int64_t x = 100000000000000ll;
    int64_t y = 282435;
    int64_t z = 564870;
    int64_t k = 8;
    int64_t res1 = D(x, y, z, k, Li(x), threads);
    int64_t res2 = 2518169986968ll;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test D(1e15) and compare with known correct value
    int128_t x = 1000000000000000ll;
    int64_t y = 737200;
    int64_t z = 1474400;
    int64_t k = 8;
    int128_t res1 = D(x, y, z, k, Li(x), threads);
    int128_t res2 = 23628309295271ll;

    std::cout << "D(" << x << ", " << y << ", " << z << ", " << k << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

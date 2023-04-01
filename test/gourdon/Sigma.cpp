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
    // Test Sigma(1e15) and compare with known correct value
    int64_t x = 1000000000000000ll;
    int64_t y = 737200;
    int64_t res1 = Sigma(x, y, threads);
    int64_t res2 = 982718153395ll;

    std::cout << "Sigma(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
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

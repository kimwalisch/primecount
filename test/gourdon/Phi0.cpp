///
/// @file   Phi0.cpp
/// @brief  Test the Phi0 function used in Gourdon's algorithm.
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
    // Test Phi0(1e15) and compare with known correct value
    int64_t x = 1000000000000000ll;
    int64_t y = 737200;
    int64_t z = 1474400;
    int64_t k = 8;
    int64_t res1 = Phi0(x, y, z, k, threads);
    int64_t res2 = 9230903137263ll;

    std::cout << "Phi0(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test Phi0(1e20) and compare with known correct value
    int128_t x = ((int128_t) 10000000000) * ((int128_t) 10000000000);
    int64_t y = 100615703;
    int64_t z = 201231406;
    int64_t k = 8;
    int128_t res1 = Phi0(x, y, z, k, threads);
    int128_t res2 = 633772346752344505ll;

    std::cout << "Phi0(" << x << ", " << y << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

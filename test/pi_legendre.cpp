///
/// @file   pi_legendre.cpp
/// @brief  Test the pi_legendre(x) function.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <PiTable.hpp>

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

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(0, 1 << 20);

  {
    int64_t x = -1;
    int64_t res = pi_legendre(x, threads);
    std::cout << "pi_legendre(" << x << ") = " << res;
    check(res == 0);
  }

  for (int64_t x = 0; x <= PiTable::max_cached(); x++)
  {
    int64_t res1 = pi_legendre(x, threads);
    int64_t res2 = pi_cache(x);
    std::cout << "pi_legendre(" << x << ") = " << res1;
    check(res1 == res2);
  }

  for (int i = 0; i < 100; i++)
  {
    int64_t x = dist(gen);
    int64_t res1 = pi_legendre(x, threads);
    int64_t res2 = pi_primesieve(x);
    std::cout << "pi_legendre(" << x << ") = " << res1;
    check(res1 == res2);
  }

  {
    // Test one larger computation: pi(1e11)
    int64_t x = 100000000000ll;
    int64_t res = pi_legendre(x, threads);
    std::cout << "pi_legendre(" << x << ") = " << res;
    check(res == 4118054813ll);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

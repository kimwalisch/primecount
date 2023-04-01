///
/// @file   P2_xa.cpp
/// @brief  Test the 2nd partial sieve function P2(x, a)
///         that counts the numbers <= x that have exactly
///         2 prime factors each exceeding the a-th prime.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using std::size_t;
using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(50000, 70000);

  int threads = 1;
  int64_t x = dist(gen);
  auto primes = generate_primes<int64_t>(x);

  for (int a = 1; primes[a] <= isqrt(x); a++)
  {
    int64_t p2 = 0;

    for (size_t b = a + 1; b < primes.size(); b++)
      for (size_t c = b; c < primes.size(); c++)
        if (primes[b] * primes[c] <= x)
          p2++;

    std::cout << "P2(" << x << ", " << a << ") = " << p2;
    check(p2 == P2(x, primes[a], a, threads));
  }

  threads = get_num_threads();

  {
    // Test B(1e13) and compare with known correct value
    int64_t x = 10000000000000ll;
    int64_t y = 178815;
    int64_t a = 16229;
    int64_t res1 = P2(x, y, a, threads);
    int64_t res2 = 113111712222ll;

    std::cout << "P2(" << x << ", " << y << ", " << a << ") = " << res1;
    check(res1 == res2);
  }

#ifdef HAVE_INT128_T
  {
    // Test P2(1e14) and compare with known correct value
    int128_t x = 100000000000000ll;
    int64_t y = 494134;
    int64_t a = 41080;
    int128_t res1 = P2(x, y, a, threads);
    int128_t res2 = 1026583290763ll;

    std::cout << "P2(" << x << ", " << y << ", " << a << ") = " << res1;
    check(res1 == res2);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

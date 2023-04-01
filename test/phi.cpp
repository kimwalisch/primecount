///
/// @file   phi.cpp
/// @brief  Test the partial sieve function phi(x, a)
///         which counts the numbers <= x that are not divisible
///         by any of the first a primes.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using std::size_t;
using namespace primecount;

void check(size_t x, size_t a, size_t phi_xa, size_t cnt)
{
  bool OK = (phi_xa == cnt);
  std::cout << "phi(" << x << ", " << a << ") = " << phi_xa;
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

void check2(size_t x, size_t a, size_t phi_xa, size_t cnt)
{
  if (phi_xa != cnt)
  {
    bool OK = (phi_xa == cnt);
    std::cout << "phi(" << x << ", " << a << ") = " << phi_xa;
    std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
    std::exit(1);
  }
  // Reduce logging because it is slow
  else if (a % 101 == 0)
  {
    bool OK = (phi_xa == cnt);
    std::cout << "phi(" << x << ", " << a << ") = " << phi_xa;
    std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  }
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());

  {
    std::uniform_int_distribution<size_t> dist(20000000, 30000000);
    size_t size = dist(gen);
    size_t x = size - 1;
    size_t cnt = size - 1;
    primesieve::iterator it;
    std::vector<char> sieve(size, 1);

    // test with small a values
    for (size_t a = 1;; a++)
    {
      auto prime = it.next_prime();
      if (prime * prime > x)
        break;

      // remove primes[a] and its multiples
      for (auto j = prime; j <= x; j += prime)
      {
        cnt -= (sieve[j] == 1);
        sieve[j] = 0;
      }

      int64_t phi_xa = phi(x, a);
      check(x, a, phi_xa, cnt);
    }
  }

  {
    std::uniform_int_distribution<size_t> dist(100000, 200000);
    size_t size = dist(gen);
    size_t x = size - 1;
    size_t cnt = size - 1;
    primesieve::iterator it;
    std::vector<char> sieve(size, 1);

    // test with large a values
    for (size_t a = 1;; a++)
    {
      auto prime = it.next_prime();
      if (prime > x)
        break;

      // remove primes[a] and its multiples
      for (auto j = prime; j <= x; j += prime)
      {
        cnt -= (sieve[j] == 1);
        sieve[j] = 0;
      }

      int64_t phi_xa = phi(x, a);
      check2(x, a, phi_xa, cnt);
    }
  }

  {
    std::cout << "Testing phi(x, a) multi-threading" << std::endl;

    int64_t iters = 500;
    int64_t sum1 = 0;
    int64_t sum2 = 0;

    #pragma omp parallel for reduction(+: sum1)
    for (int64_t i = 0; i < iters; i++)
      sum1 += pi_legendre(10000000 + i, 1);

    for (int64_t i = 0; i < iters; i++)
      sum2 += pi_legendre(10000000 + i, 1);

    if (sum1 == sum2)
    {
      std::cout << "Multi-thread sum: " << sum1 << " == Single-thread sum: " << sum2 << "   OK" << std::endl;
      std::cout << "phi(x, a) multi-threading: no data races detected!" << std::endl;
    }
    else
    {
      std::cout << "Multi-thread sum: " << sum1 << " != Single-thread sum: " << sum2 << "   ERROR" << std::endl;
      std::exit(1);
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

///
/// @file   P3.cpp
/// @brief  Test the 3rd partial sieve function P3(x, a)
///         that counts the numbers <= x that have exactly
///         3 prime factors each exceeding the a-th prime.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

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
  // Test small values of x
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(2, 1000);

    for (int i = 0; i < 100; i++)
    {
      int threads = 1;
      int64_t x = dist(gen);
      auto primes = generate_primes<int64_t>(x);

      for (int64_t a = 1; primes[a] <= iroot<3>(x); a++)
      {
        int64_t p3 = 0;

        for (size_t b = a + 1; b < primes.size(); b++)
        {
          for (size_t c = b; c < primes.size(); c++)
          {
            for (size_t d = c; d < primes.size(); d++)
            {
              if (primes[b] * primes[c] * primes[d] <= x)
                p3++;
              else
                break;
            }
          }
        }

        std::cout << "P3(" << x << ", " << a << ") = " << p3;
        check(p3 == P3(x, primes[a], a, threads));
      }
    }
  }

  // Test medium values of x
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(10000, 20000);

    for (int i = 0; i < 10; i++)
    {
      int threads = 1;
      int64_t x = dist(gen);
      auto primes = generate_primes<int64_t>(x);

      for (int64_t a = 1; primes[a] <= iroot<3>(x); a++)
      {
        int64_t p3 = 0;

        for (size_t b = a + 1; b < primes.size(); b++)
        {
          for (size_t c = b; c < primes.size(); c++)
          {
            for (size_t d = c; d < primes.size(); d++)
            {
              if (primes[b] * primes[c] * primes[d] <= x)
                p3++;
              else
                break;
            }
          }
        }

        std::cout << "P3(" << x << ", " << a << ") = " << p3;
        check(p3 == P3(x, primes[a], a, threads));
      }
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

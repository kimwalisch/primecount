///
/// @file   generate_lpf.cpp
/// @brief  Test least prime factor function
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <generate.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
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
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(200000, 300000);

  auto max = dist(gen);
  auto lpf = generate_lpf(max);
  auto primes = generate_primes<int32_t>(max);

  for (int i = 2; i <= max; i++)
  {
    int factor = i;
    int sqrt = isqrt(i);

    // find smallest prime factor
    for (int j = 1; primes[j] <= sqrt; j++)
    {
      if (i % primes[j] == 0)
      {
        factor = primes[j];
        break;
      }
    }

    std::cout << "lpf(" << i << ") = " << lpf[i];
    check(lpf[i] == factor);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

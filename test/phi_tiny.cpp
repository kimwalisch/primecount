///
/// @file   phi_tiny.cpp
/// @brief  Test the partial sieve function phi_tiny(x, a)
///         which counts the numbers <= x that are not divisible
///         by any of the first a primes with a <= 8.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>
#include <generate.hpp>

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

// Count the number of unsieved elements
int count(std::vector<char>& sieve)
{
  int cnt = 0;

  for (size_t i = 1; i < sieve.size(); i++)
    cnt += sieve[i];

  return cnt;
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(10000000, 20000000);

  int64_t max_a = PhiTiny::max_a();
  int64_t size = dist(gen);
  int64_t x = size - 1;

  auto primes = generate_n_primes<int32_t>(max_a);
  std::vector<char> sieve(size, 1);

  for (int a = 1; a <= max_a; a++)
  {
    // remove primes[a] and its multiples
    for (int j = primes[a]; j <= x; j += primes[a])
      sieve[j] = 0;

    std::cout << "phi_tiny(" << x << ", " << a << ") = " << phi_tiny(x, a);
    check(phi_tiny(x, a) == count(sieve));
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

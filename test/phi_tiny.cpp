///
/// @file   phi_tiny.cpp
/// @brief  Test the partial sieve function phi_tiny(x, a)
///         which counts the numbers <= x that are not divisible
///         by any of the first a primes with a <= 8.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>
#include <generate_primes.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
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
uint64_t count(const Vector<char>& sieve)
{
  uint64_t cnt = 0;

  for (size_t i = 1; i < sieve.size(); i++)
    cnt += sieve[i];

  return cnt;
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint64_t> dist(10000000, 20000000);

  uint64_t max_a = PhiTiny::max_a();
  uint64_t size = dist(gen);
  uint64_t x = size - 1;

  const Array<uint32_t, 9> primes = { 0, 2, 3, 5, 7, 11, 13, 17, 19 };
  ASSERT(max_a < primes.size());
  Vector<char> sieve(size);
  std::fill(sieve.begin(), sieve.end(), 1);

  for (uint64_t a = 1; a <= max_a; a++)
  {
    // remove primes[a] and its multiples
    for (uint64_t j = primes[a]; j <= x; j += primes[a])
      sieve[j] = 0;

    std::cout << "phi_tiny(" << x << ", " << a << ") = " << phi_tiny(x, a);
    check(phi_tiny(x, a) == count(sieve));
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

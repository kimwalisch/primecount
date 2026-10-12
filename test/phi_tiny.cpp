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
#include <algorithm>
#include <cstdlib>
#include <iostream>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  uint64_t max_a = PhiTiny::max_a();
  uint64_t max_x = PhiTiny::prime_products.back() * 2;

  const Array<uint32_t, 9> primes = { 0, 2, 3, 5, 7, 11, 13, 17, 19 };
  ASSERT(max_a < primes.size());
  Vector<char> sieve(max_x + 1);
  std::fill(sieve.begin(), sieve.end(), 1);

  for (uint64_t a = 1; a <= max_a; a++)
  {
    // remove primes[a] and its multiples
    for (uint64_t j = primes[a]; j <= max_x; j += primes[a])
      sieve[j] = 0;

    // Check all entries of PhiTiny's lookup tables
    uint64_t cnt = 0;
    bool OK = (phi_tiny(0, a) == 0);

    for (uint64_t i = 1; i <= max_x; i++)
    {
      cnt += sieve[i];
      OK &= (phi_tiny(i, a) == cnt);
    }

    std::cout << "phi_tiny(i, " << a << ") for i <= " << max_x;
    check(OK);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

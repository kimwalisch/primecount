///
/// @file  PhiTiny.cpp
/// @brief phi_tiny(x, a) calculates the partial sieve function in
///        constant time for small values of a <= 6 using lookup
///        tables. Let pp = prime_products_[a]:
///        phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a).
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>

#include <stdint.h>
#include <vector>
#include <cassert>

namespace {

const int32_t primes[7] = { 0, 2, 3, 5, 7, 11, 13 };

/// prime_products[n] = \prod_{i=1}^{n} primes[i]
const int32_t prime_products[7] = { 1, 2, 6, 30, 210, 2310, 30030 };

/// totients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const int32_t totients[7] = { 1, 1, 2, 8, 48, 480, 5760 };

const primecount::PhiTiny phiTiny;

}

namespace primecount {

int64_t phi_tiny(int64_t x, int64_t a)
{
  return phiTiny.phi(x, a);
}

PhiTiny::PhiTiny()
{
  phi_cache_[0].push_back(0);

  // Initialize the phi_cache_ lookup tables
  for (int a = 1; a <= max_a(); a++)
  {
    int size = prime_products[a];
    std::vector<int16_t>& cache = phi_cache_[a];
    cache.reserve(size);

    for (int x = 0; x < size; x++)
    {
      int16_t phixa = (int16_t) (phi(x, a - 1) - phi(x / primes[a], a - 1));
      cache.push_back(phixa);
    }
  }
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
/// @pre a <= 6.
///
int64_t PhiTiny::phi(int64_t x, int64_t a) const
{
  assert(x >= 0);
  assert(a <= max_a());
  return (x / prime_products[a]) * totients[a] + phi_cache_[a][x % prime_products[a]];
}

} // namespace primecount

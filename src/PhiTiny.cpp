///
/// @file  PhiTiny.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny.hpp"

#include <stdint.h>
#include <vector>
#include <cassert>

namespace {

/// Thread-safe singleton without locking
/// @note Meyer's singleton (static const Singleton instance) causes a
///       race condition with MSVC 2013 and OpenMP
///
const primecount::PhiTiny phiTiny;

}

namespace primecount {

int64_t phi_tiny(int64_t x, int64_t a)
{
  return phiTiny.phi(x, a);
}

const int64_t PhiTiny::MAX_A = 6;

const int32_t PhiTiny::primes_[7] = { 0, 2, 3, 5, 7, 11, 13 };

/// prime_products_[n] = \prod_{i=1}^{n} primes_[i]
const int32_t PhiTiny::prime_products_[7] = { 1, 2, 6, 30, 210, 2310, 30030 };

/// totients_[n] = \prod_{i=1}^{n} (primes_[i] - 1)
const int32_t PhiTiny::totients_[7] = { 1, 1, 2, 8, 48, 480, 5760 };

PhiTiny::PhiTiny()
{
  phi_cache_[0].push_back(0);

  // Initialize the phi_cache_ lookup tables
  for (int a = 1; a <= MAX_A; a++)
  {
    int size = prime_products_[a];
    std::vector<int16_t>& cache = phi_cache_[a];
    cache.reserve(size);

    for (int x = 0; x < size; x++)
    {
      int16_t phixa = static_cast<int16_t>(phi(x, a - 1) - phi(x / primes_[a], a - 1));
      cache.push_back(phixa);
    }
  }
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
/// @pre is_cached(a) == true.
///
int64_t PhiTiny::phi(int64_t x, int64_t a) const
{
  assert(x >= 0);
  assert(a <= MAX_A);
  return (x / prime_products_[a]) * totients_[a] + phi_cache_[a][x % prime_products_[a]];
}

} // namespace primecount

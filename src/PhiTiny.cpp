///
/// @file  PhiTiny.cpp
/// @brief phi_tiny(x, a) calculates the partial sieve function in
///        constant time (using lookup tables) for small values of
///        a <= 6 using the formula below:
///
///        phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
///        with pp = 2 * 3 * ... * prime[a]
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <array>
#include <vector>

namespace primecount {

const std::array<int, 7> PhiTiny::primes = { 0, 2, 3, 5, 7, 11, 13 };

// prime_products[n] = \prod_{i=1}^{n} primes[i]
const std::array<int, 7> PhiTiny::prime_products = { 1, 2, 6, 30, 210, 2310, 30030 };

// totients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const std::array<int, 7> PhiTiny::totients = { 1, 1, 2, 8, 48, 480, 5760 };

// Number of primes below x
const std::array<int, 13> PhiTiny::pi = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5 };

// Singleton
const PhiTiny phiTiny;

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
      int64_t phi_xa = phi(x, a - 1) - phi(x / primes[a], a - 1);
      cache.push_back((int16_t) phi_xa);
    }
  }
}

} // namespace

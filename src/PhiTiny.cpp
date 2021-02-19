///
/// @file  PhiTiny.cpp
/// @brief phi(x, a) counts the numbers <= x that are not divisible
///        by any of the first a primes. PhiTiny computes phi(x, a)
///        in constant time for a <= 6 using lookup tables.
///
///        phi(x, a) = (x / pp) * φ(a) + phi(x % pp, a)
///        pp = 2 * 3 * ... * prime[a]
///        φ(a) = \prod_{i=1}^{a} (prime[i] - 1)
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>

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
  // initialize phi(x % pp, a) lookup tables
  for (int a = 0; a <= max_a(); a++)
  {
    int pp = prime_products[a];
    phi_[a].resize(pp);
    phi_[a][0] = 0;

    for (int x = 1; x < pp; x++)
    {
      auto phi_xa = phi(x, a - 1) - phi(x / primes[a], a - 1);
      phi_[a][x] = (int16_t) phi_xa;
    }
  }
}

} // namespace

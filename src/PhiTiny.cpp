///
/// @file  PhiTiny.cpp
/// @brief phi(x, a) counts the numbers <= x that are not divisible
///        by any of the first a primes. PhiTiny computes phi(x, a)
///        in constant time for a <= 8 using lookup tables.
///
///        phi(x, a) = (x / pp) * φ(a) + phi(x % pp, a)
///        pp = 2 * 3 * ... * prime[a]
///        φ(a) = \prod_{i=1}^{a} (prime[i] - 1)
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>
#include <popcnt.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <array>
#include <vector>

namespace primecount {

const std::array<int, 9> PhiTiny::primes = { 0, 2, 3, 5, 7, 11, 13, 17, 19 };

// prime_products[n] = \prod_{i=1}^{n} primes[i]
const std::array<int, 9> PhiTiny::prime_products = { 1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690 };

// totients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const std::array<int, 9> PhiTiny::totients = { 1, 1, 2, 8, 48, 480, 5760, 92160, 1658880 };

// Number of primes below x
const std::array<int, 23> PhiTiny::pi = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8 };

// Singleton
const PhiTiny phiTiny;

PhiTiny::PhiTiny()
{
  // Initialize phi(x % pp, a) lookup tables
  for (int a = 0; a <= max_a(); a++)
  {
    // For primes <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (primes[a] <= 5)
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
    else
    {
      // For primes > 5 we use a compressed phi(x % pp, a) lookup
      // table. Each bit of the sieve array corresponds to
      // an integer that is not divisible by 2, 3 and 5. Hence
      // the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].

      uint64_t max_x = prime_products[a];
      uint64_t max_x_size = ceil_div(max_x, 240);
      max_x = max_x_size * 240 - 1;
      sieve_[a].resize(max_x_size);

      for (uint64_t i = pi[7]; i <= a; i++)
        for (uint64_t n = primes[i]; n <= max_x; n += primes[i] * 2)
          sieve_[a][n / 240].bits &= unset_bit_[n % 240];

      // Fill an array with the cumulative 1 bit counts.
      // sieve[i][j] contains the count of numbers < j * 240 that
      // are not divisible by any of the first i primes.
      uint64_t count = 0;
      for (auto& sieve : sieve_[a])
      {
        sieve.count = (uint32_t) count;
        count += popcnt64(sieve.bits);
      }
    }
  }
}

} // namespace

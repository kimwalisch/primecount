///
/// @file  PhiTiny.cpp
/// @brief phi(x, a) counts the numbers <= x that are not divisible
///        by any of the first a primes. PhiTiny computes phi(x, a)
///        in constant time for a <= 7 using lookup tables.
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
#include <limits>
#include <vector>

namespace primecount {

const std::array<uint32_t, 8> PhiTiny::primes = { 0, 2, 3, 5, 7, 11, 13, 17 };

// prime_products[n] = \prod_{i=1}^{n} primes[i]
const std::array<uint32_t, 8> PhiTiny::prime_products = { 1, 2, 6, 30, 210, 2310, 30030, 510510 };

// totients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const std::array<uint32_t, 8> PhiTiny::totients = { 1, 1, 2, 8, 48, 480, 5760, 92160 };

// Number of primes <= primes.back()
const std::array<uint8_t, 18> PhiTiny::pi = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7 };

// Singleton
const PhiTiny phiTiny;

PhiTiny::PhiTiny()
{
  // The pi[x] lookup table should contain the
  // number of primes <= primes.back().
  assert(pi.back() == primes.size() - 1);
  assert(primes.back() == pi.size() - 1);

  for (uint64_t a = 0; a < sieve_.size(); a++)
  {
    // For primes <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (a < phi_.size())
    {
      uint64_t pp = prime_products[a];
      phi_[a].resize(pp);
      phi_[a][0] = 0;

      for (uint64_t x = 1; x < pp; x++)
      {
        uint64_t phi_xa = phi(x, a - 1) - phi(x / primes[a], a - 1);
        assert(phi_xa <= std::numeric_limits<uint8_t>::max());
        phi_[a][x] = (uint8_t) phi_xa;
      }
    }
    else
    {
      // For primes > 5 we use a compressed phi(x % pp, a) lookup
      // table. Each bit of the sieve array corresponds to
      // an integer that is not divisible by 2, 3 and 5. Hence
      // the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].

      uint64_t pp = prime_products[a];
      uint64_t size = ceil_div(pp, 240);
      sieve_[a].resize(size);

      for (uint64_t i = pi[7]; i <= a; i++)
        for (uint64_t n = primes[i]; n < pp; n += primes[i] * 2)
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

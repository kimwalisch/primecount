///
/// @file  PhiTiny.cpp
/// @brief phi_tiny(x, a) counts the numbers <= x that are not
///        divisible by any of the first a primes. phi_tiny(x, a)
///        computes phi(x, a) in constant time for a <= 8 using
///        lookup tables and the formula below.
///
///        phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
///        with pp = 2 * 3 * ... * prime[a]
///        φ(pp) = \prod_{i=1}^{a} (prime[i] - 1)
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PhiTiny.hpp>
#include <BitSieve240.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <popcnt.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <algorithm>

namespace primecount {
namespace PhiTiny {

const Array<uint64_t, 8> primes = { 0, 2, 3, 5, 7, 11, 13, 17 };

// prime_products[n] = \prod_{i=1}^{n} primes[i]
const Array<uint64_t, 8> prime_products = { 1, 2, 6, 30, 210, 2310, 30030, 510510 };

// totients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const Array<uint64_t, 8> totients = { 1, 1, 2, 8, 48, 480, 5760, 92160 };

// Number of primes <= next_prime(primes.back())
const Array<uint8_t, 20> pi = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8 };

} // namespace
} // namespace

namespace {

using primecount::PhiTiny::primes;
using primecount::PhiTiny::prime_products;
using primecount::PhiTiny::totients;
using primecount::PhiTiny::pi;
using namespace primecount;

class PhiTinyImpl : public BitSieve240
{
public:
  PhiTinyImpl()
  {
    // The pi[x] lookup table must contain the number
    // of primes <= next_prime(primes.back()).
    ASSERT(pi.back() == primes.size());
    ASSERT(phi_.size() - 1 == (uint64_t) pi[5]);
    ASSERT(sieve_.size() == primes.size());
    static_assert(prime_products.size() == primes.size(), "Invalid prime_products size!");
    static_assert(totients.size() == primes.size(), "Invalid totients size!");

    // a = 0
    phi_[0].resize(1);
    phi_[0][0] = 0;

    for (uint64_t a = 1; a < sieve_.size(); a++)
    {
      // For prime[a] <= 5 our phi(x % pp, a) lookup table
      // is a simple two dimensional array.
      if (a < phi_.size())
      {
        uint64_t pp = prime_products[a];
        phi_[a].resize(pp);
        phi_[a][0] = 0;

        for (uint64_t x = 1; x < pp; x++)
        {
          uint64_t phi_xa = phi(x, a - 1) - phi(x / primes[a], a - 1);
          ASSERT(phi_xa <= pstd::numeric_limits<uint8_t>::max());
          phi_[a][x] = (uint8_t) phi_xa;
        }
      }
      else
      {
        // For prime[a] > 5 we use a compressed phi(x % pp, a)
        // lookup table. Each bit of the sieve array corresponds
        // to an integer that is not divisible by 2, 3 and 5.
        // Hence the 8 bits of each byte correspond to the offsets
        // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
        uint64_t pp = prime_products[a];
        uint64_t size = ceil_div(pp, 240);
        sieve_[a].resize(size);
        std::fill_n(sieve_[a].begin(), size, sieve_t{0, ~0ull});

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

  /// Uses at most one level of phi(x, a) recursion
  /// to ensure that the runtime is O(1).
  template <typename T>
  ALWAYS_INLINE T phi_tiny(T x, uint64_t a) const
  {
    if (a < PhiTiny::max_a())
      return phi(x, a);
    else
    {
      ASSERT(a == 8);
      // This code path will be executed most of the time.
      // In phi7(x) the variable a has been hardcoded to 7
      // which makes it run slightly faster than phi(x, a).
      // phi(x, 8) = phi(x, 7) - phi(x / prime[8], 7)
      return phi7(x) - phi7(x / 19);
    }
  }

private:

  template <typename T>
  ALWAYS_INLINE T phi(T x, uint64_t a) const
  {
    auto pp = prime_products[a];
    auto remainder = (uint64_t)(x % pp);
    T xpp = x / pp;
    T sum = xpp * totients[a];

    // For prime[a] <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (a < phi_.size())
      sum += phi_[a][remainder];
    else
    {
      // For prime[a] > 5 we use a compressed phi(x % pp, a)
      // lookup table. Each bit of the sieve array corresponds
      // to an integer that is not divisible by 2, 3 and 5.
      // Hence the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
      uint64_t count = sieve_[a][remainder / 240].count;
      uint64_t bits = sieve_[a][remainder / 240].bits;
      uint64_t bitmask = unset_larger_[remainder % 240];
      sum += count + popcnt64(bits & bitmask);
    }

    return sum;
  }

  /// In phi7(x) the variable a has been hardcoded to 7.
  /// phi7(x) uses division by a constant instead of regular
  /// integer division and hence phi7(x) is expected to run
  /// faster than the phi(x, a) implementation above.
  ///
  template <typename T>
  ALWAYS_INLINE T phi7(T x) const
  {
    constexpr uint32_t a = 7;
    constexpr uint32_t pp = 510510;
    constexpr uint32_t totient = 92160;
    auto remainder = (uint64_t)(x % pp);
    T xpp = x / pp;
    T sum = xpp * totient;

    // For prime[a] > 5 we use a compressed phi(x % pp, a)
    // lookup table. Each bit of the sieve array corresponds
    // to an integer that is not divisible by 2, 3 and 5.
    // Hence the 8 bits of each byte correspond to the offsets
    // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
    ASSERT(sieve_.size() - 1 == a);
    uint64_t count = sieve_[a][remainder / 240].count;
    uint64_t bits = sieve_[a][remainder / 240].bits;
    uint64_t bitmask = unset_larger_[remainder % 240];
    sum += count + popcnt64(bits & bitmask);

    return sum;
  }

  /// Packing sieve_t increases the cache's capacity by 25%
  /// which improves performance by up to 10%.
  #pragma pack(push, 1)
  struct sieve_t
  {
    uint32_t count;
    uint64_t bits;
  };

  #pragma pack(pop)

  /// sieve[a] contains only numbers that are not divisible
  /// by any of the the first a primes. sieve[a][i].count
  /// contains the count of numbers < i * 240 that are not
  /// divisible by any of the first a primes.
  Array<Vector<sieve_t>, 8> sieve_;
  Array<Vector<uint8_t>, 4> phi_;
};

// Singleton
const PhiTinyImpl phiTinyImpl;

} // namespace

namespace primecount {
namespace PhiTiny {

uint64_t phi_tiny(uint64_t x, uint64_t a)
{
  return phiTinyImpl.phi_tiny(x, a);
}

#if defined(HAVE_INT128_T)

uint128_t phi_tiny(uint128_t x, uint64_t a)
{
  return phiTinyImpl.phi_tiny(x, a);
}

#endif

} // namespace
} // namespace

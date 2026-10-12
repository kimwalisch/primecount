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
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny_lookup_tables.hpp"

#include <PhiTiny.hpp>
#include <BitSieve240.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <popcnt.hpp>
#include <Vector.hpp>

#include <stdint.h>

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

using namespace primecount;

using PhiTiny::primes;
using PhiTiny::prime_products;
using PhiTiny::totients;

static_assert(prime_products.size() == primes.size(), "Invalid prime_products size!");
static_assert(totients.size() == primes.size(), "Invalid totients size!");
static_assert(sieve_tables.size() == primes.size(), "Invalid sieve_tables size!");

class PhiTinyImpl : public BitSieve240
{
public:
  /// Uses at most one level of phi(x, a) recursion
  /// to ensure that the runtime is O(1).
  template <typename T>
  ALWAYS_INLINE static T phi_tiny(T x, uint64_t a)
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
  ALWAYS_INLINE static T phi(T x, uint64_t a)
  {
    auto pp = prime_products[a];
    auto remainder = (uint64_t)(x % pp);
    T xpp = x / pp;
    T sum = xpp * totients[a];

    // For prime[a] <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (a < phi_tables.size())
    {
      ASSERT(remainder < phi_table_sizes[a]);
      sum += phi_tables[a][remainder];
    }
    else
    {
      // For prime[a] > 5 we use a compressed phi(x % pp, a)
      // lookup table. Each bit of the sieve array corresponds
      // to an integer that is not divisible by 2, 3 and 5.
      // Hence the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
      ASSERT(remainder / 240 < sieve_table_sizes[a]);
      const sieve_t& sieve = sieve_tables[a][remainder / 240];
      uint64_t count = sieve.count;
      uint64_t bits = sieve.bits;
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
  ALWAYS_INLINE static T phi7(T x)
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
    ASSERT(sieve_tables.size() - 1 == a);
    ASSERT(remainder / 240 < sieve_table_sizes[a]);
    const sieve_t& sieve = sieve_tables[a][remainder / 240];
    uint64_t count = sieve.count;
    uint64_t bits = sieve.bits;
    uint64_t bitmask = unset_larger_[remainder % 240];
    sum += count + popcnt64(bits & bitmask);

    return sum;
  }
};

} // namespace

namespace primecount {
namespace PhiTiny {

uint64_t phi_tiny(uint64_t x, uint64_t a)
{
  return PhiTinyImpl::phi_tiny(x, a);
}

#if defined(HAVE_INT128_T)

uint128_t phi_tiny(uint128_t x, uint64_t a)
{
  return PhiTinyImpl::phi_tiny(x, a);
}

#endif

} // namespace
} // namespace

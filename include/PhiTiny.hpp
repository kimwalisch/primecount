///
/// @file  PhiTiny.hpp
/// @brief phi_tiny(x, a) counts the numbers <= x that are not
///        divisible by any of the first a primes. phi_tiny(x, a)
///        computes phi(x, a) in constant time for a <= 8 using
///        lookup tables and the formula below.
///
///        phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
///        with pp = 2 * 3 * ... * prime[a]
///        φ(pp) = \prod_{i=1}^{a} (prime[i] - 1)
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <BitSieve240.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <pod_vector.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <limits>
#include <type_traits>

namespace primecount {

class PhiTiny : public BitSieve240
{
public:
  PhiTiny();

  /// Uses at most one level of phi(x, a) recursion
  /// to ensure that the runtime is O(1).
  template <typename T>
  T phi_recursive(T x, uint64_t a) const
  {
    // Unsigned integer division is usually
    // faster than signed integer division,
    // especially for int128_t.
    using UT = typename std::make_unsigned<T>::type;

    if (a < max_a())
      return phi((UT) x, a);
    else
    {
      ASSERT(a == 8);
      // This code path will be executed most of the time.
      // In phi7(x) the variable a has been hardcoded to 7
      // which makes it run slightly faster than phi(x, a).
      // phi(x, 8) = phi(x, 7) - phi(x / prime[8], 7)
      return phi7((UT) x) - phi7((UT) x / 19);
    }
  }

  template <typename T>
  T phi(T x, uint64_t a) const
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
      sum += (T)(count + popcnt64(bits & bitmask));
    }

    return sum;
  }

  /// In phi7(x) the variable a has been hardcoded to 7.
  /// phi7(x) uses division by a constant instead of regular
  /// integer division and hence phi7(x) is expected to run
  /// faster than the phi(x, a) implementation above.
  ///
  template <typename T>
  T phi7(T x) const
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
    sum += (T)(count + popcnt64(bits & bitmask));

    return sum;
  }

  static uint64_t get_c(uint64_t y)
  {
    if (y < pi.size())
      return pi[y];
    else
      return max_a();
  }

  /// In Xavier Gourdon's algorithm the small
  /// constant is named k instead of c.
  /// k <= PrimePi[min(x_star, sqrt(x / y))]
  ///
  template <typename T>
  static uint64_t get_k(T x)
  {
    return get_c(iroot<4>(x));
  }

  static constexpr uint64_t max_a()
  {
    return primes.size();
  }

private:
  static const pod_array<uint32_t, 8> primes;
  static const pod_array<uint32_t, 8> prime_products;
  static const pod_array<uint32_t, 8> totients;
  static const pod_array<uint8_t, 20> pi;

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
  pod_array<pod_vector<sieve_t>, 8> sieve_;
  pod_array<pod_vector<uint8_t>, 4> phi_;
};

extern const PhiTiny phiTiny;

inline bool is_phi_tiny(uint64_t a)
{
  return a <= PhiTiny::max_a();
}

template <typename T>
typename std::enable_if<(sizeof(T) == sizeof(typename make_smaller<T>::type)), T>::type
phi_tiny(T x, uint64_t a)
{
  return phiTiny.phi_recursive(x, a);
}

template <typename T>
typename std::enable_if<(sizeof(T) > sizeof(typename make_smaller<T>::type)), T>::type
phi_tiny(T x, uint64_t a)
{
  using smaller_t = typename make_smaller<T>::type;

  // If possible use smaller integer type
  // to speed up integer division.
  if (x <= std::numeric_limits<smaller_t>::max())
    return phiTiny.phi_recursive((smaller_t) x, a);
  else
    return phiTiny.phi_recursive(x, a);
}

} // namespace

#endif

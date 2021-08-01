///
/// @file  PhiTiny.hpp
/// @brief phi(x, a) counts the numbers <= x that are not divisible
///        by any of the first a primes. PhiTiny computes phi(x, a) in
///        constant time for a <= 8 using lookup tables.
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

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <BitSieve240.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <array>
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

namespace primecount {

class PhiTiny : public BitSieve240
{
public:
  PhiTiny();

  // Uses at most 1 recursion level
  template <typename T>
  T phi_recursive(T x, int64_t a) const
  {
    // Unsigned integer division is usually
    // faster than signed integer division,
    // especially for int128_t.
    using UT = typename std::make_unsigned<T>::type;
    assert(a <= max_a());

    if (a < max_a())
      return phi((UT) x, a);
    else // a == max_a()
      return phi((UT) x, a - 1) - phi((UT) x / primes[a], a - 1);
  }

  template <typename T>
  T phi(T x, uint64_t a) const
  {
    assert(a < prime_products.size());
    auto pp = prime_products[a];
    auto remainder = (uint64_t)(x % pp);
    T xpp = x / pp;
    T sum = xpp * totients[a];

    // For primes <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (a < phi_.size())
      sum += phi_[a][remainder];
    else
    {
      // For primes > 5 we use a compressed phi(x % pp, a) lookup
      // table. Each bit of the sieve array corresponds to
      // an integer that is not divisible by 2, 3 and 5. Hence
      // the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
      assert(a < sieve_.size());
      assert(remainder / 240 < sieve_[a].size());
      uint64_t count = sieve_[a][remainder / 240].count;
      uint64_t bits = sieve_[a][remainder / 240].bits;
      uint64_t bitmask = unset_larger_[remainder % 240];
      sum += (T)(count + popcnt64(bits & bitmask));
    }

    return sum;
  }

  static int64_t get_c(uint64_t y)
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
  static int64_t get_k(T x)
  {
    return get_c(iroot<4>(x));
  }

  static int64_t max_a()
  {
    return primes.size() - 1;
  }

private:
  static const std::array<uint32_t, 9> primes;
  static const std::array<uint32_t, 8> prime_products;
  static const std::array<uint32_t, 8> totients;
  static const std::array<uint8_t, 20> pi;

  /// Packing sieve_t increases the cache's capacity by 25%
  /// which improves performance by up to 10%.
  #pragma pack(push, 1)
  struct sieve_t
  {
    uint32_t count = 0;
    uint64_t bits = ~0ull;
  };

  #pragma pack(pop)

  /// sieve[a] contains only numbers that are not divisible
  /// by any of the the first a primes. sieve[a][i].count
  /// contains the count of numbers < i * 240 that are not
  /// divisible by any of the first a primes.
  std::array<std::vector<sieve_t>, 8> sieve_;
  std::array<std::vector<uint8_t>, 4> phi_;
};

extern const PhiTiny phiTiny;

inline bool is_phi_tiny(int64_t a)
{
  return a <= PhiTiny::max_a();
}

template <typename T>
typename std::enable_if<(sizeof(T) == sizeof(typename make_smaller<T>::type)), T>::type
phi_tiny(T x, int64_t a)
{
  return phiTiny.phi_recursive(x, a);
}

template <typename T>
typename std::enable_if<(sizeof(T) > sizeof(typename make_smaller<T>::type)), T>::type
phi_tiny(T x, int64_t a)
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

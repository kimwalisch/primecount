///
/// @file  PhiTiny.hpp
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

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <fast_div.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <array>
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

namespace primecount {

class PhiTiny
{
public:
  PhiTiny();

  template <typename T>
  T phi(T x, int64_t a) const
  {
    assert(a <= max_a());

    T pp = prime_products[a];
    return (x / pp) * totients[a] + phi_[a][x % pp];
  }

  static int64_t get_c(int64_t y)
  {
    assert(y >= 0);

    if (y >= primes.back())
      return max_a();
    else
      return pi[y];
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
  std::array<std::vector<int16_t>, 7> phi_;
  static const std::array<int, 7> primes;
  static const std::array<int, 7> prime_products;
  static const std::array<int, 7> totients;
  static const std::array<int, 13> pi;
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
  // Unsigned integer division is usually
  // faster than signed integer division.
  using UT = typename std::make_unsigned<T>::type;
  return phiTiny.phi((UT) x, a);
}

template <typename T>
typename std::enable_if<(sizeof(T) > sizeof(typename make_smaller<T>::type)), T>::type
phi_tiny(T x, int64_t a)
{
  using smaller_t = typename make_smaller<T>::type;

  // If possible use smaller integer type
  // to speed up integer division.
  if (x <= std::numeric_limits<smaller_t>::max())
    return phiTiny.phi((smaller_t) x, a);
  else
  {
    // Unsigned integer division is usually
    // faster than signed integer division.
    using UT = typename std::make_unsigned<T>::type;
    return phiTiny.phi((UT) x, a);
  }
}

} // namespace

#endif

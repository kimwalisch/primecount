///
/// @file  PhiTiny.hpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <int128_t.hpp>

#include <stdint.h>
#include <array>
#include <cassert>
#include <vector>

namespace primecount {

class PhiTiny
{
public:
  PhiTiny();

  static int64_t max_a()
  {
    return primes.size() - 1;
  }

  /// Partial sieve function (a.k.a. Legendre-sum).
  /// phi(x, a) counts the numbers <= x that are not divisible
  /// by any of the first a primes.
  /// @pre a <= max_a()
  ///
  template <typename X, typename A>
  X phi(X x, A a) const
  {
    assert(a <= max_a());

    // phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
    // with pp = 2 * 3 * ... * prime[a]
    X pp = prime_products[a];
    return (x / pp) * totients[a] + phi_cache_[a][x % pp];
  }

  static int64_t get_c(int64_t y)
  {
    assert(y >= 0);

    if (y >= primes.back())
      return max_a();
    else
      return pi[y];
  }
private:
  std::array<std::vector<int16_t>, 7> phi_cache_;
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

template <typename X, typename A>
typename prt::make_signed<X>::type phi_tiny(X x, A a)
{
  return phiTiny.phi(x, a);
}

} // namespace

#endif

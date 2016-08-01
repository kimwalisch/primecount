///
/// @file  PhiTiny.hpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <int128.hpp>

#include <stdint.h>
#include <cassert>
#include <vector>

namespace primecount {

class PhiTiny {
public:
  PhiTiny();
  static int64_t max_a() { return 6; }
  static bool is_tiny(int64_t a) { return a <= max_a(); }

  /// Partial sieve function (a.k.a. Legendre-sum).
  /// phi(x, a) counts the numbers <= x that are not divisible
  /// by any of the first a primes.
  /// @pre is_tiny(a).
  ///
  template <typename X, typename A>
  X phi(X x, A a) const
  {
    assert(is_tiny(a));
    // phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a)
    // with pp = 2 * 3 * ... * prime[a]
    X pp = prime_products[a];
    return (x / pp) * totients[a] + phi_cache_[a][x % pp];
  }

  static int64_t get_c(int64_t y)
  {
    assert(y >= 0);

    if (y >= primes[max_a()])
      return max_a();
    else
      return pi[y];
  }
private:
  std::vector<int16_t> phi_cache_[7];
  static const int pi[20];
  static const int primes[7];
  static const int prime_products[7];
  static const int totients[7];
};

inline bool is_phi_tiny(int64_t a)
{
  return PhiTiny::is_tiny(a);
}

#if __cplusplus >= 201103L

template <typename X, typename A>
typename prt::make_signed<X>::type phi_tiny(X x, A a)
{
  extern const PhiTiny phiTiny;
  return phiTiny.phi(x, a);
}

#else /* C++98 */

template <typename X, typename A>
X phi_tiny(X x, A a)
{
  extern const PhiTiny phiTiny;
  return phiTiny.phi(x, a);
}

#endif

} // namespace

#endif /* PHITINY_HPP */

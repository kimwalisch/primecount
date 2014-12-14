///
/// @file  PhiTiny.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <inttypes.hpp>

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
    X x_div_pp = x / pp;
    X x_mod_pp = x - pp * x_div_pp;
    return x_div_pp * totients[a] + phi_cache_[a][x_mod_pp];
  }

  private:
  std::vector<int16_t> phi_cache_[7];
  static const int primes[7];
  static const int prime_products[7];
  static const int totients[7];
};

inline bool is_phi_tiny(int64_t a)
{
  return PhiTiny::is_tiny(a);
}

} // namespace

#if __cplusplus >= 201103L

#include <type_traits>

namespace primecount {

/// Workaround for std::make_signed.
/// std::make_signed fails for int128_t and uint128_t.
///
template <typename T>
struct make_signed_workaround
{
#ifndef HAVE_INT128_T
  typedef typename std::make_signed<T>::type type;
#else
  typedef typename std::conditional<std::is_same<T, uint8_t>::value, int8_t,
          typename std::conditional<std::is_same<T, uint16_t>::value, int16_t,
          typename std::conditional<std::is_same<T, uint32_t>::value, int32_t,
          typename std::conditional<std::is_same<T, uint64_t>::value, int64_t,
          typename std::conditional<std::is_same<T, uint128_t>::value, int128_t,
          T>::type>::type>::type>::type>::type type;
#endif
};

template <typename X, typename A>
typename make_signed_workaround<X>::type phi_tiny(X x, A a)
{
  extern const PhiTiny phiTiny;
  return phiTiny.phi(x, a);
}

} // namespace

#else /* C++98 */

namespace primecount {

template <typename X, typename A>
X phi_tiny(X x, A a)
{
  extern const PhiTiny phiTiny;
  return phiTiny.phi(x, a);
}

} // namespace

#endif

#endif /* PHITINY_HPP */

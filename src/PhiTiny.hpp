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

#include <stdint.h>
#include <vector>
#include <cassert>

namespace primecount {

/// This class calculates phi(x, a) in constant time for small values
/// of a < 7 using lookup tables. Let pp = prime_products_[a]:
/// phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a).
///
class PhiTiny {
public:
  PhiTiny();
  static bool is_cached(int64_t a)
  {
    return a < 7;
  }
  /// Partial sieve function (a.k.a. Legendre-sum).
  /// phi(x, a) counts the numbers <= x that are not divisible
  /// by any of the first a primes.
  /// @pre is_cached(a) == true.
  ///
  int64_t phi(int64_t x, int64_t a) const
  {
    assert(is_cached(a));
    assert(x >= 0);
    return (x / prime_products_[a]) * totients_[a] + phi_cache_[a][x % prime_products_[a]];
  }
private:
  std::vector<int16_t> phi_cache_[7];
  static const int32_t primes_[7];
  static const int32_t prime_products_[7];
  static const int32_t totients_[7];
};

} // namespace primecount

#endif

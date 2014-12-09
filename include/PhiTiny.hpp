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
  template <typename T>
  T phi(T x, int64_t a) const
  {
    assert(x >= 0);
    assert(a <= max_a());
    return (x / prime_products[a]) * totients[a] + phi_cache_[a][x % prime_products[a]];
  }

  private:
  std::vector<int16_t> phi_cache_[7];
  static const int32_t primes[7];
  static const int32_t prime_products[7];
  static const int32_t totients[7];
};

inline bool is_phi_tiny(int64_t a)
{
  return PhiTiny::is_tiny(a);
}

int64_t phi_tiny(int64_t x, int64_t a);

#ifdef HAVE_INT128_T

int128_t phi_tiny(uint128_t x, int64_t a);

#endif

} // namespace primecount

#endif

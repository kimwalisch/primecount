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

namespace primecount {

/// This class calculates phi(x, a) in constant time for small values
/// of a < 7 using lookup tables. Let pp = prime_products_[a]:
/// phi(x, a) = (x / pp) * φ(pp) + phi(x % pp, a).
///
class PhiTiny {
public:
  PhiTiny();
  int64_t phi(int64_t x, int64_t a) const;
  enum { MAX_A = 6 };
private:
  std::vector<int16_t> phi_cache_[7];
  static const int32_t primes_[7];
  static const int32_t prime_products_[7];
  static const int32_t totients_[7];
};

inline bool is_phi_tiny(int64_t a)
{
  return a <= PhiTiny::MAX_A;
}

int64_t phi_tiny(int64_t x, int64_t a);

} // namespace primecount

#endif

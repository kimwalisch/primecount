///
/// @file  PhiTiny.hpp
/// @brief phi(x, a) counts the numbers <= x that are not
///        divisible by any of the first a primes.
///        PhiTiny computes phi(x, a) in constant time
///        for a <= 6 using lookup tables.
///
///        phi(x, a) = (x / pp) * φ(a) + phi(x % pp, a)
///        pp = 2 * 3 * ... * prime[a]
///        φ(a) = \prod_{i=1}^{a} (prime[i] - 1)
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
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
#include <limits>
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
typename prt::make_signed<T>::type phi_tiny(T x, int64_t a)
{
  if (x <= std::numeric_limits<uint32_t>::max())
    return phiTiny.phi((uint32_t) x, a);
  else
    return phiTiny.phi(x, a);
}

} // namespace

#endif

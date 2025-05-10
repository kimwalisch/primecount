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
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PHITINY_HPP
#define PHITINY_HPP

#include <imath.hpp>
#include <int128_t.hpp>
#include <Vector.hpp>

#include <stdint.h>

namespace primecount {
namespace PhiTiny {

extern const Array<uint32_t, 8> primes;
extern const Array<uint32_t, 8> prime_products;
extern const Array<uint32_t, 8> totients;
extern const Array<uint8_t, 20> pi;

extern uint64_t phi_tiny(uint64_t x, uint64_t a);

#if defined(HAVE_INT128_T)
  extern uint128_t phi_tiny(uint128_t x, uint64_t a);
#endif

inline constexpr uint64_t max_a()
{
  return primes.size();
}

inline uint64_t get_c(uint64_t y)
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
uint64_t get_k(T x)
{
  return get_c(iroot<4>(x));
}

} // namespace

inline bool is_phi_tiny(uint64_t a)
{
  return a <= PhiTiny::max_a();
}

template <typename T>
typename std::enable_if<(sizeof(T) <= sizeof(uint64_t)), T>::type
phi_tiny(T x, uint64_t a)
{
  return PhiTiny::phi_tiny(uint64_t(x), a);
}

#if defined(HAVE_INT128_T)

template <typename T>
typename std::enable_if<(sizeof(T) >= sizeof(uint128_t)), T>::type
phi_tiny(T x, uint64_t a)
{
  // If possible use smaller integer type
  // to speed up integer division.
  if (x <= pstd::numeric_limits<uint64_t>::max())
    return PhiTiny::phi_tiny(uint64_t(x), a);
  else
    return PhiTiny::phi_tiny(uint128_t(x), a);
}

#endif

} // namespace

#endif

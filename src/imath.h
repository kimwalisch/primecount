///
/// @file  imath.h
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef IMATH_PRIMECOUNT_H
#define IMATH_PRIMECOUNT_H

#include <stdint.h>
#include <cmath>

namespace primecount {

inline int64_t isquare(int64_t x)
{
  return x * x;
}

/// Raise to power using template meta-programming
template <int N>
struct ipow_helper
{
  static int64_t ipow(int64_t x)
  {
    return x * ipow_helper<N - 1>::ipow(x);
  }
};

/// Raise to power using template meta-programming
template <>
struct ipow_helper<0>
{
  static int64_t ipow(int64_t)
  {
    return 1;
  }
};

/// Raise to power using template meta-programming
template <int N>
inline int64_t ipow(int64_t x)
{
  return ipow_helper<N>::ipow(x);
}

/// Integer suare root
inline int32_t isqrt(int64_t x)
{
  int32_t r = static_cast<int32_t>(std::sqrt(static_cast<double>(x)));
  // correct rounding error
  while (ipow<2>(r) > x)
    r--;
  while (ipow<2>(r + 1) <= x)
    r++;
  return r;
}

/// Integer nth root
template <int N>
inline int32_t iroot(int64_t x)
{
  int32_t r = static_cast<int32_t>(std::pow(static_cast<double>(x), 1.0 / N));
  // correct rounding error
  while (ipow<N>(r) > x)
    r--;
  while (ipow<N>(r + 1) <= x)
    r++;
  return r;
}

} // namespace primecount

#endif

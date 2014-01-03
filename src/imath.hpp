///
/// @file  imath.hpp
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
#include <vector>
#include <limits>

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

/// Initialize a vector with Möbius function values.
/// This implementation is based on code by Rick Sladkey
/// posted here: http://mathoverflow.net/a/99545
///
inline void init_moebius(std::vector<int32_t>& mu, int64_t max)
{
  mu.resize(max + 1, 1);

  for (int32_t i = 2; i * i <= max; i++)
  {
    if (mu[i] == 1)
    {
      for (int32_t j = i; j <= max; j += i)
        mu[j] *= -i;
      for (int32_t j = i * i; j <= max; j += i * i)
        mu[j] = 0;
    }
  }
  for (int32_t i = 2; i <= max; i++)
  {
    if (mu[i] == i)
      mu[i] = 1;
    else if (mu[i] == -i)
      mu[i] = -1;
    else if (mu[i] < 0)
      mu[i] = 1;
    else if (mu[i] > 0)
      mu[i] = -1;
  }
}

/// Initialize a vector with the least prime
/// factors of the integers <= max.
///
inline void init_least_factor(std::vector<int32_t>& lpf, int64_t max)
{
  lpf.resize(max + 1, 1);

  // phi(x / 1, c) contributes to the sum,
  // lpf[1] = MAX in order to pass
  // if (least_factor[1] > primes[c])
  if (lpf.size() > 1)
    lpf[1] = std::numeric_limits<int32_t>::max();

  for (int32_t i = 2; i * i <= max; i++)
    if (lpf[i] == 1)
      for (int32_t j = i * 2; j <= max; j += i)
        if (lpf[j] == 1)
          lpf[j] = i;

  for (int32_t i = 2; i <= max; i++)
    if (lpf[i] == 1)
      lpf[i] = i;
}

} // namespace primecount

#endif

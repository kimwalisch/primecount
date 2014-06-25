///
/// @file  pmath.hpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PMATH_HPP
#define PMATH_HPP

#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

namespace primecount {

/// Generate a vector with Möbius function values.
std::vector<int32_t> make_moebius(int64_t max);

/// Generate a vector with the least prime
/// factors of the integers <= max.
///
std::vector<int32_t> make_least_prime_factor(int64_t max);

/// Generate a vector with the prime counts below max
/// using the sieve of Eratosthenes.
///
std::vector<int32_t> make_pi(int64_t max);

inline int64_t isquare(int64_t x)
{
  return x * x;
}

template <typename T1, typename T2, typename T3>
inline T2 in_between(T1 min, T2 x, T3 max)
{
  if (x < min)
    return static_cast<T2>(min);
  if (x > max)
    return static_cast<T2>(max);

  return x;
}

template <typename T>
inline T number_of_bits(T)
{
  return static_cast<T>(sizeof(T) * 8);
}

/// @brief  Check if an integer is a power of 2.
/// @see    Book "Hacker's Delight".
///
template <typename T>
inline bool is_power_of_2(T x)
{
  return (x != 0 && (x & (x - 1)) == 0);
}

/// @brief  Round up to the next power of 2.
/// @see    Book "Hacker's Delight".
///
template <typename T>
inline T next_power_of_2(T x)
{
  if (x == 0)
    return 1;

  x--;
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return ++x;
}

template <typename T>
inline int ilog(T x)
{
  return static_cast<int>(std::log(static_cast<double>(x)));
}

/// Raise to power
template <int N>
inline int64_t ipow(int64_t x)
{
  int64_t r = 1;
  for (int i = 0; i < N; i++)
    r *= x;

  return r;
}

/// Integer suare root
template <typename T>
inline T isqrt(T x)
{
  T r = static_cast<T>(std::sqrt(static_cast<double>(x)));
  // correct rounding error
  while (isquare(r) > x)
    r--;
  while (isquare(r + 1) <= x)
    r++;
  return r;
}

/// Integer nth root
template <int N, typename T>
inline T iroot(T x)
{
  T r = static_cast<T>(std::pow(static_cast<double>(x), 1.0 / N));
  // correct rounding error
  while (ipow<N>(r) > x)
    r--;
  while (ipow<N>(r + 1) <= x)
    r++;
  return r;
}

/// Calculate the number of primes below x using binary search.
/// @pre primes[1] = 2, primes[3] = 3, ...
/// @pre x <= primes.back()
///
template <typename T1, typename T2>
inline T2 pi_bsearch(const std::vector<T1>& primes, T2 x)
{
  // primecount uses 1-indexing
  assert(primes[0] == 0);
  return static_cast<T2>(std::upper_bound(primes.begin() + 1, primes.end(), x) - (primes.begin() + 1));
}

/// Calculate the number of primes below x using binary search.
/// @pre primes[1] = 2, primes[3] = 3, ...
/// @pre x <= primes.back()
///
template <typename T1, typename T2, typename T3>
inline T3 pi_bsearch(const std::vector<T1>& primes, T2 len, T3 x)
{
  // primecount uses 1-indexing
  assert(primes[0] == 0);
  return static_cast<T3>(std::upper_bound(primes.begin() + 1, primes.begin() + len + 1, x) - (primes.begin() + 1));
}

} // namespace

#endif

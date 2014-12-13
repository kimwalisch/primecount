///
/// @file  Li.cpp
/// @brief Logarithmic integral approximation.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <inttypes.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

namespace {

/// Calculate the logarithmic integral using Ramanujan's fast
/// converging formula (accurate up to 10^17).
/// @see http://mathworld.wolfram.com/LogarithmicIntegral.html (15)
///
long double li(long double x)
{
  long double GAMMA = 0.57721566490153286061;
  long double logx = log(x);
  long double sum = 0;
  long double inner_sum = 0;
  long double factorial = 1;
  long double p = -1;
  long double power2 = 1;

  int k = 0;
  int terms = (int) (logx * 2 + 10);

  for (int n = 1; n < terms; n++)
  {
    factorial *= n;
    p *= -logx;
    long double q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0 / (2 * k + 1);
    sum += (p / q) * inner_sum;
  }

  return GAMMA + log(logx) + sqrt(x) * sum;
}

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes below x.
/// @post Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
template <typename T>
T Li(T x)
{
  if (x < 2)
    return 0;

  long double n = (long double) x;
  return (T) (li(n) - /* li(2) = */ 1.04516);
}

/// Calculate the inverse logarithmic integral Li^-1(x) which is a
/// very accurate approximation of the nth prime.
/// @post Li_inverse(x) < nth_prime(x) for 7 <= x <= 10^316
///
template <typename T>
T Li_inverse(T x)
{
  if (x < 1)
    return 0;

  long double logx = max(1.0, log((double) x));
  T first = (T) (x * logx);
  T last  = (T) (x * logx * 2 + 2);

  // overflow
  if (first < x) return -1;
  if (last  < x) last = numeric_limits<T>::max();

  // Find Li^-1(x) using binary search
  while (first < last)
  {
    T mid = first + (last - first) / 2;
    if (Li(mid) < x)
      first = mid + 1;
    else
      last = mid;
  }

  return first;
}

} // namespace

namespace primecount {

int64_t Li(int64_t x)
{
  return ::Li(x);
}

int64_t Li_inverse(int64_t x)
{
  return ::Li_inverse(x);
}

#ifdef HAVE_INT128_T

int128_t Li(int128_t x)
{
  return ::Li(x);
}

int128_t Li_inverse(int128_t x)
{
  return ::Li_inverse(x);
}

#endif

} // namespace primecount

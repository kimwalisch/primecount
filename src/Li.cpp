///
/// @file  Li.cpp
/// @brief Logarithmic integral approximation.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cmath>
#include <limits>

using namespace std;
using namespace primecount;

namespace {
namespace Li {

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation
///
long double li(long double x)
{
  long double gamma = 0.57721566490153286061;
  long double sum = 0;
  long double inner_sum = 0;
  long double prev_sum = 1;
  long double factorial = 1;
  long double p = -1;
  long double q = 0;
  long double power2 = 1;

  int k = 0;
  int n = 1;

  while (abs(sum - prev_sum) > numeric_limits<double>::epsilon())
  {
    p *= -log(x);
    factorial *= n;
    q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0L / (2 * k + 1);
    prev_sum = sum;
    sum += (p / q) * inner_sum;
    n++;
  }

  return gamma + log(log(x)) + sqrt(x) * sum;
}

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes below x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
template <typename T>
T Li(T x)
{
  if (x < 2)
    return 0;

  long double n = (long double) x;
  return (T) (li(n) - /* li(2) = */ 1.04516);
}

/// Calculate the inverse logarithmic integral Li^-1(x) which
/// is a very accurate approximation of the nth prime.
/// Li^-1(x) < nth_prime(x) for 7 <= x <= 10^316
///
template <typename T>
T Li_inverse(T x)
{
  if (x < 1)
    return 0;

  T first = 1;
  T last = prt::numeric_limits<T>::max();

  // Find using binary search
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

} // namespace Li
} // namespace primecount

namespace primecount {

int64_t Li(int64_t x)
{
  return Li::Li(x);
}

int64_t Li_inverse(int64_t x)
{
  return Li::Li_inverse(x);
}

#ifdef HAVE_INT128_T

int128_t Li(int128_t x)
{
  return Li::Li(x);
}

int128_t Li_inverse(int128_t x)
{
  return Li::Li_inverse(x);
}

#endif

} // namespace

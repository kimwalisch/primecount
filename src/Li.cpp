///
/// @file  Li.cpp
/// @brief Logarithmic integral and Riemann R functions.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <generate.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

namespace primecount {

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation
///
long double li(long double x)
{
  long double gamma = 0.57721566490153286061;
  long double sum = 0;
  long double inner_sum = 0;
  long double factorial = 1;
  long double p = -1;
  long double q = 0;
  long double power2 = 1;
  int k = 0;

  if (x < 1.001)
    return 0;

  for (int n = 1; n < 200; n++)
  {
    p *= -log(x);
    factorial *= n;
    q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0L / (2 * k + 1);
    long double old = sum;
    sum += (p / q) * inner_sum;
    if (abs(sum - old) < numeric_limits<double>::epsilon())
      break;
  }

  return gamma + log(log(x)) + sqrt(x) * sum;
}

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
long double Li(long double x)
{
  long double li2 = 1.04516378011749278484;

  if (x < 1.001)
    return 0;

  return li(x) - li2;
}

/// Calculate the inverse offset logarithmic integral which
/// is a very accurate approximation of the nth prime.
/// Li^-1(x) < nth_prime(x) for 7 <= x <= 10^316
///
long double Li_inverse(long double x)
{
  long double t = x * log(x);

  if (x < 2)
    return 0;

  for (int i = 0; i < 100; i++)
  {
    long double old = t;
    t -= (Li(t) - x) * log(t);
    if (abs(t - old) < numeric_limits<double>::epsilon() * max(abs(t), abs(old)))
      break;
  }

  return t;
}

/// Calculate the Riemann R function which is a very accurate
/// approximation of the number of primes below x.
///
long double Ri(long double x)
{
  int terms = 200;
  auto mu = generate_moebius(terms);
  long double sum = 0;

  for (int n = 1; n < terms; n++)
  {
    long double root = pow(x, 1.0L / n);
    long double Li_root = Li(root);
    sum += (Li_root * mu[n]) / n;
    if (abs(Li_root) < numeric_limits<double>::epsilon())
      break;
  }

  return sum;
}

/// Calculate the inverse Riemann R function which is a very
/// accurate approximation of the nth prime.
///
long double Ri_inverse(long double x)
{
  long double t = x * log(x);

  if (x < 2)
    return 0;

  for (int i = 0; i < 100; i++)
  {
    long double old = t;
    t -= (Ri(t) - x) * log(t);
    if (abs(t - old) < numeric_limits<double>::epsilon() * max(abs(t), abs(old)))
      break;
  }

  return t;
}

int64_t Li(int64_t x)
{
  return (int64_t) Li((long double) x);
}

int64_t Li_inverse(int64_t x)
{
  return (int64_t) Li_inverse((long double) x);
}

#ifdef HAVE_INT128_T

int128_t Li(int128_t x)
{
  return (int128_t) Li((long double) x);
}

int128_t Li_inverse(int128_t x)
{
  return (int128_t) Li_inverse((long double) x);
}

#endif

int64_t Ri(int64_t x)
{
  return (int64_t) Ri((long double) x);
}

int64_t Ri_inverse(int64_t x)
{
  return (int64_t) Ri_inverse((long double) x);
}

#ifdef HAVE_INT128_T

int128_t Ri(int128_t x)
{
  return (int128_t) Ri((long double) x);
}

int128_t Ri_inverse(int128_t x)
{
  return (int128_t) Ri_inverse((long double) x);
}

#endif

} // namespace

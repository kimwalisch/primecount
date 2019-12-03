///
/// @file  Li.cpp
/// @brief Logarithmic integral and Riemann R function.
///        Both the Logarithmic integral and the Riemann R function
///        are very accurate approximations of PrimePi(x). The inverse
///        Logarithmic integral and the inverse Riemann R function are
///        very accurate approximations of the nth prime.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <generate.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cassert>
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
  assert(x >= 2);

  long double gamma = 0.57721566490153286061L;
  long double sum = 0;
  long double inner_sum = 0;
  long double factorial = 1;
  long double p = -1;
  long double q = 0;
  long double power2 = 1;
  long double logx = log(x);
  int k = 0;

  for (int n = 1; n < 200; n++)
  {
    p *= -logx;
    factorial *= n;
    q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0L / (2 * k + 1);
    long double term = (p / q) * inner_sum;
    sum += term;
    if (abs(term) < numeric_limits<long double>::epsilon())
      break;
  }

  return gamma + log(logx) + sqrt(x) * sum;
}

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
long double Li(long double x)
{
  if (x < 2)
    return 0;

  long double li2 = 1.04516378011749278484L;

  return li(x) - li2;
}

/// Calculate the inverse offset logarithmic integral which
/// is a very accurate approximation of the nth prime.
/// Li^-1(x) < nth_prime(x) for 7 <= x <= 10^316
///
long double Li_inverse(long double x)
{
  if (x < 2)
    return 0;

  long double t = x * log(x);
  long double old_term = numeric_limits<long double>::infinity();

  for (int i = 0; i < 100; i++)
  {
    long double term = (Li(t) - x) * log(t);
    t -= term;
    // not converging anymore
    if (abs(term) >= abs(old_term))
      break;
    old_term = term;
  }

  return t;
}

/// Calculate the Riemann R function which is a very accurate
/// approximation of the number of primes below x.
///
long double Ri(long double x)
{
  if (x < 2)
    return 0;

  int terms = 200;
  auto mu = generate_moebius(terms);
  long double sum = 0;

  for (int n = 1; n < terms; n++)
  {
    if (mu[n])
    {
      long double root = pow(x, 1.0L / n);
      long double term = (Li(root) * mu[n]) / n;
      sum += term;
      if (abs(term) < numeric_limits<long double>::epsilon())
        break;
    }
  }

  return sum;
}

/// Calculate the inverse Riemann R function which is a very
/// accurate approximation of the nth prime.
///
long double Ri_inverse(long double x)
{
  if (x < 2)
    return 0;

  long double t = Li_inverse(x);
  long double old_term = numeric_limits<long double>::infinity();

  for (int i = 0; i < 100; i++)
  {
    long double term = (Ri(t) - x) * log(t);
    t -= term;
    // not converging anymore
    if (abs(term) >= abs(old_term))
      break;
    old_term = term;
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

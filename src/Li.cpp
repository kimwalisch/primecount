///
/// @file  Li.cpp
/// @brief Logarithmic integral and Riemann R function.
///        Both the Logarithmic integral and the Riemann R function
///        are very accurate approximations of PrimePi(x). The inverse
///        Logarithmic integral and the inverse Riemann R function are
///        very accurate approximations of the nth prime.
///
///        Note that these implementations are only accurate up to
///        about 10^19 (if the compiler supports the long double type).
///        We also include implementations based on libquadmath and the
///        non standard __float128 type. However it is currently
///        impossible to compile this code without warnings
///        (using -Wpedantic) because of a GCC bug and hence this code
///        is disabled for the time being.
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

#if defined(HAVE_QUADMATH)
  #include <quadmath.h>
#endif

using namespace std;

namespace primecount {

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation
///
long double li(long double x)
{
  assert(x >= 2);

  long double gamma = 0.577215664901532860606512090082402431L;
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

  long double li2 = 1.045163780117492784844588889194613136L;

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
    // Not converging anymore
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
    // Not converging anymore
    if (abs(term) >= abs(old_term))
      break;
    old_term = term;
  }

  return t;
}

#if defined(HAVE_QUADMATH)

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation
///
__float128 li(__float128 x)
{
  assert(x >= 2);

  __float128 gamma = 0.577215664901532860606512090082402431Q;
  __float128 sum = 0;
  __float128 inner_sum = 0;
  __float128 factorial = 1;
  __float128 p = -1;
  __float128 q = 0;
  __float128 power2 = 1;
  __float128 logx = logq(x);
  int k = 0;

  for (int n = 1; n < 300; n++)
  {
    p *= -logx;
    factorial *= n;
    q = factorial * power2;
    power2 *= 2;
    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0Q / (2 * k + 1);
    __float128 term = (p / q) * inner_sum;
    sum += term;
    if (fabsq(term) < FLT128_EPSILON)
      break;
  }

  return gamma + logq(logx) + sqrtq(x) * sum;
}

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
__float128 Li(__float128 x)
{
  if (x < 2)
    return 0;

  __float128 li2 = 1.045163780117492784844588889194613136Q;

  return li(x) - li2;
}

/// Calculate the inverse offset logarithmic integral which
/// is a very accurate approximation of the nth prime.
/// Li^-1(x) < nth_prime(x) for 7 <= x <= 10^316
///
__float128 Li_inverse(__float128 x)
{
  if (x < 2)
    return 0;

  __float128 t = x * logq(x);
  __float128 old_term = FLT128_MAX;

  for (int i = 0; i < 200; i++)
  {
    __float128 term = (Li(t) - x) * logq(t);
    t -= term;
    // Not converging anymore
    if (fabsq(term) >= fabsq(old_term))
      break;
    old_term = term;
  }

  return t;
}

/// Calculate the Riemann R function which is a very accurate
/// approximation of the number of primes below x.
///
__float128 Ri(__float128 x)
{
  if (x < 2)
    return 0;

  int terms = 300;
  auto mu = generate_moebius(terms);
  __float128 sum = 0;

  for (int n = 1; n < terms; n++)
  {
    if (mu[n])
    {
      __float128 root = powq(x, 1.0Q / n);
      __float128 term = (Li(root) * mu[n]) / n;
      sum += term;
      if (fabsq(term) < FLT128_EPSILON)
        break;
    }
  }

  return sum;
}

/// Calculate the inverse Riemann R function which is a very
/// accurate approximation of the nth prime.
///
__float128 Ri_inverse(__float128 x)
{
  if (x < 2)
    return 0;

  __float128 t = Li_inverse(x);
  __float128 old_term = FLT128_MAX;

  for (int i = 0; i < 200; i++)
  {
    __float128 term = (Ri(t) - x) * logq(t);
    t -= term;
    // Not converging anymore
    if (fabsq(term) >= fabsq(old_term))
      break;
    old_term = term;
  }

  return t;
}

#endif

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

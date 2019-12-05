///
/// @file  Li.cpp
/// @brief Logarithmic integral and Riemann R function.
///        Both the Logarithmic integral and the Riemann R function
///        are very accurate approximations of PrimePi(x). The inverse
///        Logarithmic integral and the inverse Riemann R function are
///        very accurate approximations of the nth prime.
///
///        Note that these implementations are only accurate up to
///        about 10^19 (if the long double type has 80 bits). We also
///        include implementations based on the non standard __float128
///        type and libquadmath that can be enabled using
///        'cmake -DWITH_FLOAT128=ON'. Currently __float128 support
///        is disabled by default because there are warnings during
///        compilation (using -Wpedantic) although the code
///        works perfectly fine.
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

#if defined(HAVE_FLOAT128)
  #include <quadmath.h>
#endif

namespace {

using namespace std;
using namespace primecount;

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

  for (int n = 1; true; n++)
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

  while (true)
  {
    long double term = (Li(t) - x) * log(t);

    // Not converging anymore
    if (abs(term) >= abs(old_term))
      break;

    t -= term;
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

  double terms = log2(x) + 10;
  auto mu = generate_moebius((int) terms);
  long double sum = 0;

  for (size_t n = 1; true; n++)
  {
    if (n > mu.size())
      mu = generate_moebius(n * 2);

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

  while (true)
  {
    long double term = (Ri(t) - x) * log(t);

    // Not converging anymore
    if (abs(term) >= abs(old_term))
      break;

    t -= term;
    old_term = term;
  }

  return t;
}

#if defined(HAVE_FLOAT128)

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

  for (int n = 1; true; n++)
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

  while (true)
  {
    __float128 term = (Li(t) - x) * logq(t);

    // Not converging anymore
    if (fabsq(term) >= fabsq(old_term))
      break;

    t -= term;
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

  __float128 terms = log2q(x) + 10;
  auto mu = generate_moebius((int) terms);
  __float128 sum = 0;

  for (size_t n = 1; true; n++)
  {
    if (n > mu.size())
      mu = generate_moebius(n * 2);

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

  while (true)
  {
    __float128 term = (Ri(t) - x) * logq(t);

    // Not converging anymore
    if (fabsq(term) >= fabsq(old_term))
      break;

    t -= term;
    old_term = term;
  }

  return t;
}

#endif

} // namespace

namespace primecount {

int64_t Li(int64_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int64_t) ::Li((__float128) x);
#endif

  return (int64_t) ::Li((long double) x);
}

int64_t Li_inverse(int64_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int64_t) ::Li_inverse((__float128) x);
#endif

  return (int64_t) ::Li_inverse((long double) x);
}

int64_t Ri(int64_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int64_t) ::Ri((__float128) x);
#endif

  return (int64_t) ::Ri((long double) x);
}

int64_t Ri_inverse(int64_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int64_t) ::Ri_inverse((__float128) x);
#endif

  return (int64_t) ::Ri_inverse((long double) x);
}

#ifdef HAVE_INT128_T

int128_t Li(int128_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int128_t) ::Li((__float128) x);
#endif

  return (int128_t) ::Li((long double) x);
}

int128_t Li_inverse(int128_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int128_t) ::Li_inverse((__float128) x);
#endif

  return (int128_t) ::Li_inverse((long double) x);
}

int128_t Ri(int128_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int128_t) ::Ri((__float128) x);
#endif

  return (int128_t) ::Ri((long double) x);
}

int128_t Ri_inverse(int128_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int128_t) ::Ri_inverse((__float128) x);
#endif

  return (int128_t) ::Ri_inverse((long double) x);
}

#endif

} // namespace

///
/// @file  RiemannR.cpp
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
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <generate.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cmath>
#include <limits>

#if defined(HAVE_FLOAT128)
  #include <quadmath.h>
#endif

namespace {

using namespace primecount;

/// Calculate the logarithmic integral using
/// Ramanujan's formula:
/// https://en.wikipedia.org/wiki/Logarithmic_integral_function#Series_representation
///
long double li(long double x)
{
  if (x <= 1)
    return 0;

  long double gamma = 0.577215664901532860606512090082402431L;
  long double sum = 0;
  long double inner_sum = 0;
  long double factorial = 1;
  long double p = -1;
  long double q = 0;
  long double power2 = 1;
  long double logx = std::log(x);
  int k = 0;

  for (int n = 1; true; n++)
  {
    p *= -logx;
    factorial *= n;
    q = factorial * power2;
    power2 *= 2;

    for (; k <= (n - 1) / 2; k++)
      inner_sum += 1.0L / (2 * k + 1);

    auto old_sum = sum;
    sum += (p / q) * inner_sum;

    // Not converging anymore
    if (std::abs(sum - old_sum) < std::numeric_limits<long double>::epsilon())
      break;
  }

  return gamma + std::log(logx) + std::sqrt(x) * sum;
}

/// Calculate the offset logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
long double Li(long double x)
{
  long double li2 = 1.045163780117492784844588889194613136L;

  if (x <= li2)
    return 0;
  else
    return li(x) - li2;
}

/// Calculate the inverse offset logarithmic integral which
/// is a very accurate approximation of the nth prime.
/// Li^-1(x) < nth_prime(x) for 7 <= x <= 10^316
///
/// This implementation computes Li^-1(x) as the zero of the
/// function f(z) = Li(z) - x using the Newton–Raphson method.
/// Note that Li'(z) = 1 / log(z).
/// https://math.stackexchange.com/a/853192
///
/// Newton–Raphson method:
/// zn+1 = zn - (f(zn) / f'(zn)).
/// zn+1 = zn - (Li(zn) - x) / (1 / log(zn))
/// zn+1 = zn - (Li(zn) - x) * log(zn)
///
long double Li_inverse(long double x)
{
  if (x < 2)
    return 0;

  long double t = x * std::log(x);
  long double old_term = std::numeric_limits<long double>::infinity();

  while (true)
  {
    long double term = (Li(t) - x) * std::log(t);

    // Not converging anymore
    if (std::abs(term) >= std::abs(old_term))
      break;

    t -= term;
    old_term = term;
  }

  return t;
}

/// Calculate the Riemann R function which is a very accurate
/// approximation of the number of primes below x.
/// RiemannR(x) = \sum_{n=1}^{∞} μ(n)/n * li(x^(1/n))
/// http://mathworld.wolfram.com/RiemannPrimeCountingFunction.html
///
long double Ri(long double x)
{
  if (x <= 1)
    return 0;

  long double sum = 0;
  long double old_term = std::numeric_limits<long double>::infinity();
  auto terms = (int) (std::log2(x) * 2 + 10);
  auto mu = generate_moebius(terms);

  for (int n = 1; n < terms; n++)
  {
    if (mu[n])
    {
      long double root = std::pow(x, 1.0L / n);
      long double term = (li(root) * mu[n]) / n;

      // Not converging anymore
      if (std::abs(term) >= std::abs(old_term))
        break;

      sum += term;
      old_term = term;
    }
  }

  return sum;
}

/// Calculate the inverse Riemann R function which is a very
/// accurate approximation of the nth prime.
/// This implementation computes Ri^-1(x) as the zero of the
/// function f(z) = Ri(z) - x using the Newton–Raphson method.
/// Note that Ri'(z) = 1 / log(z).
/// https://math.stackexchange.com/a/853192
///
/// Newton–Raphson method:
/// zn+1 = zn - (f(zn) / f'(zn)).
/// zn+1 = zn - (Ri(zn) - x) / (1 / log(zn))
/// zn+1 = zn - (Ri(zn) - x) * log(zn)
///
long double Ri_inverse(long double x)
{
  if (x < 2)
    return 0;

  long double t = Li_inverse(x);
  long double old_term = std::numeric_limits<long double>::infinity();

  while (true)
  {
    long double term = (Ri(t) - x) * std::log(t);

    // Not converging anymore
    if (std::abs(term) >= std::abs(old_term))
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
  if (x <= 1)
    return 0;

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

    auto old_sum = sum;
    sum += (p / q) * inner_sum;

    // Not converging anymore
    if (fabsq(sum - old_sum) < FLT128_EPSILON)
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
  __float128 li2 = 1.045163780117492784844588889194613136Q;

  if (x <= li2)
    return 0;
  else
    return li(x) - li2;
}

/// Calculate the inverse offset logarithmic integral which
/// is a very accurate approximation of the nth prime.
/// Li^-1(x) < nth_prime(x) for 7 <= x <= 10^316
///
/// This implementation computes Li^-1(x) as the zero of the
/// function f(z) = Li(z) - x using the Newton–Raphson method.
/// Note that Li'(z) = 1 / log(z).
/// https://math.stackexchange.com/a/853192
///
/// Newton–Raphson method:
/// zn+1 = zn - (f(zn) / f'(zn)).
/// zn+1 = zn - (Li(zn) - x) / (1 / log(zn))
/// zn+1 = zn - (Li(zn) - x) * log(zn)
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
/// RiemannR(x) = \sum_{n=1}^{∞} μ(n)/n * li(x^(1/n))
/// http://mathworld.wolfram.com/RiemannPrimeCountingFunction.html
///
__float128 Ri(__float128 x)
{
  if (x <= 1)
    return 0;

  __float128 sum = 0;
  __float128 old_term = FLT128_MAX;
  auto terms = (int) (log2q(x) * 2 + 10);
  auto mu = generate_moebius(terms);

  for (int n = 1; n < terms; n++)
  {
    if (mu[n])
    {
      __float128 root = powq(x, 1.0Q / n);
      __float128 term = (li(root) * mu[n]) / n;

      // Not converging anymore
      if (fabsq(term) >= fabsq(old_term))
        break;

      sum += term;
      old_term = term;
    }
  }

  return sum;
}

/// Calculate the inverse Riemann R function which is a very
/// accurate approximation of the nth prime.
/// This implementation computes Ri^-1(x) as the zero of the
/// function f(z) = Ri(z) - x using the Newton–Raphson method.
/// Note that Ri'(z) = 1 / log(z).
/// https://math.stackexchange.com/a/853192
///
/// Newton–Raphson method:
/// zn+1 = zn - (f(zn) / f'(zn)).
/// zn+1 = zn - (Ri(zn) - x) / (1 / log(zn))
/// zn+1 = zn - (Ri(zn) - x) * log(zn)
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

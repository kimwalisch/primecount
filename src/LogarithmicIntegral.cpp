///
/// @file  LogarithmicIntegral.cpp
/// @brief The Logarithmic integral function is a very accurate
///        approximation of PrimePi(x). The inverse Logarithmic
///        integral function is a very accurate approximation of the
///        nth prime.
///
///        Note that these implementations are only accurate up to
///        10^18 (if the long double type has 80 bits). We also
///        include implementations based on the non standard
///        __float128 type and libquadmath that can be enabled
///        using 'cmake -DWITH_FLOAT128=ON'. Currently __float128
///        support is disabled by default because there are warnings
///        during compilation (using -Wpedantic) although the code
///        works perfectly fine.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cmath>
#include <limits>

#if defined(HAVE_FLOAT128)
  #include <quadmath.h>
#endif

namespace {

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

/// Calculate the Eulerian logarithmic integral which is a very
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

/// Calculate the inverse Eulerian logarithmic integral which
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

/// Calculate the Eulerian logarithmic integral which is a very
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

/// Calculate the inverse Eulerian logarithmic integral which
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

#endif

} // namespace

namespace primecount {

int64_t Li(int64_t x)
{
#if defined(HAVE_FLOAT128)
  double long_double_mantissa_bits = std::numeric_limits<long double>::digits;

  if (x > 1e10 && std::log2(x) >= long_double_mantissa_bits / 1.35)
    return (int64_t) ::Li((__float128) x);
#endif
  return (int64_t) ::Li((long double) x);
}

int64_t Li_inverse(int64_t x)
{
#if defined(HAVE_FLOAT128)
  double long_double_mantissa_bits = std::numeric_limits<long double>::digits;

  if (x > 1e10 && std::log2(x) >= long_double_mantissa_bits / 1.35)
  {
    __float128 res = ::Li_inverse((__float128) x);
    if (res > (__float128) std::numeric_limits<int64_t>::max())
      return std::numeric_limits<int64_t>::max();
    else
      return (int64_t) res;
  }
#endif

  long double res = ::Li_inverse((long double) x);
  if (res > (long double) std::numeric_limits<int64_t>::max())
    return std::numeric_limits<int64_t>::max();
  else
    return (int64_t) res;
}

#ifdef HAVE_INT128_T

int128_t Li(int128_t x)
{
#if defined(HAVE_FLOAT128)
  double long_double_mantissa_bits = std::numeric_limits<long double>::digits;

  if (x > 1e10 && std::log2(x) >= long_double_mantissa_bits / 1.35)
    return (int128_t) ::Li((__float128) x);
#endif
  return (int128_t) ::Li((long double) x);
}

int128_t Li_inverse(int128_t x)
{
#if defined(HAVE_FLOAT128)
  double long_double_mantissa_bits = std::numeric_limits<long double>::digits;

  if (x > 1e10 && std::log2(x) >= long_double_mantissa_bits / 1.35)
  {
    __float128 res = ::Li_inverse((__float128) x);
    if (res > (__float128) std::numeric_limits<int128_t>::max())
      return std::numeric_limits<int128_t>::max();
    else
      return (int128_t) res;
  }
#endif

  long double res = ::Li_inverse((long double) x);
  if (res > (long double) std::numeric_limits<int128_t>::max())
    return std::numeric_limits<int128_t>::max();
  else
    return (int128_t) res;
}

#endif

} // namespace

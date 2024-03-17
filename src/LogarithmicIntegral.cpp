///
/// @file  LogarithmicIntegral.cpp
/// @brief The logarithmic integral function is a very accurate
///        approximation of PrimePi(x). The inverse logarithmic
///        integral function is a very accurate approximation of the
///        nth prime.
///
///        Note that these implementations are only accurate up to
///        10^15 (if the long double type has 80 bits). We also
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
template <typename T>
T li(T x)
{
  if (x <= 1)
    return 0;

  T gamma = (T) 0.577215664901532860606512090082402431L;
  T sum = 0;
  T inner_sum = 0;
  T factorial = 1;
  T p = -1;
  T q = 0;
  T power2 = 1;
  T logx = std::log(x);
  int k = 0;

  for (int n = 1; n < 1000; n++)
  {
    p *= -logx;
    factorial *= n;
    q = factorial * power2;
    power2 *= 2;

    for (; k <= (n - 1) / 2; k++)
      inner_sum += T(1.0) / (2 * k + 1);

    auto old_sum = sum;
    sum += (p / q) * inner_sum;

    // Not converging anymore
    if (std::abs(sum - old_sum) <= std::numeric_limits<T>::epsilon())
      break;
  }

  return gamma + std::log(logx) + std::sqrt(x) * sum;
}

/// Calculate the Eulerian logarithmic integral which is a very
/// accurate approximation of the number of primes <= x.
/// Li(x) > pi(x) for 24 <= x <= ~ 10^316
///
template <typename T>
T Li(T x)
{
  T li2 = (T) 1.045163780117492784844588889194613136L;
  T lix = li(x);
  return std::max(T(0), lix - li2);
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
template <typename T>
T Li_inverse(T x)
{
  if (x < 2)
    return 0;

  T t = x * std::log(x);
  T old_term = std::numeric_limits<T>::max();

  for (int i = 0; i < 100; i++)
  {
    T term = (Li(t) - x) * std::log(t);

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

  for (int n = 1; n < 1000; n++)
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
    if (fabsq(sum - old_sum) <= FLT128_EPSILON)
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
  __float128 lix = li(x);
  return (lix < li2) ? 0 : lix - li2;
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

  for (int i = 0; i < 100; i++)
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

template <typename FLOAT, typename T>
T Li_inverse_overflow_check(T x)
{
  FLOAT res = Li_inverse((FLOAT) x);

  // Prevent integer overflow
  if (res > (FLOAT) std::numeric_limits<T>::max())
    return std::numeric_limits<T>::max();
  else
    return (T) res;
}

} // namespace

namespace primecount {

int64_t Li(int64_t x)
{
#if defined(HAVE_FLOAT128)
  // The accuracy of our implementation depends on the precision (number
  // of bits) of the long double type and the accuracy of the math
  // functions from the libc. Since there are many different libc with
  // varying accuracy, it is impossible to know the exact threshold for
  // when we should switch to __float128. But 1e14 seems to work well in
  // practice.
  if (x > 1e14)
    return (int64_t) ::Li((__float128) x);
#endif

  if (x > 1e8)
    return (int64_t) ::Li((long double) x);
  else
    return (int64_t) ::Li((double) x);
}

int64_t Li_inverse(int64_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return Li_inverse_overflow_check<__float128>(x);
#endif

  if (x > 1e8)
    return Li_inverse_overflow_check<long double>(x);
  else
    return Li_inverse_overflow_check<double>(x);
}

#ifdef HAVE_INT128_T

int128_t Li(int128_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return (int128_t) ::Li((__float128) x);
#endif

  if (x > 1e8)
    return (int128_t) ::Li((long double) x);
  else
    return (int128_t) ::Li((double) x);
}

int128_t Li_inverse(int128_t x)
{
#if defined(HAVE_FLOAT128)
  if (x > 1e14)
    return Li_inverse_overflow_check<__float128>(x);
#endif

  if (x > 1e8)
    return Li_inverse_overflow_check<long double>(x);
  else
    return Li_inverse_overflow_check<double>(x);
}

#endif

} // namespace

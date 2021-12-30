///
/// @file  pi_deleglise_rivat.cpp
/// @brief 64-bit and 128-bit parallel implementations of the
///        Deleglise-Rivat prime counting algorithm.
///
///        Deleglise-Rivat formula:
///        pi(x) = pi(y) + S1(x, a) + S2(x, a) - 1 - P2(x, a)
///        S2(x, a) = S2_trivial(x, a) + S2_easy(x, a) + S2_hard(x, a)
///        with y = alpha * x^(1/3), a = pi(y)
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <int128_t.hpp>
#include <print.hpp>
#include <S.hpp>
#include <to_string.hpp>

#include <stdint.h>
#include <string>

using namespace primecount;

namespace {

/// Calculate the contribution of the special leaves
template <typename T>
T S2(T x,
     int64_t y,
     int64_t z,
     int64_t c,
     T s2_approx,
     int threads,
     bool is_print)
{
  T s2_trivial = S2_trivial(x, y, z, c, threads, is_print);
  T s2_easy = S2_easy(x, y, z, c, threads, is_print);
  T s2_hard_approx = s2_approx - (s2_trivial + s2_easy);
  T s2_hard = S2_hard(x, y, z, c, s2_hard_approx, threads, is_print);
  T s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_deleglise_rivat_64(int64_t x,
                              int threads,
                              bool is_print)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t z = x / y;
  int64_t pi_y = pi_noprint(y, threads);
  int64_t c = PhiTiny::get_c(y);

  if (is_print)
  {
    print("");
    print("=== pi_deleglise_rivat_64(x) ===");
    print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
    print(x, y, z, c, threads);
  }

  int64_t p2 = P2(x, y, threads, is_print);
  int64_t s1 = S1(x, y, c, threads, is_print);
  int64_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int64_t s2 = S2(x, y, z, c, s2_approx, threads, is_print);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

#if defined(HAVE_INT128_T)

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int128_t pi_deleglise_rivat_128(int128_t x,
                                int threads,
                                bool is_print)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  maxint_t limit = get_max_x(alpha);

  if (x > limit)
    throw primecount_error("pi(x): x must be <= " + to_string(limit));

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t pi_y = pi_noprint(y, threads);
  int64_t c = PhiTiny::get_c(y);

  if (is_print)
  {
    print("");
    print("=== pi_deleglise_rivat_128(x) ===");
    print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
    print(x, y, z, c, threads);
  }

  int128_t p2 = P2(x, y, threads, is_print);
  int128_t s1 = S1(x, y, c, threads, is_print);
  int128_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int128_t s2 = S2(x, y, z, c, s2_approx, threads, is_print);
  int128_t phi = s1 + s2;
  int128_t sum = phi + pi_y - 1 - p2;

  return sum;
}

#endif

} // namespace

///
/// @file  pi_deleglise_rivat_parallel1.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2.hpp>

#include <stdint.h>

using namespace std;
using namespace primecount;

namespace {

/// Calculate the contribution of the special leaves
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c,
           int64_t s2_approx,
           int threads)
{
  int64_t s2_trivial = S2_trivial(x, y, z, c, threads);
  int64_t s2_easy = S2_easy(x, y, z, c, threads);
  int64_t s2_hard_approx = s2_approx - (s2_trivial + s2_easy);
  int64_t s2_hard = S2_hard(x, y, z, c, s2_hard_approx, threads);
  int64_t s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_deleglise_rivat_parallel1(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t z = x / y;
  int64_t pi_y = pi_legendre(y, threads);
  int64_t c = PhiTiny::get_c(y);

  print("");
  print("=== pi_deleglise_rivat_parallel1(x) ===");
  print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
  print(x, y, z, c, alpha, threads);

  int64_t p2 = P2(x, y, threads);
  int64_t s1 = S1(x, y, c, threads);
  int64_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int64_t s2 = S2(x, y, z, c, s2_approx, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace

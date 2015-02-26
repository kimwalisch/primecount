///
/// @file  pi_deleglise_rivat_parallel3.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. This implementation is identical to
///        pi_deleglise_rivat_parallel2(x) but uses 128-bit integers.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <int128.hpp>
#include <S1.hpp>
#include "S2.hpp"

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Calculate the contribution of the special leaves.
/// @pre y > 0 && c > 1
///
int128_t S2(int128_t x,
            int64_t y,
            int64_t z,
            int64_t c,
            int128_t s2_approx,
            int threads)
{
  threads = validate_threads(threads, z);

  int128_t s2_trivial = S2_trivial(x, y, z, c, threads);
  int128_t s2_easy = S2_easy(x, y, z, c, threads);
  int128_t s2_hard_approx = s2_approx - (s2_trivial + s2_easy);
  int128_t s2_hard = S2_hard(x, y, z, c, s2_hard_approx, threads);
  int128_t s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

/// alpha is a tuning factor which should grow like (log(x))^3
/// for the Deleglise-Rivat prime counting algorithm.
///
double compute_alpha(int128_t x)
{
  double d = (double) x;
  double alpha = (get_alpha() >= 1) ? get_alpha() : log(d) * log(d) * log(d) / 1200;
  return in_between(1, alpha, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat_parallel3(int128_t x, int threads)
{
  if (x < 2)
    return 0;

  if (x > to_maxint(primecount::max()))
    throw primecount_error("pi(x): x must be <= " + max());

  double alpha = compute_alpha(x);
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = (int64_t) (x / y);
  int64_t pi_y = pi_legendre(y, 1);
  int64_t c;

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_deleglise_rivat_parallel3(x) ===" << endl;
    cout << "pi(x) = S1 + S2 + pi(y) - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    cout << "c = " << PhiTiny::max_a() << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }

  int128_t p2 = P2(x, y, threads);
  int128_t c = min(pi_y, PhiTiny::max_a());
  int128_t s1 = S1(x, y, c, threads);
  int128_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int128_t s2 = S2(x, y, z, c, s2_approx, threads);
  int128_t phi = s1 + s2;
  int128_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace

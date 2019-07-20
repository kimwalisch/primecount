///
/// @file  pi_gourdon1.cpp
/// @brief Simple demonstration implementation of Xavier Gourdon's
///        prime counting algorithm. Xavier Gourdon's algorithm is
///        an improved version of the Deleglise-Rivat algorithm.
///
///        Xavier Gourdon formula:
///        pi(x) = A - B + C + D + phi0 + Sigma
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <generate.hpp>
#include <PhiTiny.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <vector>

using namespace std;

namespace primecount {

/// Calculate the number of primes below x using
/// Xavier Gourdon's algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_gourdon1(int64_t x, int64_t y, int64_t z, int64_t k)
{
  if (x < 2)
    return 0;

  int64_t a = A(x, y, 1);
  int64_t b = B(x, y, 1);
  int64_t c = C(x, y, z, k, 1);
  int64_t d = 0;
  int64_t phi0 = Phi0(x, y, z, k, 1);
  int64_t sigma = Sigma(x, y, 1);

  int64_t sum = a - b + c + d + phi0 + sigma;

  return sum;
}

} // namespace

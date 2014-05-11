///
/// @file  pi_meissel.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "pmath.hpp"

#include <primecount.hpp>
#include <algorithm>
#include <stdint.h>

#ifdef _OPENMP
  #include <omp.h>
  #include "get_omp_threads.hpp"
#endif

using std::max;

namespace primecount {

/// Calculate the number of primes below x using Meissel's formula.
/// Run time: O(x/(log x)^3) operations, O(x^0.5 / log x) space.
///
int64_t pi_meissel(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  int64_t x13 = iroot<3>(x);
  int64_t a = pi_legendre(x13, 1);
  int64_t phi_xa, p2;

#ifdef _OPENMP
  threads = get_omp_threads(threads);
  omp_set_nested(true);
  #pragma omp parallel sections num_threads(threads)
  {
    #pragma omp section
    phi_xa = phi(x, a, max(1, threads - 1));
    #pragma omp section
    p2 = P2(x, x13);
  }
  omp_set_nested(false);
#else
  phi_xa = phi(x, a, threads);
  p2 = P2(x, x13);
#endif

  int64_t sum = phi_xa + a - 1 - p2;
  return sum;
}

} // namespace primecount

///
/// @file  phi.cpp
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiCache.h"
#include "PhiTiny.h"
#include "pi_bsearch.h"
#include "imath.h"

#include <primesieve.hpp>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
  #include "to_omp_threads.h"
#endif

namespace primecount {

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads)
{
  if (x < 1) return 0;
  if (a > x) return 1;
  if (a < 1) return x;

  static const PhiTiny phiTiny;

  if (phiTiny.is_cached(a))
    return phiTiny.phi(x, a);

  std::vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_n_primes(a, &primes);

  if (primes.at(a) >= x)
    return 1;

  int iters = pi_bsearch(primes, a, isqrt(x));
  PhiCache cache(primes, phiTiny);
  int64_t sum = x - a + iters;

#ifdef _OPENMP
  threads = to_omp_threads(threads);
  #pragma omp parallel for firstprivate(cache) reduction(+: sum) \
      num_threads(threads) schedule(dynamic, 16)
#endif
  for (int64_t a2 = 0; a2 < iters; a2++)
    sum += cache.phi_recursive<-1>(x / primes[a2 + 1], a2);

  return sum;
}

} // namespace primecount

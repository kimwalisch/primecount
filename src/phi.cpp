///
/// @file  phi.cpp
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <pmath.hpp>
#include <PhiCache.hpp>
#include <PhiTiny.hpp>
#include <utils.hpp>

#include <algorithm>
#include <stdint.h>
#include <vector>
#include <cassert>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

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

  if (is_phi_tiny(a))
    return phi_tiny(x, a);

  std::vector<int32_t> primes = generate_n_primes(a);

  if (primes.at(a) >= x)
    return 1;

  int64_t iters = pi_bsearch(primes, a, isqrt(x));
  PhiCache cache(primes);
  int64_t sum = x - a + iters;

  #pragma omp parallel for firstprivate(cache) schedule(dynamic, 16) \
      num_threads(validate_threads(threads)) reduction(+: sum)
  for (int64_t a2 = 0; a2 < iters; a2++)
    sum += cache.phi(x / primes[a2 + 1], a2, -1);

  return sum;
}

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, PhiCache* phiCache)
{
  assert(phiCache != 0);
  return phiCache->phi(x, a);
}

} // namespace primecount

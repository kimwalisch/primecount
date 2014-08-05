///
/// @file  P3.cpp
/// @brief 3rd partial sieve function.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <generate.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace primecount {

/// 3rd partial sieve function.
/// P3(x, a) counts the numbers <= x that have exactly 3 prime
/// factors each exceeding the a-th prime.
/// Space complexity: O(pi(sqrt(x))).
///
int64_t P3(int64_t x, int64_t a, int threads)
{
  std::vector<int32_t> primes = generate_primes(isqrt(x));
  int64_t y = iroot<3>(x);
  int64_t pi_y = pi_bsearch(primes, y);
  int64_t sum = 0;

  #pragma omp parallel for reduction(+: sum) schedule(dynamic) \
          num_threads(validate_threads(threads))
  for (int64_t i = a + 1; i <= pi_y; i++)
  {
    int64_t xi = x / primes[i];
    int64_t bi = pi_bsearch(primes, isqrt(xi));

    for (int64_t j = i; j <= bi; j++)
      sum += pi_bsearch(primes, xi / primes[j]) - (j - 1);
  }

  return sum;
}

} // namespace primecount

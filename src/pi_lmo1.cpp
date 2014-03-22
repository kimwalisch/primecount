///
/// @file  pi_lmo1.cpp
/// @brief Simple implementation of the Lagarias-Miller-Odlyzko prime
///        counting algorithm. This implementation does not use
///        segmentation to calculate S2(x) and the calculation of phi
///        in S2(x) does not use a special data structure which counts
///        the number of unsieved elements in O(log N) operations.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "PhiTiny.hpp"
#include "PhiCache.hpp"
#include "pmath.hpp"
#include "Pk.hpp"

#include <primecount.hpp>
#include <primesieve.hpp>
#include <stdint.h>
#include <vector>

using namespace primecount;

namespace {

/// Calculate the contribution of the ordinary leaves.
///
int64_t S1(int64_t x,
           int64_t x13_alpha,
           int64_t c,
           std::vector<int32_t>& primes,
           std::vector<int32_t>& lpf,
           std::vector<int32_t>& mu)
{
  int64_t S1_result = 0;

  const PhiTiny phiTiny;
  PhiCache cache(primes, phiTiny);

  for (int64_t n = 1; n <= x13_alpha; n++)
    if (lpf[n] > primes[c])
      S1_result += mu[n] * phi(x / n, c, &cache);

  return S1_result;
}

/// Calculate the contribution of the special leaves.
///
int64_t S2(int64_t x,
           int64_t x13_alpha,
           int64_t x23_alpha,
           int64_t a,
           int64_t c,
           std::vector<int32_t>& primes,
           std::vector<int32_t>& lpf,
           std::vector<int32_t>& mu)
{
  int64_t S2_result = 0;
  std::vector<char> sieve(x23_alpha + 1, 1);

  // phi(y, b) nodes with b <= c do not contribute to S2, so we
  // simply sieve out the multiples of the first c primes
  //
  for (int64_t b = 1; b <= c; b++)
  {
    int64_t prime = primes[b];
    for (int64_t i = prime; i <= x23_alpha; i += prime)
      sieve[i] = 0;
  }

  for (int64_t b = c; b < a; b++)
  {
    int64_t prime = primes[b + 1];
    int64_t phi = 0;
    int64_t old = 0;

    for (int64_t m = x13_alpha; m > x13_alpha / prime; m--)
    {
      if (mu[m] != 0 && prime < lpf[m])
      {
        // We have found a special leaf, compute it's contribution
        // phi(x / (m * primes[b + 1]), b) by counting the
        // number of unsieved elements below x / (m * primes[b + 1])
        // after having removed the multiples of the first b primes
        //
        int64_t q = x / (m * prime);
        for (int64_t k = old + 1; k <= q; k++)
          phi += sieve[k];
        old = q;
        S2_result -= mu[m] * phi;
      }
    }

    assert(prime > 2);
    // Remove the multiples of (b + 1)th prime
    for (int64_t i = prime; i <= x23_alpha; i += prime * 2)
      sieve[i] = 0;
  }

  return S2_result;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(2/3)) space.
/// @note O(x^(2/3)) space is because S2(x, c) is not segmented.
///
int64_t pi_lmo1(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  // Optimization factor, see:
  // Tomás Oliveira e Silva, Computing pi(x): the combinatorial method,
  // Revista do DETUA, vol. 4, no. 6, pp. 759-768, March 2006.
  double beta = 1.1;
  double alpha = std::max(1.0, nth_log<2>(x) * beta);

  int64_t x13 = iroot<3>(x);
  int64_t x13_alpha = (int64_t)(x13 * alpha);
  int64_t x23_alpha = (int64_t) std::pow((double) x / alpha, 2.0 / 3.0);
  int64_t a = pi_lehmer(x13_alpha);
  int64_t c = (a < 6) ? a : 6;

  std::vector<int32_t> primes;
  std::vector<int32_t> lpf = make_least_prime_factor(x13_alpha);
  std::vector<int32_t> mu = make_moebius(x13_alpha);

  primes.push_back(0);
  primesieve::generate_n_primes(a, &primes);

  int64_t phi = S1(x, x13_alpha, c, primes, lpf , mu) + S2(x, x13_alpha, x23_alpha, a, c, primes, lpf , mu);
  int64_t sum = phi + a - 1 - P2(x, a, threads);

  return sum;
}

} // namespace primecount

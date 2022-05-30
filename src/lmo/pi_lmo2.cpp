///
/// @file  pi_lmo2.cpp
/// @brief Simple demonstration implementation of the
///        Lagarias-Miller-Odlyzko prime counting algorithm.
///        This implementation uses the sieve of Eratosthenes
///        (without segmentation) to calculate S2(x).
///
///        Lagarias-Miller-Odlyzko formula:
///        pi(x) = pi(y) + S1(x, a) + S2(x, a) - 1 - P2(x, a)
///        with y = x^(1/3), a = pi(y)
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <imath.hpp>
#include <generate.hpp>
#include <pod_vector.hpp>
#include <PhiTiny.hpp>
#include <S.hpp>

#include <stdint.h>

using namespace primecount;

namespace {

/// Calculate the contribution of the special leaves.
/// This implementation uses the sieve of Eratosthenes
/// (without segmentation).
/// Memory usage: O(x^(2/3) * (log x)^2)
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t c,
           int64_t pi_y,
           const pod_vector<int32_t>& primes,
           const pod_vector<int32_t>& lpf,
           const pod_vector<int32_t>& mu)
{
  int64_t limit = x / y;
  int64_t s2 = 0;
  int64_t b = 1;
  pod_vector<bool> sieve(limit);
  std::fill(sieve.begin(), sieve.end(), 1);

  // phi(y, b) nodes with b <= c do not contribute to S2, so
  // we sieve out the multiples of the first c primes
  for (; b <= c; b++)
  {
    int64_t prime = primes[b];
    for (int64_t k = prime; k < limit; k += prime)
      sieve[k] = 0;
  }

  for (; b < pi_y; b++)
  {
    int64_t prime = primes[b];
    int64_t i = 1;
    int64_t phi = 0;

    for (int64_t m = y; m > y / prime; m--)
    {
      if (mu[m] != 0 && prime < lpf[m])
      {
        // We have found a special leaf. Compute it's contribution
        // phi(x / (primes[b] * m), b - 1) by counting the number
        // of unsieved elements <= x / (primes[b] * m) after having
        // removed the multiples of the first b - 1 primes.
        //
        for (int64_t xpm = x / (prime * m); i <= xpm; i++)
          phi += sieve[i];

        s2 -= mu[m] * phi;
      }
    }

    // remove the multiples of b-th prime
    for (int64_t k = prime; k < limit; k += prime * 2)
      sieve[k] = 0;
  }

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3))
/// Memory usage: O(x^(2/3) / (log x)^2)
///
int64_t pi_lmo2(int64_t x)
{
  if (x < 2)
    return 0;

  bool threads = 1;
  double alpha = get_alpha_lmo(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t c = PhiTiny::get_c(y);
  int64_t p2 = P2(x, y, threads);

  auto primes = generate_primes<int32_t>(y);
  auto lpf = generate_lpf(y);
  auto mu = generate_moebius(y);

  int64_t pi_y = primes.size() - 1;
  int64_t s1 = S1(x, y, c, threads);
  int64_t s2 = S2(x, y, c, pi_y, primes, lpf, mu);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace

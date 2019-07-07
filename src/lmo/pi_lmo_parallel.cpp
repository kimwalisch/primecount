///
/// @file  pi_lmo_parallel.cpp
/// @brief Parallel implementation of the Lagarias-Miller-Odlyzko
///        prime counting algorithm using OpenMP. This implementation
///        uses load balancing and counts the number of unsieved
///        elements using POPCNT without using any special counting
///        tree data structure.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <Sieve.hpp>
#include <generate.hpp>
#include <generate_phi.hpp>
#include <LoadBalancer.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <PhiTiny.hpp>
#include <PiTable.hpp>
#include <print.hpp>
#include <S1.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Compute the S2 contribution of the interval
/// [low, low + segments * segment_size[
///
int64_t S2_thread(int64_t x,
                  int64_t y,
                  int64_t z,
                  int64_t c,
                  int64_t low,
                  int64_t segments,
                  int64_t segment_size,
                  PiTable& pi,
                  vector<int32_t>& primes,
                  vector<int32_t>& lpf,
                  vector<int32_t>& mu,
                  Runtime& runtime)
{
  int64_t low1 = max(low, 1);
  int64_t limit = min(low + segments * segment_size, z + 1);
  int64_t size = pi[min(isqrt(x / low1), y)] + 1;
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_y = pi[y];
  int64_t S2_thread = 0;

  if (c >= size - 1)
    return 0;

  runtime.init_start();
  Sieve sieve(low, segment_size, size - 1);
  auto phi = generate_phi(low, size - 1, primes, pi);
  runtime.init_stop();

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    low1 = max(low, 1);

    // pre-sieve multiples of first c primes
    int64_t count_low_high = sieve.pre_sieve(c, low, high);
    int64_t b = c + 1;

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy:  mu[m] != 0 && primes[b] < lpf[m], low <= (x / n) < high
    for (; b <= min(pi_sqrty, size - 1); b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low1), y);
      int64_t count = 0;
      int64_t i = 0;

      if (prime >= max_m)
        goto next_segment;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          int64_t xn = x / (prime * m);
          int64_t stop = xn - low;
          count += sieve.count(i, stop);
          i = stop + 1;
          int64_t phi_xn = phi[b] + count;
          S2_thread -= mu[m] * phi_xn;
        }
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    // For pi_sqrty < b < pi_y
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high
    for (; b < min(pi_y, size); b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min(x / (prime * low1), y)];
      int64_t min_m = max(x / (prime * high), prime);
      int64_t count = 0;
      int64_t i = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xn = x / (prime * primes[l]);
        int64_t stop = xn - low;
        count += sieve.count(i, stop);
        i = stop + 1;
        int64_t phi_xn = phi[b] + count;
        S2_thread += phi_xn;
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    next_segment:;
  }

  return S2_thread;
}

/// Calculate the contribution of the special leaves.
/// This is a parallel implementation with advanced load balancing.
/// As most special leaves tend to be in the first segments we
/// start off with a small segment size and few segments
/// per thread, after each iteration we dynamically increase
/// the segment size and the number of segments.
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c,
           int64_t s2_approx,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu,
           int threads)
{
  print("");
  print("=== S2(x, y) ===");
  print("Computation of the special leaves");

  double time = get_time();
  threads = ideal_num_threads(threads, z);
  LoadBalancer loadBalancer(x, y, z, s2_approx);
  PiTable pi(y);

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    int64_t low = 0;
    int64_t segments = 0;
    int64_t segment_size = 0;
    int64_t S2 = 0;
    Runtime runtime;

    while (loadBalancer.get_work(&low, &segments, &segment_size, S2, runtime))
    {
      runtime.start();
      S2 = S2_thread(x, y, z, c, low, segments, segment_size, pi, primes, lpf, mu, runtime);
      runtime.stop();
    }
  }

  int64_t S2 = (int64_t) loadBalancer.get_result();
  print("S2", S2, time);

  return S2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x)
/// Memory usage: O(x^(1/3) * (log x)^2)
///
int64_t pi_lmo_parallel(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha_lmo(x);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t z = x / y;
  int64_t c = PhiTiny::get_c(y);

  print("");
  print("=== pi_lmo_parallel(x) ===");
  print("pi(x) = S1 + S2 + pi(y) - 1 - P2");
  print(x, y, z, c, alpha, threads);

  int64_t p2 = P2(x, y, threads);
  auto primes = generate_primes<int32_t>(y);
  auto lpf = generate_lpf(y);
  auto mu = generate_moebius(y);

  int64_t pi_y = primes.size() - 1;
  int64_t s1 = S1(x, y, c, threads);
  int64_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int64_t s2 = S2(x, y, z, c, s2_approx, primes, lpf, mu, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace

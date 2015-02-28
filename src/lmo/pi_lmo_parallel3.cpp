///
/// @file  pi_lmo_parallel3.cpp
/// @brief Parallel implementation of the Lagarias-Miller-Odlyzko
///        prime counting algorithm using OpenMP. This implementation
///        uses improved load balancing and counts the number of
///        unsieved elements using POPCNT without using any special
///        counting tree data structure.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2LoadBalancer.hpp>
#include <S2Status.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// For each prime calculate its first multiple >= low
vector<int64_t> generate_next_multiples(int64_t low, int64_t size, vector<int32_t>& primes)
{
  vector<int64_t> next;

  next.reserve(size);
  next.push_back(0);

  for (int64_t b = 1; b < size; b++)
  {
    int64_t prime = primes[b];
    int64_t next_multiple = ceil_div(low, prime) * prime;
    next_multiple += prime * (~next_multiple & 1);
    next.push_back(next_multiple);
  }

  return next;
}

/// Cross-off the multiples of prime in the sieve array.
/// @return  Count of crossed-off elements.
///
int64_t cross_off(int64_t prime,
                  int64_t low,
                  int64_t high,
                  int64_t& next_multiple,
                  BitSieve& sieve)
{
  int64_t unset = 0;
  int64_t k = next_multiple;

  for (; k < high; k += prime * 2)
  {
    if (sieve[k - low])
    {
      sieve.unset(k - low);
      unset++;
    }
  }

  next_multiple = k;
  return unset;
}

/// Compute the S2 contribution for the interval
/// [low_process, low_process + segments * segment_size[.
/// The missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the calling (parent) S2 function.
///
int64_t S2_thread(int64_t x,
                  int64_t y,
                  int64_t c,
                  int64_t segment_size,
                  int64_t segments_per_thread,
                  int64_t thread_num,
                  int64_t low,
                  int64_t limit,
                  vector<int32_t>& pi,
                  vector<int32_t>& primes,
                  vector<int32_t>& lpf,
                  vector<int32_t>& mu,
                  vector<int64_t>& mu_sum,
                  vector<int64_t>& phi)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t size = pi[min(isqrt(x / low), y)] + 1;
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_y = pi[y];

  if (c >= size - 1)
    return 0;

  int64_t S2_thread = 0;
  BitSieve sieve(segment_size);
  vector<int64_t> next = generate_next_multiples(low, size, primes);
  phi.resize(size, 0);
  mu_sum.resize(size, 0);

  // Process the segments assigned to the current thread
  for (; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = 2;

    sieve.fill(low, high);

    // phi(y, b) nodes with b <= c do not contribute to S2, so we
    // simply sieve out the multiples of the first c primes
    for (; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime * 2)
        sieve.unset(k - low);
      next[b] = k;
    }

    int64_t count_low_high = sieve.count((high - 1) - low);

    // For c + 1 <= b < pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy:  mu[m] != 0 && primes[b] < lpf[m], low <= (x / n) < high
    for (; b < min(pi_sqrty, size); b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);
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
          mu_sum[b] -= mu[m];
        }
      }

      phi[b] += count_low_high;
      count_low_high -= cross_off(prime, low, high, next[b], sieve);
    }

    // For pi_sqrty <= b < pi_y
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high
    for (; b < min(pi_y, size); b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min(x / (prime * low), y)];
      int64_t min_m = max3(x / (prime * high), y / prime, prime);
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
        mu_sum[b]++;
      }

      phi[b] += count_low_high;
      count_low_high -= cross_off(prime, low, high, next[b], sieve);
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
/// the segment size and the segments per thread.
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t c,
           int64_t s2_approx,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu,
           int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2(x, y) ===" << endl;
    cout << "Computation of the special leaves" << endl;
  }

  int64_t S2_total = 0;
  int64_t low = 1;
  int64_t limit = x / y + 1;
  threads = validate_threads(threads, limit);

  S2Status status;
  S2LoadBalancer loadBalancer(x, limit, threads);
  int64_t segment_size = loadBalancer.get_min_segment_size();
  int64_t segments_per_thread = 1;

  double time = get_wtime();
  vector<int32_t> pi = generate_pi(y);
  vector<int64_t> phi_total(primes.size(), 0);

  while (low < limit)
  {
    int64_t segments = ceil_div(limit - low, segment_size);
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, ceil_div(segments, threads));

    aligned_vector<vector<int64_t> > phi(threads);
    aligned_vector<vector<int64_t> > mu_sum(threads);
    aligned_vector<double> timings(threads);

    #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
    for (int i = 0; i < threads; i++)
    {
      timings[i] = get_wtime();
      S2_total += S2_thread(x, y, c, segment_size, segments_per_thread,
          i, low, limit, pi, primes, lpf, mu, mu_sum[i], phi[i]);
      timings[i] = get_wtime() - timings[i];
    }

    // Once all threads have finished reconstruct and add the 
    // missing contribution of all special leaves. This must
    // be done in order as each thread (i) requires the sum of
    // the phi values from the previous threads.
    //
    for (int i = 0; i < threads; i++)
    {
      for (size_t j = 1; j < phi[i].size(); j++)
      {
        S2_total += phi_total[j] * mu_sum[i][j];
        phi_total[j] += phi[i][j];
      }
    }

    low += segments_per_thread * threads * segment_size;
    loadBalancer.update(low, threads, &segment_size, &segments_per_thread, timings);

    if (print_status())
      status.print(S2_total, s2_approx, loadBalancer.get_rsd());
  }

  if (print_status())
    print_result("S2", S2_total, time);

  return S2_total;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3) / log x) operations, O(x^(1/3) * (log x)^2) space.
///
int64_t pi_lmo_parallel3(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha(x, 1e15, 2, 400);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_lmo_parallel3(x) ===" << endl;
    cout << "pi(x) = S1 + S2 + pi(y) - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "alpha = " << fixed << setprecision(3) << alpha << endl;
    cout << "c = " << PhiTiny::max_a() << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }

  int64_t p2 = P2(x, y, threads);

  vector<int32_t> mu = generate_moebius(y);
  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> primes = generate_primes(y);

  int64_t pi_y = primes.size() - 1;
  int64_t c = min(pi_y, PhiTiny::max_a());
  int64_t s1 = S1(x, y, c, lpf, mu, threads);
  int64_t s2_approx = S2_approx(x, pi_y, p2, s1);
  int64_t s2 = S2(x, y, c, s2_approx, primes, lpf, mu, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount

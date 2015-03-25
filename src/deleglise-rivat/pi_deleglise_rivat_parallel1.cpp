///
/// @file  pi_deleglise_rivat_parallel1.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. In the Deleglise-Rivat algorithm there
///        are 3 additional types of special leaves compared to the 
///        Lagarias-Miller-Odlyzko algorithm: trivial special leaves,
///        clustered easy leaves and sparse easy leaves.
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "S2.hpp"
#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <S2LoadBalancer.hpp>
#include <tos_counters.hpp>
#include <S1.hpp>

#include <stdint.h>
#include <algorithm>
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

template <typename T>
void cross_off(int64_t prime,
               int64_t low,
               int64_t high,
               int64_t& next_multiple,
               BitSieve& sieve,
               T& counters)
{
  int64_t segment_size = sieve.size();
  int64_t k = next_multiple;

  for (; k < high; k += prime * 2)
  {
    if (sieve[k - low])
    {
      sieve.unset(k - low);
      cnt_update(counters, k - low, segment_size);
    }
  }
  next_multiple = k;
}

/// Compute the S2 contribution of the hard special leaves.
/// Each thread processes the interval
/// [low_thread, low_thread + segments * segment_size[
/// and the missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the parent S2_hard() function.
///
int64_t S2_hard_thread(int64_t x,
                        int64_t y,
                        int64_t z,
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
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t max_prime = min3(isqrt(x / low), isqrt(z), y);
  int64_t pi_max = pi[max_prime];

  if (c > pi_max)
    return 0;

  int64_t S2_thread = 0;
  BitSieve sieve(segment_size);
  vector<int32_t> counters(segment_size);
  vector<int64_t> next = generate_next_multiples(low, pi_max + 1, primes);
  phi.resize(pi_max + 1, 0);
  mu_sum.resize(pi_max + 1, 0);

  // segmeted sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    sieve.fill(low, high);

    // phi(y, i) nodes with i <= c do not contribute to S2, so we
    // simply sieve out the multiples of the first c primes
    for (int64_t i = 2; i <= c; i++)
    {
      int64_t k = next[i];
      for (int64_t prime = primes[i]; k < high; k += prime * 2)
        sieve.unset(k - low);
      next[i] = k;
    }

    // Initialize special tree data structure from sieve
    cnt_finit(sieve, counters, segment_size);

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m, with mu[m] != 0 and primes[b] < lpf[m]
    // which satisfy: low <= (x / n) < high
    for (int64_t end = min(pi_sqrty, pi_max); b <= end; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);

      if (prime >= max_m)
        goto next_segment;

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (mu[m] != 0 && prime < lpf[m])
        {
          int64_t n = prime * m;
          int64_t count = cnt_query(counters, (x / n) - low);
          int64_t phi_xn = phi[b] + count;
          S2_thread -= mu[m] * phi_xn;
          mu_sum[b] -= mu[m];
        }
      }

      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    // For pi_sqrty <= b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= pi_max; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min3(x / (prime * low), z / prime, y)];
      int64_t min_hard_leaf = max3(x / (prime * high), y / prime, prime);

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_hard_leaf; l--)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        int64_t count = cnt_query(counters, xn - low);
        int64_t phi_xn = phi[b] + count;
        S2_thread += phi_xn;
        mu_sum[b]++;
      }

      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    next_segment:;
  }

  return S2_thread;
}

/// Calculate the contribution of the special leaves which require
/// a sieve (in order to reduce the memory usage).
/// This is a parallel implementation with advanced load balancing.
/// As most special leaves tend to be in the first segments we
/// start off with a small segment size and few segments
/// per thread, after each iteration we dynamically increase
/// the segment size and the segments per thread.
///
int64_t S2_hard(int64_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 vector<int32_t>& pi,
                 vector<int32_t>& primes,
                 vector<int32_t>& lpf,
                 vector<int32_t>& mu,
                 int threads)
{
  int64_t S2_total = 0;
  int64_t low = 1;
  int64_t limit = z + 1;

  S2LoadBalancer loadBalancer(x, y, z, threads);
  int64_t segment_size = loadBalancer.get_min_segment_size();
  int64_t segments_per_thread = 1;
  vector<int64_t> phi_total(pi[min(isqrt(z), y)] + 1, 0);

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
      S2_total += S2_hard_thread(x, y, z, c, segment_size, segments_per_thread,
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
  }

  return S2_total;
}

/// Calculate the contribution of the special leaves.
/// @pre y > 0 && c > 1
///
int64_t S2(int64_t x,
           int64_t y,
           int64_t z,
           int64_t c,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu,
           int threads)
{
  int64_t limit = z + 1;
  threads = validate_threads(threads, limit);
  vector<int32_t> pi = generate_pi(y);

  int64_t s2_trivial = S2_trivial(x, y, z, c, threads);
  int64_t s2_easy = S2_easy(x, y, z, c, threads);
  int64_t s2_hard = S2_hard(x, y, z, c, pi, primes, lpf, mu, threads);
  int64_t s2 = s2_trivial + s2_easy + s2_hard;

  return s2;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat_parallel1(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = get_alpha(x, 0.0017154, -0.0508992, 0.483613, 0.0672202);
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t) (x13 * alpha);
  int64_t z = x / y;
  int64_t p2 = P2(x, y, threads);

  vector<int32_t> mu = generate_moebius(y);
  vector<int32_t> lpf = generate_least_prime_factors(y);
  vector<int32_t> primes = generate_primes(y);

  int64_t pi_y = pi_bsearch(primes, y);
  int64_t c = PhiTiny::get_c(y);
  int64_t s1 = S1(x, y, c, threads);
  int64_t s2 = S2(x, y, z, c, primes, lpf, mu, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount

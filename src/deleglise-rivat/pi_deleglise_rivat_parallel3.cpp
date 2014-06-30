///
/// @file  pi_deleglise_rivat_parallel3.cpp
/// @brief Parallel implementation of the Lagarias-Miller-Odlyzko
///        prime counting algorithm with the improvements of Deleglise
///        and Rivat. This version uses compression (see FactorTable,
///        PiTable) to reduce the memory usage.
/// 
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <balance_S2_load.hpp>
#include <BitSieve.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <tos_counters.hpp>
#include <utils.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// For each prime calculate its first multiple >= low.
template <typename T1, typename T2>
void init_next_multiples(T1& next, T2& primes, int64_t size, int64_t low)
{
  next.reserve(size);
  next.push_back(0);

  for (int64_t b = 1; b < size; b++)
  {
    int64_t prime = primes[b];
    int64_t next_multiple = ((low + prime - 1) / prime) * prime;
    next_multiple += prime * (~next_multiple & 1);
    next.push_back(next_multiple);
  }
}

template <typename T1, typename T2>
void cross_off(int64_t prime, int64_t low, int64_t high, int64_t& next_multiple, T1& sieve, T2& counters)
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

/// Compute the S2 contribution for the interval
/// [low_process, low_process + segments * segment_size[.
/// The missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the calling (parent) S2 function.
///
int64_t S2_thread(int64_t x,
                  int64_t y,
                  int64_t z,
                  int64_t c,
                  int64_t segment_size,
                  int64_t segments_per_thread,
                  int64_t thread_num,
                  int64_t low,
                  int64_t limit,
                  FactorTable& factors,
                  PiTable& pi,
                  vector<int32_t>& primes,
                  vector<int64_t>& mu_sum,
                  vector<int64_t>& phi)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t pi_y = pi(y);
  int64_t pi_sqrty = pi(isqrt(y));
  int64_t max_prime = min(isqrt(x / low), y);
  int64_t max_index = pi(max_prime);
  int64_t phi_size = pi(min(isqrt(z), max_prime)) + 1;
  int64_t S2_thread = 0;

  BitSieve sieve(segment_size);
  vector<int32_t> counters(segment_size);
  vector<int64_t> next;
  init_next_multiples(next, primes, phi_size, low);
  phi.resize(phi_size, 0);
  mu_sum.resize(phi_size, 0);

  // Process the segments assigned to the current thread
  for (; low < limit; low += segment_size)
  {
    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = c + 1;

    // check if we need the sieve
    if (c < phi_size)
    {
      sieve.memset(low);

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
    }

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m, with mu[m] != 0 and primes[b] < lpf[m]
    // which satisfy: low <= (x / n) < high
    for (int64_t end = min(pi_sqrty, max_index); b <= end; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t max_m = min(x / (prime * low), y);

      if (prime >= max_m)
        goto next_segment;

      FactorTable::to_index(&min_m);
      FactorTable::to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (prime < factors.lpf(m))
        {
          int64_t n = prime * factors.get_number(m);
          int64_t count = cnt_query(counters, (x / n) - low);
          int64_t phi_xn = phi[b] + count;
          int64_t mu_m = factors.mu(m);
          S2_thread -= mu_m * phi_xn;
          mu_sum[b] -= mu_m;
        }
      }

      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    // For pi_sqrty <= b < pi_y
    // Find all special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (int64_t end = min(pi_y, max_index + 1); b < end; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi(min(x / (prime * low), y));

      if (prime >= primes[l])
        goto next_segment;

      int64_t min_hard_leaf = max(x / (prime * high), y / prime);
      min_hard_leaf = in_between(prime, min_hard_leaf, y);
      int64_t min_trivial_leaf = min(x / (prime * prime), y);
      int64_t min_clustered_easy_leaf = min(isqrt(x / prime), y);
      int64_t min_sparse_easy_leaf = min(z / prime, y);

      min_trivial_leaf = max(min_hard_leaf, min_trivial_leaf);
      min_clustered_easy_leaf = max(min_hard_leaf, min_clustered_easy_leaf);
      min_sparse_easy_leaf = max(min_hard_leaf, min_sparse_easy_leaf);

      // Find all trivial leaves which satisfy:
      // phi(x / (primes[b] * primes[l]), b - 1) = 1
      if (primes[l] > min_trivial_leaf)
      {
        int64_t l_min = pi(min_trivial_leaf);
        S2_thread += l - l_min;
        l = l_min;
      }

      // Find all clustered easy leaves which satisfy:
      // x / n <= y such that phi(x / n, b - 1) = pi(x / n) - b + 2
      // And phi(x / n, b - 1) == phi(x / m, b - 1)
      while (primes[l] > min_clustered_easy_leaf)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        assert(xn < isquare(primes[b]));
        int64_t phi_xn = pi(xn) - b + 2;
        int64_t m = prime * primes[b + phi_xn - 1];
        int64_t xm = max(x / m, min_clustered_easy_leaf);
        int64_t l2 = pi(xm);
        S2_thread += phi_xn * (l - l2);
        l = l2;
      }

      // Find all sparse easy leaves which satisfy:
      // x / n <= y such that phi(x / n, b - 1) = pi(x / n) - b + 2
      for (; primes[l] > min_sparse_easy_leaf; l--)
      {
        int64_t n = prime * primes[l];
        int64_t xn = x / n;
        assert(xn < isquare(primes[b]));
        S2_thread += pi(xn) - b + 2;
      }

      if (b < phi_size)
      {
        // Find all hard leaves which satisfy:
        // low <= (x / n) < high
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
           int64_t z,
           int64_t c,
           vector<int32_t>& primes,
           FactorTable& factors,
           int threads)
{
  int64_t limit = z + 1;
  threads = validate_threads(threads, limit);

  int64_t S2_total = 0;
  int64_t low = 1;
  int64_t sqrt_limit = isqrt(limit);
  int64_t logx = max(1, ilog(x));
  int64_t min_segment_size = 1 << 6;
  int64_t segments_per_thread = 1;
  int64_t segment_size = next_power_of_2(sqrt_limit / (logx * threads));
  segment_size = max(segment_size, min_segment_size);

  PiTable pi(y);
  vector<int64_t> phi_total(pi(min(isqrt(z), y)) + 1, 0);
  double relative_standard_deviation = 30;

  while (low < limit)
  {
    int64_t segments = (limit - low + segment_size - 1) / segment_size;
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, (segments + threads - 1) / threads);

    aligned_vector<vector<int64_t> > phi(threads);
    aligned_vector<vector<int64_t> > mu_sum(threads);
    aligned_vector<double> timings(threads);

    #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
    for (int i = 0; i < threads; i++)
    {
      timings[i] = get_wtime();
      S2_total += S2_thread(x, y, z, c, segment_size, segments_per_thread,
          i, low, limit, factors, pi, primes, mu_sum[i], phi[i]);
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
    balance_S2_load(&segment_size, &segments_per_thread, min_segment_size,
        sqrt_limit, &relative_standard_deviation, timings);
  }

  return S2_total;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * log x) space.
///
int64_t pi_deleglise_rivat_parallel3(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  // alpha is a tuning factor
  double alpha = in_between(1, log((double) x), iroot<6>(x));
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = x / y;

  vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(y, &primes);
  FactorTable factors(y);

  int64_t pi_y = pi_bsearch(primes, y);
  int64_t c = min<int64_t>(PhiTiny::MAX_A, pi_y);
  int64_t s1 = S1(x, y, c, primes, factors);
  int64_t s2 = S2(x, y, z, c, primes, factors, threads);
  int64_t p2 = P2(x, y, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount

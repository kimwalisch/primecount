///
/// @file  pi_deleglise_rivat_parallel2.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. Compared to
///        pi_deleglise_rivat_parallel1.cpp this version uses
///        compression (FactorTable & PiTable) to reduce the memory
///        usage.
///        
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2LoadBalancer.hpp>
#include <tos_counters.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
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

/// Calculate the contribution of the trivial leaves.
///
int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   PiTable& pi,
                   vector<int32_t>& primes,
                   int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_trivial(x, y) ===" << endl;
    cout << "Computation of the trivial special leaves" << endl;
  }

  int64_t pi_y = pi(y);
  int64_t pi_sqrtz = pi(min(isqrt(z), y));
  int64_t S2_total = 0;
  double time = get_wtime();

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrtz + 1); b < pi_y; b++)
  {
    int64_t prime = primes[b];
    S2_total += pi_y - pi(max(x / (prime * prime), prime));
  }

  if (print_status())
    print_result("S2_trivial", S2_total, time);

  return S2_total;
}

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                PiTable& pi,
                vector<int32_t>& primes,
                int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_easy(x, y) ===" << endl;
    cout << "Computation of the easy special leaves" << endl;
  }

  int64_t pi_sqrty = pi(isqrt(y));
  int64_t pi_x13 = pi(iroot<3>(x));
  int64_t S2_total = 0;
  double time = get_wtime();

  #pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
  {
    int64_t prime = primes[b];
    int64_t min_trivial_leaf = x / (prime * prime);
    int64_t min_clustered_easy_leaf = isqrt(x / prime);
    int64_t min_sparse_easy_leaf = z / prime;
    int64_t min_hard_leaf = max(y / prime, prime);

    min_sparse_easy_leaf = max(min_sparse_easy_leaf, min_hard_leaf);
    min_clustered_easy_leaf = max(min_clustered_easy_leaf, min_hard_leaf);
    int64_t l = pi(min(min_trivial_leaf, y));
    int64_t S2_result = 0;

    // Find all clustered easy leaves:
    // x / n <= y and phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (primes[l] > min_clustered_easy_leaf)
    {
      int64_t n = prime * primes[l];
      int64_t xn = x / n;
      assert(xn < isquare(primes[b]));
      int64_t phi_xn = pi(xn) - b + 2;
      int64_t m = prime * primes[b + phi_xn - 1];
      int64_t xm = max(x / m, min_clustered_easy_leaf);
      int64_t l2 = pi(xm);
      S2_result += phi_xn * (l - l2);
      l = l2;
    }

    // Find all sparse easy leaves:
    // x / n <= y and phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; primes[l] > min_sparse_easy_leaf; l--)
    {
      int64_t n = prime * primes[l];
      int64_t xn = x / n;
      assert(xn < isquare(primes[b]));
      S2_result += pi(xn) - b + 2;
    }

    S2_total += S2_result;
  }

  if (print_status())
    print_result("S2_easy", S2_total, time);

  return S2_total;
}

/// Compute the S2 contribution of the special leaves that require
/// a sieve. Each thread processes the interval
/// [low_thread, low_thread + segments * segment_size[
/// and the missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the parent S2_sieve() function.
///
int64_t S2_sieve_thread(int64_t x,
                        int64_t y,
                        int64_t z,
                        int64_t c,
                        int64_t segment_size,
                        int64_t segments_per_thread,
                        int64_t thread_num,
                        int64_t low,
                        int64_t limit,
                        FactorTable<uint16_t>& factors,
                        PiTable& pi,
                        vector<int32_t>& primes,
                        vector<int64_t>& mu_sum,
                        vector<int64_t>& phi)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t pi_sqrty = pi(isqrt(y));
  int64_t max_prime = min3(isqrt(x / low), isqrt(z), y);
  int64_t pi_max = pi(max_prime);
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

    // check if we need the sieve
    if (c <= pi_max)
    {
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
    }

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

      factors.to_index(&min_m);
      factors.to_index(&max_m);

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

    // For pi_sqrty <= b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= pi_max; b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi(min3(x / (prime * low), z / prime, y));
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
int64_t S2_sieve(int64_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 PiTable& pi,
                 vector<int32_t>& primes,
                 FactorTable<uint16_t>& factors,
                 int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_sieve(x, y) ===" << endl;
    cout << "Computation of the special leaves requiring a sieve" << endl;
  }

  int64_t S2_total = 0;
  int64_t low = 1;
  int64_t limit = z + 1;
  double time = get_wtime();

  S2LoadBalancer loadBalancer(x, limit, threads);
  int64_t segment_size = loadBalancer.get_min_segment_size();
  int64_t segments_per_thread = 1;
  vector<int64_t> phi_total(pi(min(isqrt(z), y)) + 1, 0);

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
      S2_total += S2_sieve_thread(x, y, z, c, segment_size, segments_per_thread,
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
    loadBalancer.update(low, threads, &segment_size, &segments_per_thread, timings);
  }

  if (print_status())
    print_result("S2_sieve", S2_total, time);

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
           FactorTable<uint16_t>& factors,
           int threads)
{
  int64_t S2_total = 0;
  int64_t limit = z + 1;
  threads = validate_threads(threads, limit);
  PiTable pi(y);

  S2_total += S2_trivial(x, y, z, c, pi, primes, threads);
  S2_total += S2_easy(x, y, z, c, pi, primes, threads);
  S2_total += S2_sieve(x, y, z, c, pi, primes, factors, threads);

  return S2_total;
}

/// alpha is a tuning factor which should grow like (log(x))^3
/// for the Deleglise-Rivat prime counting algorithm.
///
double compute_alpha(int64_t x)
{
  double d = (double) x;
  double alpha = log(d) * log(d) * log(d) / 1200;
  return in_between(1, alpha, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int64_t pi_deleglise_rivat_parallel2(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double alpha = compute_alpha(x);
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = x / y;

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_deleglise_rivat_parallel2(x) ===" << endl;
    cout << "pi(x) = S1 + S2 + pi(y) - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << PhiTiny::max_a() << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }

  int64_t p2 = P2(x, y, threads);
  vector<int32_t> primes = generate_primes(y);
  FactorTable<uint16_t> factors(y);

  int64_t pi_y = pi_bsearch(primes, y);
  int64_t c = min(pi_y, PhiTiny::max_a());
  int64_t s1 = S1(x, y, c, primes[c], factors, threads);
  int64_t s2 = S2(x, y, z, c, primes, factors, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount

///
/// @file  pi_deleglise_rivat_parallel3.cpp
/// @brief Parallel implementation of the Deleglise-Rivat prime
///        counting algorithm. This implementation is identical to
///        pi_deleglise_rivat_parallel2(x) but uses 128-bit integers.
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
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <aligned_vector.hpp>
#include <BitSieve.hpp>
#include <generate.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <int128.hpp>
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
template <typename T>
vector<int64_t> generate_next_multiples(int64_t low, int64_t size, vector<T>& primes)
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
template <typename P>
int128_t S2_trivial(uint128_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   PiTable& pi,
                   vector<P>& primes,
                   int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_trivial(x, y) ===" << endl;
    cout << "Computation of the trivial special leaves" << endl;
  }

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[min(isqrt(z), y)];
  int128_t S2_total = 0;
  double time = get_wtime();

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrtz + 1); b < pi_y; b++)
  {
    uint128_t prime = primes[b];
    uint64_t xn = (uint64_t) max(x / (prime * prime), prime);
    S2_total += pi_y - pi(xn);
  }

  if (print_status())
    print_result("S2_trivial", S2_total, time);

  return S2_total;
}

/// Compute the S2 contribution of the special leaves that require
/// a sieve. Each thread processes the interval
/// [low_thread, low_thread + segments * segment_size[
/// and the missing special leaf contributions for the interval
/// [1, low_process[ are later reconstructed and added in
/// the parent S2_sieve() function.
///
template <typename P, typename F>
int128_t S2_sieve_thread(uint128_t x,
                         int64_t y,
                         int64_t z,
                         int64_t c,
                         int64_t segment_size,
                         int64_t segments_per_thread,
                         int64_t thread_num,
                         int64_t low,
                         int64_t limit,
                         FactorTable<F>& factors,
                         PiTable& pi,
                         vector<P>& primes,
                         vector<int64_t>& mu_sum,
                         vector<int64_t>& phi)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t max_prime = min(min(isqrt(x / low), y), isqrt(z));
  int64_t pi_max = pi[max_prime];
  int128_t S2_thread = 0;

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
      int128_t prime128 = primes[b];
      int64_t prime = primes[b];
      int64_t min_m = max(min(x / (prime128 * high), y), y / prime);
      int64_t max_m = min(x / (prime128 * low), y);

      if (prime >= max_m)
        goto next_segment;

      factors.to_index(&min_m);
      factors.to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        if (prime < factors.lpf(m))
        {
          int64_t n = prime * factors.get_number(m);
          int64_t xn = (int64_t) (x / n);
          int64_t count = cnt_query(counters, xn - low);
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
      int128_t prime128 = primes[b];
      int64_t prime = primes[b];
      int64_t l = pi(min(min(x / (prime128 * low), y), z / prime));
      int64_t min_hard_leaf = max3(min(x / (prime128 * high), y), y / prime, prime);

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
template <typename P, typename F>
int128_t S2_sieve(int128_t x,
                  int64_t y,
                  int64_t z,
                  int64_t c,
                  PiTable& pi,
                  vector<P>& primes,
                  FactorTable<F>& factors,
                  int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_sieve(x, y) ===" << endl;
    cout << "Computation of the special leaves requiring a sieve" << endl;
  }

  int128_t S2_total = 0;
  int64_t low = 1;
  int64_t limit = z + 1;
  double time = get_wtime();

  S2LoadBalancer loadBalancer(x, limit, threads);
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
        S2_total += phi_total[j] * (int128_t) mu_sum[i][j];
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
template <typename P, typename F>
int128_t S2(int128_t x,
            int64_t y,
            int64_t z,
            int64_t c,
            vector<P>& primes,
            FactorTable<F>& factors,
            int threads)
{
  int128_t S2_total = 0;

  threads = validate_threads(threads, z);
  PiTable pi(y);

  S2_total += S2_trivial(x, y, z, c, pi, primes, threads);
  S2_total += S2_easy(x, y, z, c, pi, primes, threads);
  S2_total += S2_sieve(x, y, z, c, pi, primes, factors, threads);

  return S2_total;
}

/// alpha is a tuning factor which should grow like (log(x))^3
/// for the Deleglise-Rivat prime counting algorithm.
///
double compute_alpha(int128_t x)
{
  double d = (double) x;
  double alpha = log(d) * log(d) * log(d) / 1000;
  return in_between(1, alpha, iroot<6>(x));
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Deleglise-Rivat algorithm.
/// Run time: O(x^(2/3) / (log x)^2) operations, O(x^(1/3) * (log x)^3) space.
///
int128_t pi_deleglise_rivat_parallel3(int128_t x, int threads)
{
  if (x < 2)
    return 0;

  if (x > to_maxint(primecount::max()))
    throw primecount_error("pi(x): x must be <= " + max());

  double alpha = compute_alpha(x);
  int64_t y = (int64_t) (alpha * iroot<3>(x));
  int64_t z = (int64_t) (x / y);
  int64_t pi_y;
  int64_t c;

  if (print_status())
  {
    cout << endl;
    cout << "=== pi_deleglise_rivat_parallel3(x) ===" << endl;
    cout << "pi(x) = S1 + S2 + pi(y) - 1 - P2" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
    cout << "c = " << PhiTiny::max_a() << endl;
    cout << "threads = " << validate_threads(threads) << endl;
  }

  int128_t p2 = P2(x, y, threads);
  int128_t s1, s2;

  if (y <= FactorTable<uint16_t>::max())
  {
    // if y < 2^32 we can use 32-bit primes and a 16-bit FactorTable
    // which uses ~ (y / 2) bytes of memory

    vector<uint32_t> primes = generate_primes<uint32_t>(y);
    FactorTable<uint16_t> factors(y);

    pi_y = primes.size() - 1;
    c = min(pi_y, PhiTiny::max_a());
    s1 = S1(x, y, c, primes[c], factors, threads);
    s2 = S2(x, y, z, c, primes, factors, threads);
  }
  else
  {
    // if y >= 2^32 we need to use 64-bit primes and a 32-bit
    // FactorTable which uses ~ y bytes of memory

    vector<int64_t> primes = generate_primes<int64_t>(y);
    FactorTable<uint32_t> factors(y);

    pi_y = primes.size() - 1;
    c = min(pi_y, PhiTiny::max_a());
    s1 = S1(x, y, c, primes[c], factors, threads);
    s2 = S2(x, y, z, c, primes, factors, threads);
  }

  int128_t phi = s1 + s2;
  int128_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace

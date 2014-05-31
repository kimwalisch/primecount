///
/// @file  pi_lmo_parallel4.cpp
/// @brief Parallel implementation of the Lagarias-Miller-Odlyzko
///        prime counting algorithm using OpenMP.
///
///        In this version the special leaves for c + 1 <= b < pi_sqrty
///        have been split up into 2 categories:
///           1) The special leaves that are a product of 2 primes.
///           2) The special leaves that are a product of a prime and a
///              square free integer (which must not be prime).
///        Although this split up can give up to 15 percent speed
///        improvement it uses considerably more memory.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <pmath.hpp>
#include <PhiTiny.hpp>
#include <tos_counters.hpp>
#include <utils.hpp>

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

/// For each prime calculate its first multiple >= low.
template <typename T1, typename T2>
void initialize_next_multiples(T1* next, T2& primes, int64_t size, int64_t low)
{
  next->reserve(size);
  next->push_back(0);

  for (int64_t b = 1; b < size; b++)
  {
    int64_t prime = primes[b];
    int64_t next_multiple = ((low + prime - 1) / prime) * prime;
    next->push_back(next_multiple);
  }
}

template <typename T1, typename T2>
void cross_off(int64_t prime, int64_t low, int64_t high, int64_t& next_multiple, T1& sieve, T2& counters)
{
  int64_t segment_size = sieve.size();
  int64_t k = next_multiple + prime * (~next_multiple & 1);

  for (; k < high; k += prime * 2)
  {
    if (sieve[k - low])
    {
      sieve[k - low] = 0;
      cnt_update(counters, k - low, segment_size);
    }
  }
  next_multiple = k;
}

void init_square_free_iters(vector<vector<int32_t>::iterator >* square_free_iters,
                            vector<vector<int32_t> >& square_free_candidates,
                            vector<int32_t>& primes,
                            int64_t c,
                            int64_t x,
                            int64_t y,
                            int64_t low)
{
  for (size_t i = c + 1; i < square_free_iters->size(); i++)
  {
    int64_t max_m = min(x / (primes[i] * low), y);
    (*square_free_iters)[i] = upper_bound(square_free_candidates[i].begin(),
        square_free_candidates[i].end(),
        max_m) - 1;
  }
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
                  int64_t pi_sqrty,
                  int64_t pi_y,
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
                  vector<int64_t>& phi,
                  vector<vector<int32_t> >& square_free_candidates)
{
  low += segment_size * segments_per_thread * thread_num;
  limit = min(low + segment_size * segments_per_thread, limit);
  int64_t size = pi[min(isqrt(x / low), y)] + 1;
  int64_t S2_thread = 0;

  if (c >= size - 1)
    return 0;

  vector<char> sieve(segment_size);
  vector<int32_t> counters(segment_size);
  vector<int64_t> next;
  vector<vector<int32_t>::iterator > square_free_iters(pi_sqrty);
  init_square_free_iters(&square_free_iters, square_free_candidates, primes, c, x, y, low);
  initialize_next_multiples(&next, primes, size, low);
  phi.resize(size, 0);
  mu_sum.resize(size, 0);

  // Process the segments corresponding to the current thread
  for (; low < limit; low += segment_size)
  {
    fill(sieve.begin(), sieve.end(), 1);

    // Current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t b = 1;

    for (; b <= c; b++)
    {
      int64_t k = next[b];
      for (int64_t prime = primes[b]; k < high; k += prime)
        sieve[k - low] = 0;
      next[b] = k;
    }

    // Initialize special tree data structure from sieve
    cnt_finit(sieve, counters, segment_size);

    // For c + 1 <= b < pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy:  mu[m] != 0 && primes[b] < lpf[m], low <= (x / n) < high
    for (; b < pi_sqrty; b++)
    {
      int64_t prime = primes[b];
      int64_t min_m = max(x / (prime * high), y / prime);
      int64_t l = pi[min(x / (prime * low), y)];
      min_m = in_between(prime, min_m, y);
      int64_t l_min = pi[min_m];
      vector<int32_t>::iterator iter = square_free_iters[b];

      if (prime >= primes[l])
        goto next_segment;

      // Special leaves which are a product of 2 primes
      for (; l > l_min; l--)
      {
        int64_t n = prime * primes[l];
        int64_t count = cnt_query(counters, (x / n) - low);
        int64_t phi_xn = phi[b] + count;
        S2_thread += phi_xn;
        mu_sum[b]++;
      }

      // Special leaves which are a product of a prime and a
      // square_free integer which must satisfy:
      // !is_prime(square_free) && prime < least_prime_factor[square_free]
      for (; *iter > min_m; iter--)
      {
        int64_t square_free = *iter;
        int64_t n = prime * square_free;
        int64_t count = cnt_query(counters, (x / n) - low);
        int64_t phi_xn = phi[b] + count;
        S2_thread -= mu[square_free] * phi_xn;
        mu_sum[b] -= mu[square_free];
      }

      square_free_iters[b] = iter;
      phi[b] += cnt_query(counters, (high - 1) - low);
      cross_off(prime, low, high, next[b], sieve, counters);
    }

    // For pi_sqrty <= b < pi_y
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high
    for (; b < min(pi_y, size); b++)
    {
      int64_t prime = primes[b];
      int64_t l = pi[min(x / (prime * low), y)];
      int64_t min_m = max(x / (prime * high), y / prime);
      min_m = in_between(prime, min_m, y);
      int64_t min_l = pi[min_m];

      if (prime >= primes[l])
        goto next_segment;

      for (; l > min_l; l--)
      {
        int64_t n = prime * primes[l];
        int64_t count = cnt_query(counters, (x / n) - low);
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
           int64_t pi_y,
           int64_t c,
           vector<int32_t>& primes,
           vector<int32_t>& lpf,
           vector<int32_t>& mu,
           int threads)
{
  threads = validate_threads(threads);

  int64_t S2_total = 0;
  int64_t low = 1;
  int64_t limit = x / y + 1;
  int64_t sqrt_limit = isqrt(limit);
  int64_t logx = max(1, ilog(x));
  int64_t min_segment_size = 1 << 6;
  int64_t segment_size = next_power_of_2(sqrt_limit / (logx * threads));
  segment_size = max(segment_size, min_segment_size);
  int64_t segments_per_thread = 1;
  int64_t pi_sqrty = pi_bsearch(primes, isqrt(y));

  vector<int32_t> pi = make_pi(y);
  vector<int64_t> phi_total(primes.size(), 0);
  vector<vector<int32_t> > square_free_candidates = 
      generate_square_free_candidates(c, y, lpf, mu, pi, primes);

  while (low < limit)
  {
    int64_t segments = (limit - low + segment_size - 1) / segment_size;
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, (segments + threads - 1) / threads);
    double seconds = get_wtime();

    vector<vector<int64_t> > phi(threads);
    vector<vector<int64_t> > mu_sum(threads);

    #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
    for (int i = 0; i < threads; i++)
      S2_total += S2_thread(x, y, c, pi_sqrty, pi_y, segment_size, segments_per_thread,
          i, low, limit, pi, primes, lpf, mu, mu_sum[i], phi[i],
              square_free_candidates);

    seconds = get_wtime() - seconds;
    low += segments_per_thread * threads * segment_size;

    // Dynamically increase segment_size or segments_per_thread
    // if the running time is less than a certain threshold.
    // We start off with a small segment size and few segments
    // per thread as most special leaves are in the first segments
    // whereas later on there are very few special leaves.
    //
    if (low > sqrt_limit && seconds < 10)
    {
      if (segment_size < sqrt_limit)
        segment_size <<= 1;
      else
        segments_per_thread *= 2;
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
  }

  return S2_total;
}

} // namespace

namespace primecount {

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3) * log log x) space.
///
int64_t pi_lmo_parallel4(int64_t x, int threads)
{
  if (x < 2)
    return 0;

  double beta = 1.0;
  double alpha = in_between(1, log(log((double) x)) * beta, iroot<6>(x));
  int64_t x13 = iroot<3>(x);
  int64_t y = (int64_t)(x13 * alpha);

  vector<int32_t> mu = make_moebius(y);
  vector<int32_t> lpf = make_least_prime_factor(y);
  vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(y, &primes);

  int64_t pi_y = primes.size() - 1;
  int64_t c = min<int64_t>(PhiTiny::MAX_A, pi_y);
  int64_t s1 = S1(x, y, c, primes, lpf , mu);
  int64_t s2 = S2(x, y, pi_y, c, primes, lpf , mu, threads);
  int64_t p2 = P2(x, y, threads);
  int64_t phi = s1 + s2;
  int64_t sum = phi + pi_y - 1 - p2;

  return sum;
}

} // namespace primecount

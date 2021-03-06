///
/// @file  B.cpp
/// @brief The B formula is a partial computation of the P2(x, a)
///        formula from the Lagarias-Miller-Odlyzko and Deleglise-Rivat
///        prime counting algorithms. P2(x, a) counts the numbers <= x
///        that have exactly 2 prime factors each exceeding the a-th
///        prime. Both P2 and B have a runtime complexity of
///        O(z log log z) and use O(z^(1/2)) memory, with z = x / y.
///
///        B(x, y) formula:
///        \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <int128_t.hpp>
#include <LoadBalancerP2.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace primecount;

namespace {

template <typename T>
struct ThreadResult
{
  T sum;
  int64_t pix;
  int64_t iters;
};

/// Count primes inside [prime, stop]
int64_t count_primes(primesieve::iterator& it,
                     int64_t& prime,
                     int64_t stop)
{
  int64_t count = 0;
  int64_t p = prime;

  for (; p <= stop; count++)
    p = it.next_prime();

  prime = p;
  return count;
}

template <typename T>
ThreadResult<T>
B_thread(T x,
         int64_t y,
         int64_t z,
         int64_t low,
         int64_t thread_num,
         int64_t thread_dist)
{
  T sum = 0;
  int64_t pix = 0;
  int64_t iters = 0;
  low += thread_dist * thread_num;

  if (low < z)
  {
    // thread sieves [low, z[
    z = min(low + thread_dist, z);
    int64_t start = (int64_t) max(x / z, y);
    int64_t stop = (int64_t) min(x / low, isqrt(x));

    primesieve::iterator it(low - 1, z);
    primesieve::iterator rit(stop + 1, start);
    int64_t next = it.next_prime();
    int64_t prime = rit.prev_prime();

    // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
    for (; prime > start; iters++)
    {
      int64_t xp = (int64_t)(x / prime);
      pix += count_primes(it, next, xp);
      prime = rit.prev_prime();
      sum += pix;
    }

    // Count the remaining primes
    pix += count_primes(it, next, z - 1);
  }

  return { sum, pix, iters };
}

/// \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
/// Run time: O(z log log z)
/// Memory usage: O(z^(1/2))
///
template <typename T>
T B_OpenMP(T x, int64_t y, int threads)
{
  if (x < 4)
    return 0;

  T sum = 0;
  int64_t low = 2;
  int64_t pi_low_minus_1 = 0;
  int64_t z = (int64_t)(x / max(y, 1));
  LoadBalancerP2 loadBalancer(z, threads);
  threads = loadBalancer.get_threads();
  aligned_vector<ThreadResult<T>> res(threads);

  while (low < z)
  {
    int64_t thread_dist = loadBalancer.get_thread_dist(low);

    #pragma omp parallel for num_threads(threads)
    for (int64_t i = 0; i < threads; i++)
      res[i] = B_thread(x, y, z, low, i, thread_dist);

    // The threads above have computed the sum of:
    // PrimePi(n) - PrimePi(thread_low - 1)
    // for many different values of n. However we actually want to
    // compute the sum of PrimePi(n). In order to get the complete
    // sum we now have to calculate the missing sum contributions in
    // sequential order as each thread depends on values from the
    // previous thread. The missing sum contribution for each thread
    // can be calculated using pi_low_minus_1 * iters.
    for (int64_t i = 0; i < threads; i++)
    {
      T thread_sum = res[i].sum;
      thread_sum += (T) pi_low_minus_1 * res[i].iters;
      sum += thread_sum;
      pi_low_minus_1 += res[i].pix;
    }

    low += thread_dist * threads;

    if (is_print())
    {
      int precision = get_status_precision(x);
      cout << "\rStatus: " << fixed << setprecision(precision)
           << get_percent(low, z) << '%' << flush;
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t B(int64_t x, int64_t y, int threads)
{
#ifdef ENABLE_MPI
  if (mpi_num_procs() > 1)
    return B_mpi(x, y, threads);
#endif

  print("");
  print("=== B(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int64_t sum = B_OpenMP((uint64_t) x, y, threads);

  print("B", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t B(int128_t x, int64_t y, int threads)
{
#ifdef ENABLE_MPI
  if (mpi_num_procs() > 1)
    return B_mpi(x, y, threads);
#endif

  print("");
  print("=== B(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int128_t sum = B_OpenMP((uint128_t) x, y, threads);

  print("B", sum, time);
  return sum;
}

#endif

} // namespace

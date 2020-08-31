///
/// @file  B_mpi.cpp
/// @brief Implementation of the B formula (from Xavier Gourdon's
///        algorithm) that has been distributed using MPI (Message
///        Passing Interface) and multi-threaded using OpenMP.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <mpi_reduce_sum.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace primecount;

namespace {

/// Calculate the thread sieving distance. The idea
/// is to gradually increase the thread_dist in
/// order to keep all CPU cores busy.
///
class LoadBalancer
{
public:
  LoadBalancer(int64_t z)
    : z_(z)
  { }

  void get_work(int64_t low,
                int* threads,
                int64_t* thread_dist)
  {
    double start_time = time_;
    time_ = get_time();
    double seconds = time_ - start_time;

    if (start_time > 0)
    {
      if (seconds < 60)
        thread_dist_ *= 2;
      if (seconds > 60)
        thread_dist_ /= 2;
    }

    int64_t dist = z_ - min(low, z_);
    int64_t max_dist = ceil_div(dist, *threads);
    thread_dist_ = in_between(min_dist_, thread_dist_, max_dist);
    *threads = ideal_num_threads(*threads, dist, thread_dist_);
    *thread_dist = thread_dist_;
  }

private:
  double time_ = -1;
  int64_t min_dist_ = 1 << 20;
  int64_t thread_dist_ = min_dist_;
  int64_t z_;
};

/// Count the primes inside [prime, stop]
int64_t count_primes(primesieve::iterator& it, int64_t& prime, int64_t stop)
{
  int64_t count = 0;

  for (; prime <= stop; count++)
    prime = it.next_prime();

  return count;
}

template <typename T>
struct ThreadResult
{
  T sum;
  int64_t pix;
  int64_t iters;
};

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

    primesieve::iterator rit(stop + 1, start);
    primesieve::iterator it(low - 1, z);
    int64_t next = it.next_prime();
    int64_t prime = rit.prev_prime();

    // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
    while (prime > start)
    {
      int64_t xp = (int64_t)(x / prime);
      if (xp >= z) break;
      pix += count_primes(it, next, xp);
      iters++;
      sum += pix;
      prime = rit.prev_prime();
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
  int proc_id = mpi_proc_id();
  int procs = mpi_num_procs();

  int64_t low = 2;
  int64_t thread_dist = 0;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t proc_dist = ceil_div(z, procs);
  low += proc_dist * proc_id;
  int64_t pi_low_minus_1 = pi_simple(low - 1, threads);
  z = min(low + proc_dist, z);

  LoadBalancer loadBalancer(z);
  loadBalancer.get_work(low, &threads, &thread_dist);
  aligned_vector<ThreadResult<T>> res(threads);

  while (low < z)
  {
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
    loadBalancer.get_work(low, &threads, &thread_dist);

    if (is_print())
    {
      double percent = get_percent(low, z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
          << percent << '%' << flush;
    }
  }

  sum = mpi_reduce_sum(sum);

  return sum;
}

} // namespace

namespace primecount {

int64_t B_mpi(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== B_mpi(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int64_t sum = B_OpenMP((uint64_t) x, y, threads);

  print("B", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t B_mpi(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== B_mpi(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int128_t sum = B_OpenMP((uint128_t) x, y, threads);

  print("B", sum, time);
  return sum;
}

#endif

} // namespace

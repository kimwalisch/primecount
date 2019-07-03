///
/// @file  P2_mpi.cpp
/// @brief 2nd partial sieve function. P2(x, y) counts the
///        numbers <= x that have exactly 2 prime factors
///        each exceeding the a-th prime.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
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

int64_t count_primes(primesieve::iterator& it, int64_t& prime, int64_t stop)
{
  int64_t count = 0;

  for (; prime <= stop; count++)
    prime = it.next_prime();

  return count;
}

/// Calculate the thread sieving distance. The idea is to
/// gradually increase the thread_distance in order to
/// keep all CPU cores busy.
///
void balanceLoad(int64_t* thread_distance, 
                 int64_t low,
                 int64_t z,
                 int threads,
                 double start_time)
{
  double seconds = get_time() - start_time;

  int64_t min_distance = 1 << 23;
  int64_t max_distance = ceil_div(z - low, threads);

  if (seconds < 60)
    *thread_distance *= 2;
  if (seconds > 60)
    *thread_distance /= 2;

  *thread_distance = in_between(min_distance, *thread_distance, max_distance);
}

template <typename T>
T P2_thread(T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t thread_num,
            int64_t thread_distance,
            int64_t& pix,
            int64_t& pix_count)
{
  T p2 = 0;
  pix = 0;
  pix_count = 0;

  low += thread_distance * thread_num;
  z = min(low + thread_distance, z);
  int64_t start = (int64_t) max(x / z, y);
  int64_t stop = (int64_t) min(x / low, isqrt(x));

  primesieve::iterator rit(stop + 1, start);
  primesieve::iterator it(low - 1, z);

  int64_t next = it.next_prime();
  int64_t prime = rit.prev_prime();

  // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i])
  while (prime > start)
  {
    int64_t xp = (int64_t)(x / prime);
    if (xp >= z) break;
    pix += count_primes(it, next, xp);
    pix_count++;
    p2 += pix;
    prime = rit.prev_prime();
  }

  pix += count_primes(it, next, z - 1);

  return p2;
}

/// P2(x, y) counts the numbers <= x that have exactly 2
/// prime factors each exceeding the a-th prime.
/// Memory usage: O(z^(1/2))
///
template <typename T>
T P2_mpi_master(T x, int64_t y, int threads)
{
  static_assert(prt::is_signed<T>::value,
                "P2(T x, ...): T must be signed integer type");

  if (x < 4)
    return 0;

  T a = pi_legendre(y, threads);
  T b = pi_legendre((int64_t) isqrt(x), threads);

  if (a >= b)
    return 0;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;

  int proc_id = mpi_proc_id();
  int procs = mpi_num_procs();

  int64_t distance = z - low;
  int64_t proc_distance = ceil_div(distance, procs);
  low += proc_distance * proc_id;
  z = min(low + proc_distance, z);

  T p2 = 0;
  T pix_total = pi_legendre(low - 1, threads);

  if (is_mpi_master_proc())
    p2 = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;

  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double time = get_time();

    #pragma omp parallel for num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_thread(x, y, z, low, i, thread_distance, pix[i], pix_counts[i]);

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, time);

    // add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      p2 += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    if (is_print())
    {
      double percent = get_percent(low, z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  p2 = mpi_reduce_sum(p2);

  return p2;
}

} // namespace

namespace primecount {

int64_t P2_mpi(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== P2_mpi(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_time();
  int64_t p2 = P2_mpi_master(x, y, threads);

  print("P2", p2, time);
  return p2;
}

#ifdef HAVE_INT128_T

int128_t P2_mpi(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== P2_mpi(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_time();
  int128_t p2 = P2_mpi_master(x, y, threads);

  print("P2", p2, time);
  return p2;
}

#endif

} // namespace

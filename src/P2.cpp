///
/// @file  P2.cpp
/// @brief 2nd partial sieve function.
///        P2(x, y) counts the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Calculate the segment_size per thread.
/// The idea is to gradually increase the segment_size per
/// thread in order to keep all CPU cores busy.
///
int64_t balanceLoad(int64_t segment_size, double start_time)
{
  double seconds = get_wtime() - start_time;

  if (seconds < 30)
    segment_size *= 2;
  else if (segment_size >= 4)
    segment_size -= segment_size / 4;

  return segment_size;
}

template <typename T>
int64_t count_primes(primesieve::iterator& it,
                     int64_t& prime,
                     T stop)
{
  int64_t primes = 0;

  for (; prime <= stop; primes++)
    prime = it.next_prime();

  return primes;
}

template <typename T>
T P2_OpenMP_thread(T x,
                   int64_t y,
                   int64_t z,
                   int64_t segment_size,
                   int64_t thread_num,
                   int64_t low,
                   int64_t& pix,
                   int64_t& pix_count)
{
  pix = 0;
  pix_count = 0;
  low += thread_num * segment_size;
  z = min(low + segment_size, z);
  int64_t start = (int64_t) max(x / z + 1, y);
  int64_t stop = (int64_t) min(x / low, isqrt(x));

  primesieve::iterator rit(stop + 1, start);
  primesieve::iterator it(low - 1, z);

  int64_t next = it.next_prime();
  int64_t prime = rit.previous_prime();
  T P2_thread = 0;

  // \sum_{i = pi[start]}^{pi[stop]} pi(x / primes[i])
  while (prime >= start &&
         x / prime < z)
  {
    pix += count_primes(it, next, x / prime);
    P2_thread += pix;
    pix_count++;
    prime = rit.previous_prime();
  }

  pix += count_primes(it, next, z - 1);

  return P2_thread;
}

/// P2(x, y) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime, a = pi(y).
/// Space complexity: O((x / y)^(1/2)).
///
template <typename T>
T P2_OpenMP_master(T x, int64_t y, int threads)
{
#if __cplusplus >= 201103L
  static_assert(prt::is_signed<T>::value,
                "P2(T x, ...): T must be signed integer type");
#endif

  if (x < 4)
    return 0;

  T a = pi_legendre(y, threads);
  T b = pi_legendre((int64_t) isqrt(x), threads);

  if (a >= b)
    return 0;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t segment_size = 1 << 20;
  threads = validate_threads(threads, z);

  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  // \sum_{i=a+1}^{b} pi(x / primes[i]) - (i - 1)
  T p2 = 0;
  T pix_total = 0;

  // \sum_{i=a+1}^{b} -(i - 1)
  p2 = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < z)
  {
    int64_t segments = ceil_div(z - low, segment_size);
    threads = in_between(1, threads, segments);
    double time = get_wtime();

    #pragma omp parallel for num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_OpenMP_thread(x, y, z, segment_size,
        i, low, pix[i], pix_counts[i]);

    low += segment_size * threads;
    segment_size = balanceLoad(segment_size, time);

    // Add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      p2 += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    if (print_status())
    {
      double percent = get_percent((double) low, (double) z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  return p2;
}

} // namespace

namespace primecount {

int64_t P2(int64_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int64_t p2 = P2_OpenMP_master(x, y, threads);

  print("P2", p2, time);
  return p2;
}

#ifdef HAVE_INT128_T

int128_t P2(int128_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int128_t p2 = P2_OpenMP_master(x, y, threads);

  print("P2", p2, time);
  return p2;
}

#endif

} // namespace primecount

///
/// @file  B.cpp
/// @brief The B formula is a partial computation of the P2(x, a)
///        formula from the Lagarias-Miller-Odlyzko and Deleglise-
///        Rivat prime counting algorithms. P2(x, a) counts the
///        numbers <= x that have exactly 2 prime factors each
///        exceeding the a-th prime. Both P2 and B have a runtime
///        complexity of O(z log log z) and use O(z^(1/2)) memory,
///        with z = x / y. This implementation is a simplified
///        version of P2.cpp.
///
///        B(x, y) formula:
///        \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <tuple>

using namespace std;
using namespace primecount;

namespace {

/// Count the primes inside [prime, stop]
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
std::tuple<T, int64_t, int64_t>
B_thread(T x,
         int64_t y,
         int64_t z,
         int64_t low,
         int64_t thread_num,
         int64_t thread_distance)
{
  T sum = 0;
  int64_t pix = 0;
  int64_t pix_count = 0;

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
    sum += pix;
    prime = rit.prev_prime();
  }

  pix += count_primes(it, next, z - 1);
  auto res = make_tuple(sum, pix, pix_count);

  return res;
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
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;
  int64_t pix_low = 0;

  // prevents CPU false sharing
  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);
  
  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double time = get_time();

    #pragma omp parallel for num_threads(threads) reduction(+: sum)
    for (int i = 0; i < threads; i++)
    {
      auto res = B_thread(x, y, z, low, i, thread_distance);
      pix[i] = std::get<1>(res);
      pix_counts[i] = std::get<2>(res);
      sum += std::get<0>(res);
    }

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, time);

    // Add the missing sum contributions in sequential order.
    // The sum from the parallel for loop above is the sum
    // of the prime counts inside [thread_start, thread_stop].
    // For each such prime count we now have to add the
    // missing part from [0, thread_start - 1].
    for (int i = 0; i < threads; i++)
    {
      T count = pix_counts[i];
      sum += pix_low * count;
      pix_low += pix[i];
    }

    if (is_print())
    {
      double percent = get_percent(low, z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t B(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== B(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int64_t sum = B_OpenMP((intfast64_t) x, y, threads);

  print("B", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t B(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== B(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int128_t sum = B_OpenMP((intfast128_t) x, y, threads);

  print("B", sum, time);
  return sum;
}

#endif

} // namespace

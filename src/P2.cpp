///
/// @file  P2.cpp
/// @brief P2(x, a) is the 2nd partial sieve function.
///        P2(x, a) counts the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime. This implementation
///        uses the primesieve library for quickly iterating over
///        primes using next_prime() and prev_prime() which greatly
///        simplifies the implementation.
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
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
#include <noinline.hpp>
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
NOINLINE std::tuple<T, int64_t, int64_t>
P2_thread(T x,
          int64_t y,
          int64_t z,
          int64_t low,
          int64_t thread_num,
          int64_t thread_distance)
{
  T sum = 0;
  int64_t pix = 0;
  int64_t pix_count = 0;

  // thread sieves [low, z[
  low += thread_distance * thread_num;
  z = min(low + thread_distance, z);
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
    pix_count++;
    sum += pix;
    prime = rit.prev_prime();
  }

  // prime count [low, z[
  pix += count_primes(it, next, z - 1);
  auto res = make_tuple(sum, pix, pix_count);

  return res;
}

/// P2(x, y) counts the numbers <= x that have exactly 2
/// prime factors each exceeding the a-th prime.
/// Run-time: O(z log log z)
///
template <typename T>
T P2_OpenMP(T x, int64_t y, int threads)
{
  static_assert(prt::is_signed<T>::value,
                "T must be signed integer type");

  if (x < 4)
    return 0;

  T a = pi_simple(y, threads);
  T b = pi_simple((int64_t) isqrt(x), threads);

  if (a >= b)
    return 0;

  // \sum_{i=a+1}^{b} -(i - 1)
  T sum = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;
  int64_t pi_low_minus_1 = 0;

  // prevents CPU false sharing
  using res_t = std::tuple<T, int64_t, int64_t>;
  aligned_vector<res_t> res(threads);

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double time = get_time();

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++)
      res[i] = P2_thread(x, y, z, low, i, thread_distance);

    // The threads above have computed the sum of:
    // PrimePi(n) - PrimePi(thread_low - 1)
    // for many different values of n. However we actually want to
    // compute the sum of PrimePi(n). In order to get the complete
    // sum we now have to calculate the missing sum contributions in
    // sequential order as each thread depends on values from the
    // previous thread. The missing sum contribution for each thread
    // can be calculated using pi_low_minus_1 * thread_count.
    for (int i = 0; i < threads; i++)
    {
      auto thread_sum = std::get<0>(res[i]);
      auto thread_pix = std::get<1>(res[i]);
      auto thread_count = std::get<2>(res[i]);
      thread_sum += (T) pi_low_minus_1 * thread_count;
      sum += thread_sum;
      pi_low_minus_1 += thread_pix;
    }

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, time);

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

int64_t P2(int64_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print_vars(x, y, threads);

  double time = get_time();
  int64_t sum = P2_OpenMP(x, y, threads);

  print("P2", sum, time);
  return sum;
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
  print_vars(x, y, threads);

  double time = get_time();
  int128_t sum = P2_OpenMP(x, y, threads);

  print("P2", sum, time);
  return sum;
}

#endif

} // namespace

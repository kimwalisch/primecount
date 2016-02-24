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
#include <BitSieve.hpp>
#include <generate.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <Wheel.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Calculate the segments per thread.
/// The idea is to gradually increase the segments per thread (based
/// on elapsed time) in order to keep all CPU cores busy.
///
int64_t balanceLoad(int64_t segments_per_thread, double start_time)
{
  double seconds = get_wtime() - start_time;

  if (seconds < 30)
    segments_per_thread *= 2;
  else if (segments_per_thread >= 4)
    segments_per_thread -= segments_per_thread / 4;

  return segments_per_thread;
}

/// Cross-off the multiples inside [low, high[
/// of the primes <= sqrt(high - 1).
///
void cross_off(BitSieve& sieve,
               vector<int32_t> primes,
               Wheel& wheel,
               int64_t c,
               int64_t low,
               int64_t high)
{
  int64_t pi_sqrt_high = pi_bsearch(primes, isqrt(high - 1));

  for (int64_t i = c + 1; i <= pi_sqrt_high; i++)
  {
    int64_t prime = primes[i];
    int64_t m = wheel[i].next_multiple;
    int64_t wheel_index = wheel[i].wheel_index;

    // cross-off the multiples of prime inside [low, high[
    for (; m < high; m += prime * Wheel::next_multiple_factor(&wheel_index))
      sieve.unset(m - low);

    wheel[i].set(m, wheel_index);
  }
}

template <typename T>
T P2_OpenMP_thread(T x,
                   int64_t y,
                   int64_t segment_size,
                   int64_t segments_per_thread,
                   int64_t thread_num,
                   int64_t low,
                   int64_t limit,
                   int64_t& pix,
                   int64_t& pix_count,
                   vector<int32_t>& primes)
{
  pix = 0;
  pix_count = 0;
  low += thread_num * segments_per_thread * segment_size;
  limit = min(low + segments_per_thread * segment_size, limit);
  int64_t size = pi_bsearch(primes, isqrt(limit)) + 1;
  int64_t start = (int64_t) max(x / limit + 1, y);
  int64_t stop  = (int64_t) min(x / low, isqrt(x));
  T P2_thread = 0;

  // P2_thread = \sum_{i = pi[start]}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
  // We use a reverse prime iterator to calculate P2_thread
  primesieve::iterator pi(stop + 1, start);
  int64_t prime = pi.previous_prime();

  bool sieve_primes = true;
  Wheel wheel(primes, size, low, sieve_primes);
  BitSieve sieve(segment_size);

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t x_div_prime = 0;
    int64_t j = 0;
    int64_t c = 6;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low, sieve_primes);

    // cross-off the multiples of the primes <= sqrt(high - 1)
    cross_off(sieve, primes, wheel, c, low, high);

    while (prime >= start &&
           (x_div_prime = (int64_t) (x / prime)) < high)
    {
      int64_t next_count = x_div_prime - low;
      pix += sieve.count(j, next_count);
      j = next_count + 1;
      pix_count++;
      P2_thread += pix;
      prime = pi.previous_prime();
    }

    pix += sieve.count(j, (high - 1) - low);
  }

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
  int64_t limit = (int64_t)(x / max(y, 1));
  int64_t sqrt_limit = isqrt(limit);
  int64_t segment_size = max(sqrt_limit, 1 << 12);
  int64_t segments_per_thread = 64;
  threads = validate_threads(threads, limit);

  vector<int32_t> primes = generate_primes(sqrt_limit);
  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  // \sum_{i=a+1}^{b} pi(x / primes[i]) - (i - 1)
  T p2 = 0;
  T pix_total = 0;

  // \sum_{i=a+1}^{b} -(i - 1)
  p2 = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < limit)
  {
    int64_t segments = ceil_div(limit - low, segment_size);
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, ceil_div(segments, threads));
    double time = get_wtime();

    #pragma omp parallel for \
        num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_OpenMP_thread(x, y, segment_size, segments_per_thread,
         i, low, limit, pix[i], pix_counts[i], primes);

    low += segments_per_thread * threads * segment_size;
    segments_per_thread = balanceLoad(segments_per_thread, time);

    // Add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      p2 += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    if (print_status())
    {
      double percent = get_percent((double) low, (double) limit);
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

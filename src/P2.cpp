///
/// @file  P2.cpp
/// @brief 2nd partial sieve function.
///        P2(x, y) counts the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
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
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace P2 {

class ReversePrimeIterator
{
public:
  ReversePrimeIterator(int64_t stop, int64_t start) :
    iter_(stop, start),
    prime_(stop)
  { }
  int64_t previous_prime()
  {
    if (prime_ <= 2)
      return -1;
    prime_ = iter_.previous_prime();
    return prime_;
  }
private:
  primesieve::iterator iter_;
  int64_t prime_;
};

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
    int64_t k = wheel[i].next_multiple;
    int64_t wheel_index = wheel[i].wheel_index;

    // cross-off the multiples of prime inside [low, high[
    for (; k < high; k += prime * Wheel::next_multiple_factor(&wheel_index))
      sieve.unset(k - low);

    wheel[i].next_multiple = k;
    wheel[i].wheel_index = wheel_index;
  }
}

/// Calculate the segments per thread.
/// The idea is to gradually increase the segments per thread (based
/// on elapsed time) in order to keep all CPU cores busy. 
///
int64_t balanceLoad(int64_t segments_per_thread, double seconds1, double time1)
{
  double time2 = get_wtime();
  double seconds = time2 - seconds1;
  double time = time2 - time1;
  double increase_threshold = in_between(0.5, time / 10, 20);

  if (seconds < increase_threshold)
    segments_per_thread += segments_per_thread * 3;
  else if (segments_per_thread >= 4)
    segments_per_thread -= segments_per_thread / 4;

  return segments_per_thread;
}

template <typename T>
T P2_thread(T x,
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

  // P2_thread = \sum_{i=pi[start]}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
  // We use a reverse prime iterator to calculate P2_thread
  ReversePrimeIterator prime_iter(stop + 1, start);
  int64_t prime = prime_iter.previous_prime();
  int64_t xp = (int64_t) (x / prime);

  bool sieve_primes = true;
  Wheel wheel(primes, size, low, sieve_primes);
  BitSieve sieve(segment_size);

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t c = 6;
    int64_t j = 0;

    // pre-sieve the multiples of the first c primes
    sieve.pre_sieve(c, low, sieve_primes);

    // cross-off the multiples of the primes <= sqrt(high - 1)
    cross_off(sieve, primes, wheel, c, low, high);

    while (prime >= start && 
           xp < high)
    {
      pix += sieve.count(j, xp - low);
      j = xp - low + 1;
      pix_count++;
      P2_thread += pix;
      prime = prime_iter.previous_prime();
      xp = (int64_t) (x / prime);
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
T P2(T x, int64_t y, int threads)
{
#if __cplusplus >= 201103L
  static_assert(prt::is_signed<T>::value,
                "P2(T x, ...): T must be signed integer type");
#endif

  T a = pi_legendre(y, 1);
  T b = pi_legendre((int64_t) isqrt(x), 1);

  if (x < 4 || a >= b)
    return 0;

  int64_t low = 2;
  int64_t limit = (int64_t)(x / max(y, 1));
  int64_t segment_size = max(isqrt(limit), 1 << 12);
  int64_t segments_per_thread = 1;
  threads = validate_threads(threads, limit);

  vector<int32_t> primes = generate_primes(isqrt(limit));
  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);
  double time = get_wtime();

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
    double seconds = get_wtime();

    #pragma omp parallel for \
        num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_thread(x, y, segment_size, segments_per_thread, i,
         low, limit, pix[i], pix_counts[i], primes);

    low += segments_per_thread * threads * segment_size;
    segments_per_thread = balanceLoad(segments_per_thread, seconds, time);

    // Add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      p2 += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    if (print_status())
      cout << "\rStatus: " << get_percent(low, limit) << '%' << flush;
  }

  return p2;
}

} // namespace P2
} // namespace

namespace primecount {

int64_t P2(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int64_t p2 = P2::P2(x, y, threads);

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
  int128_t p2 = P2::P2(x, y, threads);

  print("P2", p2, time);
  return p2;
}

#endif

/// P2_lehmer(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime. This implementation is
/// optimized for small values of a < pi(x^(1/3)) which requires
/// sieving up to a large limit (x / primes[a]). Sieving is done in
/// parallel using primesieve (segmented sieve of Eratosthenes).
/// Space complexity: O(pi(sqrt(x))).
///
int64_t P2_lehmer(int64_t x, int64_t a, int threads)
{
  print("");
  print("=== P2_lehmer(x, a) ===");
  print("Computation of the 2nd partial sieve function");

  double time = get_wtime();
  vector<int32_t> primes = generate_primes(isqrt(x));
  vector<int64_t> counts(primes.size());

  int64_t b = pi_bsearch(primes, isqrt(x));
  int64_t p2 = 0;
  int64_t pix = 0;

  #pragma omp parallel for schedule(dynamic) \
      num_threads(validate_threads(threads, b, 1000))
  for (int64_t i = b; i > a; i--)
  {
    int64_t prev = (i == b) ? 0 : x / primes[i + 1] + 1;
    int64_t xi = x / primes[i];
    counts[i] = primesieve::count_primes(prev, xi);
  }

  for (int64_t i = b; i > a; i--)
  {
    pix += counts[i];
    p2 += pix - (i - 1);
  }

  print("P2", p2, time);
  return p2;
}

} // namespace primecount

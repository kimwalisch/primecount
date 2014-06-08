///
/// @file  P2.cpp
/// @brief 2nd partial sieve function.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <bit_sieve.hpp>
#include <pmath.hpp>
#include <utils.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// @return previous_prime or -1 if old <= 2
inline int64_t get_previous_prime(primesieve::iterator* iter, int64_t old)
{
  return (old > 2) ? iter->previous_prime() : -1;
}

/// For each prime calculate its first multiple >= low
/// which is not a multiple of 2.
void init_next_multiples(vector<int64_t>& next, vector<int32_t>& primes, int64_t size, int64_t low)
{
  next.reserve(size);
  next.push_back(0);

  for (int64_t b = 1; b < size; b++)
  {
    int64_t prime = primes[b];
    int64_t next_multiple = ((low + prime - 1) / prime) * prime;
    next_multiple += prime * (~next_multiple & 1);
    next.push_back(max(isquare(prime), next_multiple));
  }
}

int64_t P2_thread(int64_t x,
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
  int64_t sqrtx = isqrt(x);
  int64_t start = max(x / limit + 1, y);
  int64_t stop = min(x / low, sqrtx);
  int64_t P2_thread = 0;

  bit_sieve sieve(segment_size);
  vector<int64_t> pi_input;
  vector<int64_t> next;
  init_next_multiples(next, primes, size, low);

  // P2_thread = \sum_{i=pi[start]}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
  // We use a reverse prime iterator to calculate P2_thread.
  primesieve::iterator iter(stop + 1, start);
  int64_t previous_prime = get_previous_prime(&iter, stop + 1);
  int64_t xp = x / previous_prime;

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t sqrt = isqrt(high - 1);
    int64_t j = 0;

    sieve.fill(low);

    // cross-off multiples
    for (int64_t i = 2; i < size && primes[i] <= sqrt; i++)
    {
      int64_t k;
      int64_t p2 = primes[i] * 2;
      for (k = next[i]; k < high; k += p2)
        sieve.unset(k - low);
      next[i] = k;
    }

    while (previous_prime >= start && xp < high)
    {
      pix += sieve.count(j, xp - low);
      j = xp - low + 1;
      pix_count++;
      P2_thread += pix;
      previous_prime = get_previous_prime(&iter, previous_prime);
      xp = x / previous_prime;
    }

    pix += sieve.count(j, (high - 1) - low);
  }

  return P2_thread;
}

} // namespace

namespace primecount {

/// 2nd partial sieve function.
/// P2(x, y) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime, a = pi(y).
/// Space complexity: O((x / y)^(1/2)).
///
int64_t P2(int64_t x, int64_t y, int threads)
{
  int64_t a = pi_legendre(y, 1);
  int64_t b = pi_legendre(isqrt(x), 1);

  if (x < 4 || a >= b)
    return 0;

  threads = validate_threads(threads);

  // \sum_{i=a+1}^{b} pi(x / primes[i]) - (i - 1)
  // initialize with \sum_{i=a+1}^{b} -i + 1
  int64_t sum = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;
  int64_t low = 2;
  int64_t pix_total = 1; // We sieve using primes > 2
  int64_t limit = x / max<int64_t>(1, y);
  int64_t segment_size = max<int64_t>(64, isqrt(limit));
  int64_t segments_per_thread = 1;

  vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(isqrt(limit), &primes);
  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  while (low < limit)
  {
    int64_t segments = (limit - low + segment_size - 1) / segment_size;
    threads = in_between(1, threads, segments);
    segments_per_thread = in_between(1, segments_per_thread, (segments + threads - 1) / threads);
    double seconds = get_wtime();

    #pragma omp parallel for num_threads(threads) reduction(+: sum)
    for (int i = 0; i < threads; i++)
      sum += P2_thread(x, y, segment_size, segments_per_thread, i, low, limit, 
          pix[i], pix_counts[i], primes);

    seconds = get_wtime() - seconds;
    low += segments_per_thread * threads * segment_size;

    // Adjust thread load balancing
    if (seconds < 10)
      segments_per_thread *= 2;
    else if (seconds > 30 && segments_per_thread > 1)
      segments_per_thread /= 2;

    // Add the missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      sum += pix_total * pix_counts[i];
      pix_total += pix[i];
    }
  }

  return sum;
}

/// 2nd partial sieve function.
/// P2_lehmer(x, a) counts the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime. This implementation is
/// optimized for small values of a < pi(x^(1/3)) which requires
/// sieving up to a large limit (x / primes[a]). Sieving is done in
/// parallel using primesieve (segmented sieve of Eratosthenes).
/// Space complexity: O(pi(sqrt(x))).
///
int64_t P2_lehmer(int64_t x, int64_t a, int threads)
{
  vector<int32_t> primes;
  vector<int64_t> counts;
  primes.push_back(0);
  primesieve::generate_primes(isqrt(x), &primes);
  counts.resize(primes.size());

  int64_t b = pi_bsearch(primes, isqrt(x));
  int64_t sum = 0;
  int64_t pix = 0;

  #pragma omp parallel for num_threads(validate_threads(threads)) schedule(dynamic)
  for (int64_t i = b; i > a; i--)
  {
    int64_t prev = (i == b) ? 0 : x / primes[i + 1] + 1;
    int64_t xi = x / primes[i];
    counts[i] = primesieve::count_primes(prev, xi);
  }

  for (int64_t i = b; i > a; i--)
  {
    pix += counts[i];
    sum += pix - (i - 1);
  }

  return sum;
}

} // namespace primecount

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
#include <pi_bsearch.hpp>
#include <pmath.hpp>
#include <utils.hpp>

#include <stdint.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace {

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
  int64_t P2_thread = 0;
  int64_t start = max(x / limit + 1, y);
  int64_t stop = in_between(2, x / low, isqrt(x));

  vector<char> sieve(segment_size);
  vector<int64_t> pi_input;
  vector<int64_t> next;
  next.push_back(0);
  next.reserve(size);

  // P2_thread = \sum_{i=pi[start]}^{pi[stop]} pi(x / primes[i])
  primesieve::generate_primes(start, stop, &pi_input);

  // reverse iterator
  auto iter = pi_input.rbegin();
  auto rend = pi_input.rend();

  for (size_t i = 0; i < pi_input.size(); i++)
    pi_input[i] = x / pi_input[i];

  // initialize next multiples
  for (int64_t b = 1; b < size; b++)
  {
    int64_t prime = primes[b];
    int64_t next_multiple = ((low + prime - 1) / prime) * prime;
    next_multiple = max(isquare(prime), next_multiple + (~next_multiple & 1) * prime);
    next.push_back(next_multiple);
  }

  // segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    fill(sieve.begin(), sieve.end(), 1);

    // current segment = interval [low, high[
    int64_t high = min(low + segment_size, limit);
    int64_t sqrt = isqrt(high - 1);
    int64_t j = ~low & 1;

    // cross-off multiples
    for (int64_t i = 2; i < size && primes[i] <= sqrt; i++)
    {
      int64_t k;
      int64_t p2 = primes[i] * 2;
      for (k = next[i]; k < high; k += p2)
        sieve[k - low] = 0;
      next[i] = k;
    }

    for (; iter != rend && *iter < high; iter++)
    {
      for (int64_t xil = *iter - low; j <= xil; j += 2)
        pix += sieve[j];
      // P2_thread += pi(x / primes[i])
      P2_thread += pix;
      pix_count++;
    }

    for (; j < high - low; j += 2)
      pix += sieve[j];
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
  int64_t pix_total = 1;
  int64_t low = 3;
  int64_t limit = (y > 0) ? x / y : x;
  int64_t sqrt_limit = isqrt(limit);
  int64_t min_segment_size = 64;
  int64_t segment_size = max(min_segment_size, sqrt_limit);
  int64_t segments_per_thread = 1;

  vector<int32_t> primes;
  primes.push_back(0);
  primesieve::generate_primes(sqrt_limit, &primes);
  vector<int64_t> pix_counts(threads);
  vector<int64_t> pix(threads);

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

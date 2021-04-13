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
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>
#include <LoadBalancerP2.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <cassert>

using namespace std;
using namespace primecount;

namespace {

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
T B_thread(T x,
           int64_t y,
           int64_t low,
           int64_t high)
{
  assert(low > 0);
  assert(low < high);
  int threads = 1;
  int64_t pix = pi_noprint(low - 1, threads);

  // thread sieves [low, high[
  int64_t sqrtx = isqrt(x);
  int64_t start = max(y, min(x / high, sqrtx));
  int64_t stop = min(x / low, sqrtx);

  primesieve::iterator it(low - 1, high);
  primesieve::iterator rit(stop + 1, start);
  int64_t next = it.next_prime();
  int64_t prime = rit.prev_prime();
  T sum = 0;

  // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i])
  while (prime > start)
  {
    int64_t xp = (int64_t)(x / prime);
    pix += count_primes(it, next, xp);
    prime = rit.prev_prime();
    sum += pix;
  }

  return sum;
}

/// \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
/// Run time: O(z log log z)
/// Memory usage: O(z^(1/2))
///
template <typename T>
T B_OpenMP(T x,
           int64_t y,
           int threads,
           bool is_print)
{
  if (x < 4)
    return 0;

  T sum = 0;
  int64_t z = (int64_t)(x / max(y, 1));
  LoadBalancerP2 loadBalancer(x, z, threads, is_print);
  threads = loadBalancer.get_threads();

  // for (low = sqrt(x); low < z; low += dist)
  #pragma omp parallel num_threads(threads) reduction(+:sum)
  {
    int64_t low, high;
    while (loadBalancer.get_work(low, high))
      sum += B_thread(x, y, low, high);
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t B(int64_t x,
          int64_t y,
          int threads,
          bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== B(x, y) ===");
    print_gourdon_vars(x, y, threads);
  }

  double time = get_time();
  int64_t sum = B_OpenMP((uint64_t) x, y, threads, is_print);

  if (is_print)
    print("B", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t B(int128_t x,
           int64_t y,
           int threads,
           bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== B(x, y) ===");
    print_gourdon_vars(x, y, threads);
  }

  double time = get_time();
  int128_t sum = B_OpenMP((uint128_t) x, y, threads, is_print);

  if (is_print)
    print("B", sum, time);

  return sum;
}

#endif

} // namespace

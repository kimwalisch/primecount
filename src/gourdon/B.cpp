///
/// @file  B.cpp
/// @brief The B formula is a partial computation of the P2(x, a)
///        formula from the Lagarias-Miller-Odlyzko and Deleglise-Rivat
///        prime counting algorithms. P2(x, a) counts the numbers <= x
///        that have exactly 2 prime factors each exceeding the a-th
///        prime. Both P2 and B have a runtime complexity of
///        O(n log log n) and use O(n^(1/2)) memory, with n = x / y.
///
///        B(x, y) formula:
///        \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>
#include <LoadBalancerP2.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>

using namespace primecount;

namespace {

/// Thread sieves [low, high[
template <typename T>
T B_thread(T x,
           int64_t y,
           int64_t low,
           int64_t high)
{
  ASSERT(low > 0);
  ASSERT(low < high);
  int64_t sqrtx = isqrt(x);
  int64_t start = max(y, min(x / high, sqrtx));
  int64_t stop = min(x / low, sqrtx);
  primesieve::iterator it1(stop + 1, start);
  int64_t prime = it1.prev_prime();

  if (prime <= start)
    return 0;

  // The first iteration requires computing pi(x / prime)
  // using the prime counting function.
  int threads = 1;
  uint64_t xp = (uint64_t)(x / prime);
  int64_t pi_xp = pi_noprint(xp, threads);
  T sum = pi_xp;
  prime = it1.prev_prime();

  // All other iterations compute pi(x / prime)
  // using a prime sieve.
  primesieve::iterator it2(xp, high);
  it2.generate_next_primes();

  // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i])
  for (; prime > start; prime = it1.prev_prime())
  {
    xp = (uint64_t)(x / prime);

    for (; it2.primes_[it2.size_ - 1] <= xp; it2.generate_next_primes())
      pi_xp += it2.size_ - it2.i_;
    for (; it2.primes_[it2.i_] <= xp; it2.i_++)
      pi_xp += 1;

    sum += pi_xp;
  }

  return sum;
}

/// \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
/// Run time: O(n log log n), with n = x / y
/// Memory usage: O(n^(1/2))
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
  int64_t xy = (int64_t)(x / max(y, 1));
  LoadBalancerP2 loadBalancer(x, xy, threads, is_print);
  threads = loadBalancer.get_threads();

  // for (low = sqrt(x); low < x / y; low += dist)
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

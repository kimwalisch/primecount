///
/// @file  S2_hard.cpp
/// @brief Calculate the contribution of the hard special leaves using
///        a prime sieve. This is a multi-threaded implementation
///        which uses compression (PiTable & FactorTable) to reduce
///        the memory usage by about 10x.
///
///        Usually the computation of the hard special leaves
///        requires a binary indexed tree a.k.a. Fenwick tree to count
///        the number of unsieved elements in O(log n) time. But it
///        is actually much faster to simply count the number of
///        unsieved elements directly from the sieve array using the
///        POPCNT instruction. Hence this implementation does not use
///        a binary indexed tree.
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
#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <Sieve.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <generate_phi.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <LoadBalancer.hpp>
#include <min.hpp>
#include <print.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Compute the contribution of the hard special leaves using a
/// segmented sieve. Each thread processes the interval
/// [low, low + segments * segment_size[.
///
/// Note that in the Deleglise-Rivat paper it is suggested to use a
/// segment size of y. In practice however this uses too much memory
/// especially when using multi-threading. Hence we are using a
/// segment size of sqrt(z) as suggested in Xavier Gourdon's paper.
/// In primecount's implementation a segment size of sqrt(z) seems
/// ideal since slightly increasing the segment size decreases
/// performance because of cache misses and slightly decreasing the
/// segment size also decreases performance.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_thread(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int64_t low,
                 int64_t segments,
                 int64_t segment_size,
                 FactorTable& factor,
                 PiTable& pi,
                 Primes& primes,
                 Runtime& runtime)
{
  int64_t low1 = max(low, 1);
  int64_t limit = min(low + segments * segment_size, z + 1);
  int64_t max_b = pi[min(isqrt(x / low1), isqrt(z), y)];
  int64_t pi_sqrty = pi[isqrt(y)];
  T s2_hard = 0;

  if (c > max_b)
    return s2_hard;

  runtime.init_start();
  Sieve sieve(low, segment_size, max_b);
  auto phi = generate_phi(low, max_b, primes, pi);
  runtime.init_stop();

  // Segmented sieve of Eratosthenes
  for (; low < limit; low += segment_size)
  {
    // current segment [low, high[
    int64_t high = min(low + segment_size, limit);
    low1 = max(low, 1);

    // pre-sieve multiples of first c primes
    sieve.pre_sieve(c, low, high);

    int64_t count_low_high = sieve.count((high - 1) - low);
    int64_t b = c + 1;

    // For c + 1 <= b <= pi_sqrty
    // Find all special leaves: n = primes[b] * m
    // which satisfy: mu[m] != 0 && primes[b] < lpf[m] && low <= (x / n) < high
    for (int64_t end = min(pi_sqrty, max_b); b <= end; b++)
    {
      int64_t prime = primes[b];
      T x2 = x / prime;
      int64_t x2_div_high = min(fast_div(x2, high), y);
      int64_t min_m = max(x2_div_high, y / prime);
      int64_t max_m = min(fast_div(x2, low1), y);
      int64_t count = 0;
      int64_t start = 0;

      if (prime >= max_m)
        goto next_segment;

      factor.to_index(&min_m);
      factor.to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        // mu(m) != 0 && prime < lpf(m)
        if (prime < factor.mu_lpf(m))
        {
          int64_t fm = factor.get_number(m);
          int64_t xn = fast_div64(x2, fm);
          int64_t stop = xn - low;
          count += sieve.count(start, stop, low, high, count, count_low_high);
          start = stop + 1;
          int64_t phi_xn = phi[b] + count;
          int64_t mu_m = factor.mu(m);
          s2_hard -= mu_m * phi_xn;
        }
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    // For pi_sqrty < b <= pi_sqrtz
    // Find all hard special leaves: n = primes[b] * primes[l]
    // which satisfy: low <= (x / n) < high
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T x2 = x / prime;
      int64_t x2_div_low = min(fast_div(x2, low1), y);
      int64_t x2_div_high = min(fast_div(x2, high), y);
      int64_t l = pi[min(x2_div_low, z / prime)];
      int64_t min_hard = max(x2_div_high, prime);
      int64_t count = 0;
      int64_t start = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_hard; l--)
      {
        int64_t xn = fast_div64(x2, primes[l]);
        int64_t stop = xn - low;
        count += sieve.count(start, stop, low, high, count, count_low_high);
        start = stop + 1;
        int64_t phi_xn = phi[b] + count;
        s2_hard += phi_xn;
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    next_segment:;
  }

  return s2_hard;
}

/// Calculate the contribution of the hard special leaves.
///
/// This is a parallel S2_hard(x, y) implementation with advanced load
/// balancing. As most special leaves tend to be in the first segments
/// we start off with a tiny segment size and one segment per thread.
/// After each iteration we dynamically increase the segment size (until
/// it reaches some limit) or the number of segments.
///
/// S2_hard(x, y) has been parallelized using an idea devised by Xavier
/// Gourdon. The idea is to make the individual threads completely
/// independent from each other so that no thread depends on values
/// calculated by another thread. The benefit of this approach is that
/// the algorithm will scale well up to a very large number of CPU
/// cores. In order to make the threads independent from each other
/// each thread needs to precompute a lookup table of phi(x, a) values
/// (this is done in S2_hard_thread(x, y)) every time the thread starts
/// a new computation.
///
template <typename T, typename FactorTable, typename Primes>
T S2_hard_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 T s2_hard_approx,
                 Primes& primes,
                 FactorTable& factor,
                 int threads)
{
  threads = ideal_num_threads(threads, z);

  LoadBalancer loadBalancer(x, y, z, s2_hard_approx);
  int64_t max_prime = min(y, z / isqrt(y));
  PiTable pi(max_prime);

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    int64_t low = 0;
    int64_t segments = 0;
    int64_t segment_size = 0;
    T s2_hard = 0;
    Runtime runtime;

    while (loadBalancer.get_work(&low, &segments, &segment_size, s2_hard, runtime))
    {
      runtime.start();
      // Unsigned integer division is usually slightly
      // faster than signed integer division
      using UT = typename make_unsigned<T>::type;
      s2_hard = S2_hard_thread((UT) x, y, z, c, low, segments, segment_size, factor, pi, primes, runtime);
      runtime.stop();
    }
  }

  T s2_hard = (T) loadBalancer.get_result();

  return s2_hard;
}

} // namespace

namespace primecount {

int64_t S2_hard(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int64_t s2_hard_approx,
                int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return S2_hard_mpi(x, y, z, c, s2_hard_approx, threads);
#endif

  print("");
  print("=== S2_hard(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  double time = get_time();
  FactorTable<uint16_t> factor(y, threads);
  int64_t max_prime = min(y, z / isqrt(y));
  auto primes = generate_primes<int32_t>(max_prime);
  int64_t s2_hard = S2_hard_OpenMP(x, y, z, c, s2_hard_approx, primes, factor, threads);

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

#ifdef HAVE_INT128_T

int128_t S2_hard(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int128_t s2_hard_approx,
                 int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return S2_hard_mpi(x, y, z, c, s2_hard_approx, threads);
#endif

  print("");
  print("=== S2_hard(x, y) ===");
  print("Computation of the hard special leaves");
  print(x, y, c, threads);

  double time = get_time();
  int128_t s2_hard;

  // uses less memory
  if (y <= FactorTable<uint16_t>::max())
  {
    FactorTable<uint16_t> factor(y, threads);
    int64_t max_prime = min(y, z / isqrt(y));
    auto primes = generate_primes<uint32_t>(max_prime);
    s2_hard = S2_hard_OpenMP(x, y, z, c, s2_hard_approx, primes, factor, threads);
  }
  else
  {
    FactorTable<uint32_t> factor(y, threads);
    int64_t max_prime = min(y, z / isqrt(y));
    auto primes = generate_primes<int64_t>(max_prime);
    s2_hard = S2_hard_OpenMP(x, y, z, c, s2_hard_approx, primes, factor, threads);
  }

  print("S2_hard", s2_hard, time);
  return s2_hard;
}

#endif

} // namespace

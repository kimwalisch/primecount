///
/// @file  D4.cpp
/// @brief This is a highly optimized implementation of the D(x, y)
///        formula in Xavier Gourdon's prime counting algorithm. The D
///        formula is very similar to the formula of the hard special
///        leaves in the Deleglise-Rivat algorithm. Hence this
///        implementation is basically identical to S2_hard.cpp except
///        that the bounds have been changed slightly.
///
///        This implementation uses multi-threading with advanced load
///        balancing, it scales well up to a large number of CPU cores
///        because the compute threads are completely independent from
///        each other. This implementation also uses the highly
///        optimized Sieve class and the DFactorTable class which is a
///        compressed lookup table of moebius function values,
///        least prime factors and max prime factors.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <PiTable.hpp>
#include <Sieve.hpp>
#include <LoadBalancer.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <generate_phi.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <print.hpp>
#include <S2.hpp>

#include "DFactorTable.hpp"

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Compute the contribution of the hard special leaves using a
/// segmented sieve. Each thread processes the interval
/// [low, low + segments * segment_size[.
///
template <typename T, typename DFactorTable, typename Primes>
T D_thread(T x,
           int64_t x_star,
           int64_t xz,
           int64_t y,
           int64_t z,
           int64_t k,
           int64_t low,
           int64_t segments,
           int64_t segment_size,
           DFactorTable& factor,
           PiTable& pi,
           Primes& primes,
           Runtime& runtime)
{
  int64_t low1 = max(low, 1);
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t limit = min(low + segments * segment_size, xz + 1);
  int64_t max_b_prime = min3(isqrt(x / low1), isqrt(limit), x_star);
  int64_t max_b = pi[max_b_prime];
  T sum = 0;

  if (k > max_b)
    return sum;

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

    // pre-sieve multiples of first k primes
    sieve.pre_sieve(k, low, high);

    int64_t count_low_high = sieve.count((high - 1) - low);
    int64_t b = k + 1;

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves: n = primes[b] * m
    // In the interval: low <= (x / n) < high
    // Which satisfy:  mu[m] != 0 && lpf[m] > primes[b] && mpf[m] <= y
    for (int64_t end = min(pi_sqrtz, max_b); b <= end; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_div_low = min(fast_div(xp, low1), z);
      int64_t xp_div_high = min(fast_div(xp, high), z);
      int64_t min_m = max(xp_div_high, z / prime);
      int64_t max_m = min(x / ipow<T>(prime, 3), xp_div_low);

      int64_t count = 0;
      int64_t start = 0;

      if (prime >= max_m)
        goto next_segment;

      factor.to_index(&min_m);
      factor.to_index(&max_m);

      for (int64_t m = max_m; m > min_m; m--)
      {
        // mu[m] != 0 && 
        // lpf[m] > prime &&
        // mpf[m] <= y
        if (prime < factor.is_leaf(m))
        {
          int64_t fm = factor.get_number(m);
          int64_t xpm = fast_div64(xp, fm);
          int64_t stop = xpm - low;
          count += sieve.count(start, stop, low, high, count, count_low_high);
          start = stop + 1;
          int64_t phi_xpm = phi[b] + count;
          int64_t mu_m = factor.mu(m);
          sum -= mu_m * phi_xpm;
        }
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    // For pi_sqrtz < b <= pi_x_star
    // Find all special leaves: n = primes[b] * prime2
    // which satisfy: low <= (x / n) < high && prime2 <= y
    for (; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t xp_div_low = min(fast_div(xp, low1), y);
      int64_t xp_div_high = min(fast_div(xp, high), y);
      int64_t min_m = max(xp_div_high, prime);
      int64_t max_m = min(x / ipow<T>(prime, 3), xp_div_low);

      int64_t l = pi[max_m];
      int64_t count = 0;
      int64_t start = 0;

      if (prime >= primes[l])
        goto next_segment;

      for (; primes[l] > min_m; l--)
      {
        int64_t xpq = fast_div64(xp, primes[l]);
        int64_t stop = xpq - low;
        count += sieve.count(start, stop, low, high, count, count_low_high);
        start = stop + 1;
        int64_t phi_xpq = phi[b] + count;
        sum += phi_xpq;
      }

      phi[b] += count_low_high;
      count_low_high -= sieve.cross_off(b, prime);
    }

    next_segment:;
  }

  return sum;
}

/// Calculate the contribution of the hard special leaves.
///
/// This is a parallel D(x, y) implementation with advanced load
/// balancing. As most special leaves tend to be in the first segments
/// we start off with a tiny segment size and one segment per thread.
/// After each iteration we dynamically increase the segment size (until
/// it reaches some limit) or the number of segments.
///
/// D(x, y) has been parallelized using an idea devised by Xavier
/// Gourdon. The idea is to make the individual threads completely
/// independent from each other so that no thread depends on values
/// calculated by another thread. The benefit of this approach is that
/// the algorithm will scale well up to a very large number of CPU
/// cores. In order to make the threads independent from each other
/// each thread needs to precompute a lookup table of phi(x, a) values
/// (this is done in D_thread(x, y)) every time the thread starts
/// a new computation.
///
template <typename T, typename DFactorTable, typename Primes>
T D_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t k,
           T d_approx,
           Primes& primes,
           DFactorTable& factor,
           int threads)
{
  int64_t xz = x / z;
  int64_t x_star = get_x_star_gourdon(x, y);
  threads = ideal_num_threads(threads, xz);

  PiTable pi(y);
  LoadBalancer loadBalancer(x, xz, d_approx);

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    int64_t low = 0;
    int64_t segments = 0;
    int64_t segment_size = 0;
    T sum = 0;
    Runtime runtime;

    while (loadBalancer.get_work(&low, &segments, &segment_size, sum, runtime))
    {
      runtime.start();
      // Unsigned integer division is usually slightly
      // faster than signed integer division
      using UT = typename make_unsigned<T>::type;
      sum = D_thread((UT) x, x_star, xz, y, z, k, low, segments, segment_size, factor, pi, primes, runtime);
      runtime.stop();
    }
  }

  T sum = (T) loadBalancer.get_sum();

  return sum;
}

} // namespace

namespace primecount {

int64_t D(int64_t x,
          int64_t y,
          int64_t z,
          int64_t k,
          int64_t d_approx,
          int threads)
{
  print("");
  print("=== D(x, y) ===");
  print_gourdon(x, y, z, k, threads);

  double time = get_time();
  DFactorTable<uint16_t> factor(y, z, threads);
  auto primes = generate_primes<int32_t>(y);
  int64_t sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads);

  print("D", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t D(int128_t x,
           int64_t y,
           int64_t z,
           int64_t k,
           int128_t d_approx,
           int threads)
{
  print("");
  print("=== D(x, y) ===");
  print_gourdon(x, y, z, k, threads);

  double time = get_time();
  int128_t sum;

  // uses less memory
  if (z <= DFactorTable<uint16_t>::max())
  {
    DFactorTable<uint16_t> factor(y, z, threads);
    auto primes = generate_primes<uint32_t>(y);
    sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads);
  }
  else
  {
    DFactorTable<uint32_t> factor(y, z, threads);

    // uses less memory
    if (y <= numeric_limits<uint32_t>::max())
    {
      auto primes = generate_primes<uint32_t>(y);
      sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads);
    }
    else
    {
      auto primes = generate_primes<int64_t>(y);
      sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads); 
    }
  }

  print("D", sum, time);
  return sum;
}

#endif

} // namespace

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
  T sum = 0;
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t low1 = max(low, 1);
  int64_t limit = min(low + segments * segment_size, xz + 1);
  int64_t max_b = pi[min3(isqrt(x / low1), isqrt(limit), x_star)];
  int64_t min_b = pi[min(xz / limit, x_star)];
  min_b = max(k, min_b) + 1;

  if (min_b > max_b)
    return 0;

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

    // For i < min_b there are no special leaves:
    // low <= x / (primes[i] * m) < high
    sieve.pre_sieve(primes, min_b - 1, low, high);
    int64_t count_low_high = sieve.count((high - 1) - low);
    int64_t b = min_b;

    // For k + 1 <= b <= pi_sqrtz
    // Find all special leaves in the current segment that are
    // composed of a prime and a square free number:
    // low <= x / (primes[b] * m) < high
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
      count_low_high -= sieve.cross_off_count(prime, b);
    }

    // For pi_sqrtz < b <= pi_x_star
    // Find all special leaves in the current segment
    // that are composed of 2 primes:
    // low <= x / (primes[b] * primes[l]) < high
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
      count_low_high -= sieve.cross_off_count(prime, b);
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
           int threads,
           double& time)
{
  int64_t xz = x / z;
  LoadBalancer loadBalancer(x, y, z, k, xz, d_approx);
  int resume_threads = loadBalancer.get_resume_threads();
  int64_t x_star = get_x_star_gourdon(x, y);
  threads = ideal_num_threads(threads, xz);
  PiTable pi(y);

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    // 1st resume computations from backup file
    for (int j = i; j < resume_threads; j += threads)
    {
      int64_t low = 0;
      int64_t segments = 0;
      int64_t segment_size = 0;
      T sum = 0;
      Runtime runtime;

      if (loadBalancer.resume(j, low, segments, segment_size))
      {
        runtime.start();
        // Unsigned integer division is usually slightly
        // faster than signed integer division
        using UT = typename make_unsigned<T>::type;
        sum = D_thread((UT) x, x_star, xz, y, z, k, low, segments, segment_size, factor, pi, primes, runtime);
        loadBalancer.update_result(j, sum);
        runtime.stop();
      }
    }

    int64_t low = 0;
    int64_t segments = 0;
    int64_t segment_size = 0;
    T sum = 0;
    Runtime runtime;

    // 2nd get new work from loadBalancer
    while (loadBalancer.get_work(i, &low, &segments, &segment_size, sum, runtime))
    {
      runtime.start();
      // Unsigned integer division is usually slightly
      // faster than signed integer division
      using UT = typename make_unsigned<T>::type;
      sum = D_thread((UT) x, x_star, xz, y, z, k, low, segments, segment_size, factor, pi, primes, runtime);
      runtime.stop();
    }
  }

  loadBalancer.finish_backup();
  time = loadBalancer.get_wtime();
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
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int64_t xz = x / z;
  LoadBalancer loadBalancer(x, y, z, k, xz, d_approx);
  maxint_t sum = 0;

  if (!loadBalancer.resume(sum, time))
  {
    DFactorTable<uint16_t> factor(y, z, threads);
    auto primes = generate_primes<int32_t>(y);
    sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, time);
  }

  print("D", sum, time);
  return (int64_t) sum;
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
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int64_t xz = x / z;
  LoadBalancer loadBalancer(x, y, z, k, xz, d_approx);
  maxint_t sum = 0;

  if (!loadBalancer.resume(sum, time))
  {
    // uses less memory
    if (z <= DFactorTable<uint16_t>::max())
    {
      DFactorTable<uint16_t> factor(y, z, threads);
      auto primes = generate_primes<uint32_t>(y);
      sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, time);
    }
    else
    {
      DFactorTable<uint32_t> factor(y, z, threads);

      // uses less memory
      if (y <= numeric_limits<uint32_t>::max())
      {
        auto primes = generate_primes<uint32_t>(y);
        sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, time);
      }
      else
      {
        auto primes = generate_primes<int64_t>(y);
        sum = D_OpenMP(x, y, z, k, d_approx, primes, factor, threads, time);
      }
    }
  }

  print("D", sum, time);
  return sum;
}

#endif

} // namespace

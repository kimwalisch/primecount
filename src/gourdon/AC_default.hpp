///
/// @file  AC_default.hpp
/// @brief Implementation of the A + C formulas in Xavier Gourdon's
///        prime counting algorithm. In this implementation the memory
///        usage of the pi[x] lookup table has been reduced from
///        O(x^(1/2)) to O(x^(1/4)) by using a segmented pi[x] lookup
///        table. In each segment we process the leaves that satisfy:
///        low <= x / (prime * m) < high.
///
///        The A & C formulas roughly correspond to the easy special
///        leaves in the Deleglise-Rivat algorithm. Since both
///        formulas use a very similar segmented algorithm that goes
///        up to x^(1/2) it makes sense to merge the A & C formulas
///        hence reducing the runtime complexity by a factor of
///        O(x^(1/2) * ln ln x^(1/2)) and avoiding initializing some
///        data structures twice. Merging the A & C formulas also
///        improves scaling on systems with many CPU cores.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.pdf
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef AC_DEFAULT_HPP
#define AC_DEFAULT_HPP

#include "LoadBalancerAC.hpp"
#include "SegmentedPiTable.hpp"

#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>
#include <PiTable.hpp>

#include <stdint.h>

namespace {

using namespace primecount;

template <typename T,
          int MULTIPLIER,
          typename XP,
          typename Primes>
ALWAYS_INLINE T sum_pi(XP xp,
                       uint64_t i,
                       uint64_t last,
                       uint64_t b,
                       const Primes& primes,
                       const SegmentedPiTable& segmentedPi)
{
  if (i > last)
    return 0;

  T sum = 0;
  uint64_t size = last - i + 1;

  // Unroll loop to increase instruction level parallelism
  for (; i + 3 <= last; i += 4)
  {
    uint64_t xpq0 = fast_div64(xp, primes[i]);
    uint64_t xpq1 = fast_div64(xp, primes[i+1]);
    uint64_t xpq2 = fast_div64(xp, primes[i+2]);
    uint64_t xpq3 = fast_div64(xp, primes[i+3]);

    sum += segmentedPi[xpq0] +
           segmentedPi[xpq1] +
           segmentedPi[xpq2] +
           segmentedPi[xpq3];
  }

  NO_UNROLL_LOOP
  for (; i <= last; i++)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq];
  }

  return sum * MULTIPLIER + size * 2 - size * T(b);
}

/// Compute the A formula.
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename XP,
          typename Primes>
T A(T xlow,
    T xhigh,
    XP xp,
    uint64_t y,
    uint64_t b,
    const Primes& primes,
    const PiTable& pi,
    const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t prime = primes[b];
  uint64_t sqrt_xp = isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
  uint64_t i = pi[max(prime, min_2nd_prime)] + 1;
  uint64_t max_i1 = pi[min(xp / y, max_2nd_prime)];
  uint64_t max_i2 = pi[max_2nd_prime];

  // for (; i <= max_i1; i++)
  //   sum += segmentedPi[xp / primes[i]];
  sum += sum_pi<T, 1>(xp, i, max_i1, 2, primes, segmentedPi);
  i = max(i, max_i1 + 1);

  // for (; i <= max_i2; i++)
  //   sum += segmentedPi[xp / primes[i]] * 2;
  sum += sum_pi<T, 2>(xp, i, max_i2, 2, primes, segmentedPi);

  return sum;
}

/// Compute the 1st part of the C formula.
/// pi[(x/z)^(1/3)] < b <= pi[sqrt(z)]
/// x / (primes[b] * m) <= z
/// low <= x / (primes[b] * m) < high
///
/// m is either a prime <= y or a product of 2 distinct primes.
/// In both cases m is coprime to the first b primes, m <= z,
/// and its largest prime factor is <= y.
/// Since each prime factor of m is > (x / z)^(1/3) and z < sqrt(x),
/// m cannot contain more than 2 prime factors.
///
template <typename T,
          typename XP,
          typename Primes>
T C1(T xlow,
     T xhigh,
     XP xp,
     uint64_t b,
     uint64_t y,
     uint64_t z,
     const Primes& primes,
     const PiTable& pi,
     const SegmentedPiTable& segmentedPi)
{
  T sum = 0;
  uint64_t prime = primes[b];
  uint64_t max_m = min(xlow / prime, z);
  uint64_t x_div_prime3 = fast_div64(xp, prime * prime);
  uint64_t min_m = fast_div64(xhigh, prime);
  min_m = max3(min_m, x_div_prime3, z / prime);

  if (min_m >= max_m)
    return 0;

  // m = primes[i]
  uint64_t max_prime = min(y, max_m);

  if (min_m < max_prime)
  {
    uint64_t min_i = pi[min_m] + 1;
    uint64_t max_i = pi[max_prime];

    // for (i = min_i; i <= max_i; i++)
    //   sum -= segmentedPi[xp / primes[i]] - b + 2;
    sum -= sum_pi<T, 1>(xp, min_i, max_i, b, primes, segmentedPi);
  }

  // m = primes[i] * primes[j]
  uint64_t max_q = min(y, isqrt(max_m));
  uint64_t min_q = max(prime, min_m / y);

  if (min_q < max_q)
  {
    uint64_t min_q_i = pi[min_q] + 1;
    uint64_t max_q_i = pi[max_q];

    for (uint64_t i = min_q_i; i <= max_q_i; i++)
    {
      uint64_t q = primes[i];
      uint64_t min_r = max(q, min_m / q);
      uint64_t max_r = min(y, max_m / q);

      if (min_r >= max_r)
        continue;

      uint64_t min_j = pi[min_r] + 1;
      uint64_t max_j = pi[max_r];
      XP xpq = fast_div(xp, q);

      // for (j = min_j; j <= max_j; j++)
      //   sum += segmentedPi[xpq / primes[j]] - b + 2;
      sum += sum_pi<T, 1>(xpq, min_j, max_j, b, primes, segmentedPi);
    }
  }

  return sum;
}

/// Compute the 2nd part of the C formula.
/// Computes the clustered and sparse easy leaves of the C
/// formula for which the second factor is necessarily prime.
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename XP,
          typename Primes>
T C2(T xlow,
     T xhigh,
     XP xp,
     uint64_t y,
     uint64_t b,
     uint64_t pi_y,
     uint64_t max_clustered_global,
     const Primes& primes,
     const PiTable& pi,
     const SegmentedPiTable& segmentedPi)
{
  uint64_t prime = primes[b];
  uint64_t max_m = min3(xlow / prime, xp / prime, y);
  uint64_t x_div_prime3 = fast_div64(xp, prime * prime);
  uint64_t xhigh_div_prime = fast_div64(xhigh, prime);
  uint64_t min_m = max3(xhigh_div_prime, x_div_prime3, prime);

  if (min_m >= max_m)
    return 0;

  uint64_t pi_min_m = pi[min_m];
  uint64_t sqrt_xp = isqrt(xp);
  uint64_t min_clustered = in_between(min_m, sqrt_xp, max_m);
  uint64_t pi_min_clustered = pi[min_clustered];
  uint64_t min_clustered_global = max3(x_div_prime3, sqrt_xp, prime);
  min_clustered_global = min(min_clustered_global, y);
  uint64_t pi_conj_lo = pi_min_m;
  uint64_t pi_conj_hi = pi_min_m;
  uint64_t i = pi_min_m + 1;

  T sum = 0;

  // Compute the boundary correction once
  if (max_m >= max_clustered_global &&
      pi_y > pi_min_clustered)
  {
    uint64_t q_lo = fast_div64(xp, max_clustered_global);
    uint64_t q_hi = fast_div64(xp, min_clustered_global + 1);
    uint64_t pi_min_clustered_global = pi[min_clustered_global];
    sum += T(pi[q_lo]) * pi_y - T(pi[q_hi]) * pi_min_clustered_global;
    sum -= T(b - 2) * (pi_y - pi_min_clustered_global);
  }

  // Reflected range: ]pi_conj_lo, pi_conj_hi]
  if (min_clustered_global < max_clustered_global &&
      segmentedPi.low() < max_clustered_global &&
      segmentedPi.high() > min_clustered_global + 1)
  {
    uint64_t q_lo = fast_div64(xp, max_clustered_global);
    uint64_t q_hi = fast_div64(xp, min_clustered_global + 1);
    pi_conj_lo = max(pi[q_lo], pi_min_m);
    pi_conj_hi = min(pi[q_hi], pi_min_clustered);
    pi_conj_hi = max(pi_conj_hi, pi_conj_lo);
  }

  // Sparse leaves below the reflected range.
  // for (; i <= pi_conj_lo; i++)
  //   sum += segmentedPi[xp / primes[i]] - b + 2;
  sum += sum_pi<T, 1>(xp, i, pi_conj_lo, b, primes, segmentedPi);
  i = pi_conj_lo + 1;

  // Reflected leaves are counted once as a
  // sparse leaf and once as a conjugate.
  // for (; i <= pi_conj_hi; i++)
  //   sum += segmentedPi[xp / primes[i]] * 2 - b + 2;
  sum += sum_pi<T, 2>(xp, i, pi_conj_hi, b, primes, segmentedPi);
  i = pi_conj_hi + 1;

  // Sparse leaves above the reflected range.
  // for (; i <= pi_min_clustered; i++)
  //   sum += segmentedPi[xp / primes[i]] - b + 2;
  sum += sum_pi<T, 1>(xp, i, pi_min_clustered, b, primes, segmentedPi);

  return sum;
}

/// Compute A + C
template <typename T,
          typename Primes>
T AC_OpenMP_default(T x,
                    int64_t y,
                    int64_t z,
                    int64_t k,
                    int64_t x_star,
                    const PiTable& pi,
                    const Primes& primes,
                    int threads,
                    bool is_print)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t xy = x / y;
  int64_t xz = x / z;

  // These load balancing settings work well on my
  // dual-socket AMD EPYC 7642 server with 192 CPU cores.
  int64_t thread_threshold = 1000;
  int max_threads = (int) std::pow(xz, 1 / 3.7);
  threads = min(threads, max_threads);
  threads = ideal_num_threads(x13, threads, thread_threshold);
  INDETERMINATE LoadBalancerAC loadBalancer(sqrtx, y, threads, is_print);

  int64_t pi_y = pi[y];
  int64_t max_clustered_global = primes[pi_y];
  int64_t sqrtz = isqrt(z);
  int64_t pi_sqrtz = pi[sqrtz];
  int64_t pi_root3_xy = pi[iroot<3>(xy)];
  int64_t pi_root3_xz = pi[iroot<3>(xz)];

  // In order to reduce the thread creation & destruction
  // overhead we reuse the same threads throughout the
  // entire computation. The same threads are used for:
  //
  // 1) Computation of the C1 formula.
  // 2) Computation of the C2 formula.
  // 3) Computation of the A formula.
  //
  #pragma omp parallel num_threads(threads) reduction(+: sum)
  {
    // SegmentedPiTable is accessed very frequently.
    // In order to get good performance it is important that
    // SegmentedPiTable fits into the CPU's cache.
    // Hence we use a small segment_size of x^(1/4).
    SegmentedPiTable segmentedPi;
    ThreadDataAC thread;

    // for (low = 0; low < sqrt(x); low += segment_size)
    while (loadBalancer.get_work(thread))
    {
      int64_t low = thread.low;
      int64_t segment_size = thread.segment_size;
      int64_t limit = low + thread.segments * segment_size;
      limit = min(limit, sqrtx);

      for (; low < limit; low += segment_size)
      {
        // Current segment [low, high[
        int64_t high = low + segment_size;
        high = min(high, sqrtx);
        segmentedPi.init(low, high, limit);

        // We measure the thread computation time excluding the
        // first expensive initialization of the segmentedPi
        // lookup table. If the thread computation time is close
        // to 0 then we increase the number of segments in the
        // loadBalancer which should improve performance.
        if (low == thread.low)
          thread.secs = get_time();

        int64_t pi_sqrt_low = pi[isqrt(low)];
        T xlow = x / max(low, 1);
        T xhigh = x / high;

        if (low < z)
        {
          int64_t min_c1 = max(k, pi_root3_xz);
          int64_t min_c1_prime = min(xz / high, sqrtz);
          min_c1 = max3(min_c1, pi_sqrt_low, pi[min_c1_prime]) + 1;

          // C1 formula: pi[(x/z)^(1/3)] < b <= pi[sqrt(z)]
          for (int64_t b = min_c1; b <= pi_sqrtz; b++)
          {
            T xp = x / primes[b];

            if (xp <= pstd::numeric_limits<uint64_t>::max())
              sum -= C1(xlow, xhigh, uint64_t(xp), b, y, z, primes, pi, segmentedPi);
            else
              sum -= C1(xlow, xhigh, xp, b, y, z, primes, pi, segmentedPi);
          }
        }

        int64_t min_c2 = max3(k, pi_root3_xy, pi_sqrtz);
        int64_t min_c2_prime = min(xhigh / y, x_star);
        min_c2 = max3(min_c2, pi_sqrt_low, pi[min_c2_prime]) + 1;
        int64_t x_div_high2 = fast_div64(xhigh, high);
        int64_t min_a = min(x_div_high2, x13);
        min_a = pi[max(x_star, min_a)] + 1;

        // Upper bound of A & C2 formulas:
        // x / (p * q) >= low
        // p * next_prime(p) <= x / low
        // p <= sqrt(x / low)
        int64_t sqrt_xlow = isqrt(xlow);
        int64_t max_c2 = pi[min(sqrt_xlow, x_star)];
        T max_c2_prime = xlow / max(max_clustered_global, 1);
        int64_t max_c2_clustered = pi[min3(max_c2_prime, sqrt_xlow, x_star)];
        int64_t min_c2_sparse = pi[min(x_div_high2, x_star)] + 1;
        min_c2_sparse = max3(min_c2, min_c2_sparse, max_c2_clustered + 1);
        int64_t max_a = pi[min(sqrt_xlow, x13)];

        // C2 formula: pi[sqrt(z)] < b <= pi[x_star]
        for (int64_t b = min_c2; b <= max_c2_clustered; b++)
        {
          T xp = x / primes[b];

          if (xp <= pstd::numeric_limits<uint64_t>::max())
            sum += C2(xlow, xhigh, uint64_t(xp), y, b, pi_y, max_clustered_global, primes, pi, segmentedPi);
          else
            sum += C2(xlow, xhigh, xp, y, b, pi_y, max_clustered_global, primes, pi, segmentedPi);
        }

        // C2 formula: pi[sqrt(z)] < b <= pi[x_star]
        for (int64_t b = min_c2_sparse; b <= max_c2; b++)
        {
          T xp = x / primes[b];

          if (xp <= pstd::numeric_limits<uint64_t>::max())
            sum += C2(xlow, xhigh, uint64_t(xp), y, b, pi_y, max_clustered_global, primes, pi, segmentedPi);
          else
            sum += C2(xlow, xhigh, xp, y, b, pi_y, max_clustered_global, primes, pi, segmentedPi);
        }

        // A formula: pi[x_star] < b <= pi[x13]
        for (int64_t b = min_a; b <= max_a; b++)
        {
          T xp = x / primes[b];

          if (xp <= pstd::numeric_limits<uint64_t>::max())
            sum += A(xlow, xhigh, uint64_t(xp), y, b, primes, pi, segmentedPi);
          else
            sum += A(xlow, xhigh, xp, y, b, primes, pi, segmentedPi);
        }
      }
    }
  }

  return sum;
}

} // namespace

#endif

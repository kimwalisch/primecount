///
/// @file  AC_libdivide.cpp
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
///        This is an optimized version of AC(x, y) which uses
///        libdivide. libdivide allows to replace expensive integer
///        divsion instructions by a sequence of shift, add and
///        multiply instructions that will calculate the integer
///        division much faster.
///
///        In-depth description of this algorithm:
///        https://github.com/kimwalisch/primecount/blob/master/doc/Easy-Special-Leaves.md
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <LoadBalancerAC.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <gourdon.hpp>
#include <int128_t.hpp>
#include <libdivide.h>
#include <min.hpp>
#include <imath.hpp>
#include <pod_vector.hpp>
#include <print.hpp>
#include <RelaxedAtomic.hpp>

#include <stdint.h>

using std::numeric_limits;
using namespace primecount;

namespace {

/// Compute the A formula using libdivide.
/// 64-bit function: xp < 2^64
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename LibdividePrimes>
T A_64(T xlow,
       T xhigh,
       uint64_t xp,
       uint64_t y,
       uint64_t prime,
       const LibdividePrimes& primes,
       const PiTable& pi,
       const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t sqrt_xp = isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
  uint64_t i = pi[max(prime, min_2nd_prime)] + 1;
  uint64_t max_i1 = pi[min(xp / y, max_2nd_prime)];
  uint64_t max_i2 = pi[max_2nd_prime];

  // pq = primes[b] * primes[i]
  // x / pq >= y && low <= x / pq < high
  for (; i <= max_i1; i++)
  {
    uint64_t xpq = xp / primes[i];
    sum += segmentedPi[xpq];
  }

  // pq = primes[b] * primes[i]
  // x / pq < y && low <= x / pq < high
  for (; i <= max_i2; i++)
  {
    uint64_t xpq = xp / primes[i];
    sum += segmentedPi[xpq] * 2;
  }

  return sum;
}

/// Compute the A formula.
/// 128-bit function: xp >= 2^64
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename Primes>
T A_128(T xlow,
        T xhigh,
        T xp,
        uint64_t y,
        uint64_t prime,
        const Primes& primes,
        const PiTable& pi,
        const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t sqrt_xp = (uint64_t) isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
  uint64_t i = pi[max(prime, min_2nd_prime)] + 1;
  uint64_t max_i1 = pi[min(xp / y, max_2nd_prime)];
  uint64_t max_i2 = pi[max_2nd_prime];

  // pq = primes[b] * primes[i]
  // x / pq >= y && low <= x / pq < high
  for (; i <= max_i1; i++)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq];
  }

  // pq = primes[b] * primes[i]
  // x / pq < y && low <= x / pq < high
  for (; i <= max_i2; i++)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq] * 2;
  }

  return sum;
}

/// Compute the 1st part of the C formula.
/// pi[(x/z)^(1/3)] < b <= pi[sqrt(z)]
/// x / (primes[b] * m) <= z
///
/// m may be a prime <= y or a square free number <= z which is
/// coprime to the first b primes and whose largest prime factor <= y.
/// This algorithm recursively iterates over the square free numbers
/// coprime to the first b primes. This algorithm is described in
/// section 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March 2015.
///
template <int MU, 
          typename T, 
          typename Primes>
T C1(T xp,
     uint64_t b,
     uint64_t i,
     uint64_t pi_y,
     uint64_t m,
     uint64_t min_m,
     uint64_t max_m,
     const Primes& primes,
     const PiTable& pi)
{
  T sum = 0;

  for (i++; i <= pi_y; i++)
  {
    // Calculate next m
    T m128 = (T) m * primes[i];
    if (m128 > max_m)
      return sum;

    uint64_t m64 = (uint64_t) m128;

    if (m64 > min_m) {
      uint64_t xpm = fast_div64(xp, m64);
      T phi_xpm = pi[xpm] - b + 2;
      sum += phi_xpm * MU;
    }

    sum += C1<-MU>(xp, b, i, pi_y, m64, min_m, max_m, primes, pi);
  }

  return sum;
}

/// Compute the 2nd part of the C formula.
/// 64-bit function: xp < 2^64
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T, 
          typename LibdividePrimes>
T C2_64(T xlow,
        T xhigh,
        uint64_t xp,
        uint64_t y,
        uint64_t b,
        uint64_t prime,
        const LibdividePrimes& primes,
        const PiTable& pi,
        const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t max_m = min3(xlow / prime, xp / prime, y);
  T min_m128 = max3(xhigh / prime, xp / (prime * prime), prime);
  uint64_t min_m = min(min_m128, max_m);
  uint64_t i = pi[max_m];
  uint64_t pi_min_m = pi[min_m];
  uint64_t min_clustered = isqrt(xp);
  min_clustered = in_between(min_m, min_clustered, max_m);
  uint64_t pi_min_clustered = pi[min_clustered];

  // Find all clustered easy leaves where
  // successive leaves are identical.
  // pq = primes[b] * primes[i]
  // Which satisfy: low <= x / pq < high && q <= y && pq > z
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  while (i > pi_min_clustered)
  {
    uint64_t xpq = xp / primes[i];
    uint64_t pi_xpq = segmentedPi[xpq];
    uint64_t phi_xpq = pi_xpq - b + 2;
    uint64_t xpq2 = xp / primes[pi_xpq + 1];
    uint64_t imin = pi[max(xpq2, min_clustered)];
    sum += phi_xpq * (i - imin);
    i = imin;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // pq = primes[b] * primes[i]
  // Which satisfy: low <= x / pq < high && q <= y && pq > z
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  for (; i > pi_min_m; i--)
  {
    uint64_t xpq = xp / primes[i];
    sum += segmentedPi[xpq] - b + 2;
  }

  return sum;
}

/// Compute the 2nd part of the C formula.
/// 128-bit function: xp >= 2^64
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename Primes>
T C2_128(T xlow,
         T xhigh,
         T xp,
         uint64_t y,
         uint64_t b,
         const Primes& primes,
         const PiTable& pi,
         const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t prime = primes[b];
  uint64_t max_m = min3(xlow / prime, xp / prime, y);
  T min_m128 = max3(xhigh / prime, xp / (prime * prime), prime);
  uint64_t min_m = min(min_m128, max_m);
  uint64_t i = pi[max_m];
  uint64_t pi_min_m = pi[min_m];
  uint64_t min_clustered = (uint64_t) isqrt(xp);
  min_clustered = in_between(min_m, min_clustered, max_m);
  uint64_t pi_min_clustered = pi[min_clustered];

  // Find all clustered easy leaves where
  // successive leaves are identical.
  // pq = primes[b] * primes[i]
  // Which satisfy: low <= x / pq < high && q <= y && pq > z
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  while (i > pi_min_clustered)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    uint64_t pi_xpq = segmentedPi[xpq];
    uint64_t phi_xpq = pi_xpq - b + 2;
    uint64_t xpq2 = fast_div64(xp, primes[pi_xpq + 1]);
    uint64_t imin = pi[max(xpq2, min_clustered)];
    sum += phi_xpq * (i - imin);
    i = imin;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // pq = primes[b] * primes[i]
  // Which satisfy: low <= x / pq < high && q <= y && pq > z
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  for (; i > pi_min_m; i--)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq] - b + 2;
  }

  return sum;
}

/// Compute A + C
template <typename T,
          typename Primes>
T AC_OpenMP(T x,
            int64_t y,
            int64_t z,
            int64_t k,
            int64_t x_star,
            int64_t max_a_prime,
            const Primes& primes,
            int threads,
            bool is_print)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);
  LoadBalancerAC loadBalancer(sqrtx, y, threads, is_print);

  // Initialize libdivide vector from primes vector
  pod_vector<libdivide::branchfree_divider<uint64_t>> lprimes;
  lprimes.resize(primes.size());
  for (std::size_t i = 1; i < lprimes.size(); i++)
    lprimes[i] = primes[i];

  // PiTable's size = z because of the C1 formula.
  // PiTable is accessed much less frequently than
  // SegmentedPiTable, hence it is OK that PiTable's size
  // is fairly large and does not fit into the CPU's cache.
  PiTable pi(max(z, max_a_prime), threads);

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_root3_xy = pi[iroot<3>(x / y)];
  int64_t pi_root3_xz = pi[iroot<3>(x / z)];
  RelaxedAtomic<int64_t> min_c1(max(k, pi_root3_xz) + 1);

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
    int64_t low, high;

    // C1 formula: pi[(x/z)^(1/3)] < b <= pi[pi_sqrtz]
    for (int64_t b = min_c1++; b <= pi_sqrtz; b = min_c1++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t max_m = min(xp / prime, z);
      T min_m128 = max(xp / (prime * prime), z / prime);
      int64_t min_m = min(min_m128, max_m);

      sum -= C1<-1>(xp, b, b, pi_y, 1, min_m, max_m, primes, pi);
    }

    // for (low = 0; low < sqrt; low += segment_size)
    while (loadBalancer.get_work(low, high))
    {
      // Current segment [low, high[
      segmentedPi.init(low, high);
      T xlow = x / max(low, 1);
      T xhigh = x / high;

      int64_t min_c2 = max(k, pi_root3_xy);
      min_c2 = max(min_c2, pi_sqrtz);
      min_c2 = max(min_c2, pi[isqrt(low)]);
      min_c2 = max(min_c2, pi[min(xhigh / y, x_star)]);
      min_c2 += 1;

      int64_t min_a = min(xhigh / high, x13);
      min_a = pi[max(x_star, min_a)] + 1;

      // Upper bound of A & C2 formulas:
      // x / (p * q) >= low
      // p * next_prime(p) <= x / low
      // p <= sqrt(x / low)
      T sqrt_xlow = isqrt(xlow);
      int64_t max_c2 = pi[min(sqrt_xlow, x_star)];
      int64_t max_a = pi[min(sqrt_xlow, x13)];

      // C2 formula: pi[sqrt(z)] < b <= pi[x_star]
      for (int64_t b = min_c2; b <= max_c2; b++)
      {
        int64_t prime = primes[b];
        T xp = x / prime;

        if (xp <= numeric_limits<uint64_t>::max())
          sum += C2_64(xlow, xhigh, (uint64_t) xp, y, b, prime, lprimes, pi, segmentedPi);
        else
          sum += C2_128(xlow, xhigh, xp, y, b, primes, pi, segmentedPi);
      }

      // A formula: pi[x_star] < b <= pi[x13]
      for (int64_t b = min_a; b <= max_a; b++)
      {
        int64_t prime = primes[b];
        T xp = x / prime;

        if (xp <= numeric_limits<uint64_t>::max())
          sum += A_64(xlow, xhigh, (uint64_t) xp, y, prime, lprimes, pi, segmentedPi);
        else
          sum += A_128(xlow, xhigh, xp, y, prime, primes, pi, segmentedPi);
      }
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t AC(int64_t x,
           int64_t y,
           int64_t z,
           int64_t k,
           int threads,
           bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== AC(x, y) ===");
    print_gourdon_vars(x, y, z, k, threads);
  }

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  auto primes = generate_primes<uint32_t>(max_prime);

  int64_t sum = AC_OpenMP((uint64_t) x, y, z, k, x_star, max_a_prime, primes, threads, is_print);

  if (is_print)
    print("A + C", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t AC(int128_t x,
            int64_t y,
            int64_t z,
            int64_t k,
            int threads,
            bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== AC(x, y) ===");
    print_gourdon_vars(x, y, z, k, threads);
  }

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  int128_t sum;

  // uses less memory
  if (max_prime <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(max_prime);
    sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, max_a_prime, primes, threads, is_print);
  }
  else
  {
    auto primes = generate_primes<uint64_t>(max_prime);
    sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, max_a_prime, primes, threads, is_print);
  }

  if (is_print)
    print("A + C", sum, time);

  return sum;
}

#endif

} // namespace

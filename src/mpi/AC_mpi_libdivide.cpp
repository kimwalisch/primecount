///
/// @file  AC_mpi_libdivide.cpp
/// @brief Implementation of the A + C formulas (from Xavier Gourdon's
///        algorithm) that have been distributed using MPI (Message
///        Passing Interface) and multi-threaded using OpenMP.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <mpi_reduce_sum.hpp>
#include <fast_div.hpp>
#include <for_atomic.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <libdivide.h>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <StatusAC.hpp>

#include <stdint.h>
#include <atomic>
#include <vector>

using namespace std;
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

  // x / (p * q) >= y
  for (; i <= max_i1; i++)
  {
    uint64_t xpq = xp / primes[i];
    sum += segmentedPi[xpq];
  }

  // x / (p * q) < y
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

  // x / (p * q) >= y
  for (; i <= max_i1; i++)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq];
  }

  // x / (p * q) < y
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
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
  while (i > pi_min_clustered)
  {
    uint64_t xpq = xp / primes[i];
    uint64_t phi_xpq = segmentedPi[xpq] - b + 2;
    uint64_t xpq2 = xp / primes[b + phi_xpq - 1];
    uint64_t i2 = pi[max(xpq2, min_clustered)];
    sum += phi_xpq * (i - i2);
    i = i2;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
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
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
  while (i > pi_min_clustered)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    uint64_t phi_xpq = segmentedPi[xpq] - b + 2;
    uint64_t xpq2 = fast_div64(xp, primes[b + phi_xpq - 1]);
    uint64_t i2 = pi[max(xpq2, min_clustered)];
    sum += phi_xpq * (i - i2);
    i = i2;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
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
            int threads)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);
  StatusAC status;

  // Initialize libdivide vector using primes
  vector<libdivide::branchfree_divider<uint64_t>> lprimes(1);
  lprimes.insert(lprimes.end(), primes.begin() + 1, primes.end());

  // PiTable's size = z because of the C1 formula.
  // PiTable is accessed much less frequently than
  // SegmentedPiTable, hence it is OK that PiTable's size
  // is fairly large and does not fit into the CPU's cache. 
  PiTable pi(max(z, max_a_prime), threads);

  // SegmentedPiTable is accessed very frequently.
  // In order to get good performance it is important that
  // SegmentedPiTable fits into the CPU's cache.
  // Hence we use a small size of x^(1/3).
  SegmentedPiTable segmentedPi(sqrtx, x13, threads);

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_root3_xy = pi[iroot<3>(x / y)];
  int64_t pi_root3_xz = pi[iroot<3>(x / z)];
  int64_t min_c1 = max(k, pi_root3_xz) + 1;
  int proc_id = mpi_proc_id();
  int procs = mpi_num_procs();

  atomic<int64_t> atomic_a(-1);
  atomic<int64_t> atomic_c1(-1);
  atomic<int64_t> atomic_c2(-1);

  // In order to reduce the thread creation & destruction
  // overhead we reuse the same threads throughout the
  // entire computation. The same threads are used for:
  //
  // 1) Computation of the C1 formula.
  // 2) Initialization of the segmentedPi lookup table.
  // 3) Computation of the C2 formula.
  // 4) Computation of the A formula.
  //
  #pragma omp parallel num_threads(threads) reduction(+: sum)
  {
    status.print(segmentedPi);

    // C1 formula: pi[(x/z)^(1/3)] < b <= pi[pi_sqrtz]
    for_atomic_add(min_c1 + proc_id, b <= pi_sqrtz, procs, atomic_c1)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t max_m = min(xp / prime, z);
      T min_m128 = max(xp / (prime * prime), z / prime);
      int64_t min_m = min(min_m128, max_m);
      sum -= C1<-1>(xp, b, b, pi_y, 1, min_m, max_m, primes, pi);
    }

    // This computes A and the 2nd part of the C formula.
    // Find all special leaves of type:
    // x / (primes[b] * primes[i]) < x^(1/2)
    // where b is bounded by pi[z^(1/2)] < b <= pi[x^(1/3)].
    // Since we need to lookup PrimePi[n] values for n < x^(1/2)
    // we use a segmented PrimePi[n] table of size y
    // (y = O(x^(1/3) * log(x)^3)) to reduce the memory usage.
    while (segmentedPi.low() < sqrtx)
    {
      // Current segment [low, high[
      segmentedPi.init();
      status.print(segmentedPi);
      int64_t low = segmentedPi.low();
      int64_t high = segmentedPi.high();
      T xlow = x / max(low, 1);
      T xhigh = x / high;

      int64_t min_c2 = max(k, pi_root3_xy);
      min_c2 = max(min_c2, pi_sqrtz);
      min_c2 = max(min_c2, pi[isqrt(low)]);
      min_c2 = max(min_c2, pi[min(xhigh / y, x_star)]);
      min_c2 += 1;

      // Upper bound of A & C2 formulas:
      // x / (p * q) >= low
      // p * next_prime(p) <= x / low
      // p <= sqrt(x / low)
      T sqrt_xlow = isqrt(xlow);
      int64_t max_c2 = pi[min(sqrt_xlow, x_star)];
      int64_t max_b = pi[min(sqrt_xlow, x13)];

      // C2 formula: pi[sqrt(z)] < b <= pi[x_star]
      for_atomic_add(min_c2 + proc_id, b <= max_c2, procs, atomic_c2)
      {
        int64_t prime = primes[b];
        T xp = x / prime;

        if (xp <= numeric_limits<uint64_t>::max())
          sum += C2_64(xlow, xhigh, (uint64_t) xp, y, b, prime, lprimes, pi, segmentedPi);
        else
          sum += C2_128(xlow, xhigh, xp, y, b, primes, pi, segmentedPi);
      }

      // A formula: pi[x_star] < b <= pi[x13]
      for_atomic_add(pi_x_star + 1 + proc_id, b <= max_b, procs, atomic_a)
      {
        int64_t prime = primes[b];
        T xp = x / prime;

        if (xp <= numeric_limits<uint64_t>::max())
          sum += A_64(xlow, xhigh, (uint64_t) xp, y, prime, lprimes, pi, segmentedPi);
        else
          sum += A_128(xlow, xhigh, xp, y, prime, primes, pi, segmentedPi);
      }

      // Is this the last segment?
      if (high >= sqrtx)
        break;

      // Wait until all threads have finished
      // computing the current segment.
      #pragma omp barrier
      #pragma omp master
      {
        segmentedPi.next();
        atomic_a = -1;
        atomic_c2 = -1;
      }

      #pragma omp barrier
    }
  }

  sum = mpi_reduce_sum(sum);

  return sum;
}

} // namespace

namespace primecount {

int64_t AC_mpi(int64_t x,
               int64_t y,
               int64_t z,
               int64_t k,
               int threads)
{
  print("");
  print("=== AC_mpi(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  auto primes = generate_primes<uint32_t>(max_prime);

  int64_t sum = AC_OpenMP((uint64_t) x, y, z, k, x_star, max_a_prime, primes, threads);

  print("A + C", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t AC_mpi(int128_t x,
                int64_t y,
                int64_t z,
                int64_t k,
                int threads)
{
  print("");
  print("=== AC_mpi(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

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
    sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, max_a_prime, primes, threads);
  }
  else
  {
    auto primes = generate_primes<uint64_t>(max_prime);
    sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, max_a_prime, primes, threads);
  }

  print("A + C", sum, time);
  return sum;
}

#endif

} // namespace

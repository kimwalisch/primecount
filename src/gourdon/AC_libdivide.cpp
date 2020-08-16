///
/// @file  AC_libdivide.cpp
/// @brief Implementation of the A + C formulas in Xavier Gourdon's
///        prime counting algorithm. In this version the memory usage
///        has been reduced from O(x^(1/2)) to O(z) by segmenting
///        the pi[x] lookup table. In each segment we process the
///        leaves that satisfy: low <= x / (prime * m) < high.
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
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <libdivide.h>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <Status.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Compute the A formula.
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) <= x^(1/2)
///
template <typename T,
          typename LibdividePrimes>
T A(T xlow,
    T xhigh,
    uint64_t xp,
    uint64_t y,
    uint64_t b,
    uint64_t prime,
    const PiTable& pi,
    const LibdividePrimes& primes,
    const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t sqrt_xp = isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t i = pi[min_2nd_prime];
  i = max(i, b) + 1;
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
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

/// 128-bit function.
/// Compute the A formula.
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) <= x^(1/2)
///
template <typename T,
          typename Primes>
T A(T xlow,
    T xhigh,
    T xp,
    uint64_t y,
    uint64_t b,
    const PiTable& pi,
    const Primes& primes,
    const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t prime = primes[b];
  uint64_t sqrt_xp = (uint64_t) isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t i = pi[min_2nd_prime];
  i = max(i, b) + 1;
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
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
/// k < b <= pi[sqrt(z)]
/// x / (primes[b] * m) <= z
/// 
/// Recursively iterate over the square free numbers coprime
/// to the first b primes. This algorithm is described in
/// section 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March
/// 2015.
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
     const PiTable& pi,
     const Primes& primes)
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

      if (MU > 0)
        sum += pi[xpm] - b + 2;
      else
        sum -= pi[xpm] - b + 2;
    }

    sum += C1<-MU>(xp, b, i, pi_y, m64, min_m, max_m, pi, primes);
  }

  return sum;
}

/// Compute the 2nd part of the C formula.
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) <= x^(1/2)
///
template <typename T, 
          typename LibdividePrimes>
T C2(T xlow,
     T xhigh,
     uint64_t xp,
     uint64_t y,
     uint64_t b,
     uint64_t prime,
     const PiTable& pi,
     const LibdividePrimes& primes,
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
    uint64_t i2 = segmentedPi[xpq2];
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

/// 128-bit function.
/// Compute the 2nd part of the C formula.
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) <= x^(1/2)
///
template <typename T,
          typename Primes>
T C2(T xlow,
     T xhigh,
     T xp,
     uint64_t y,
     uint64_t b,
     const PiTable& pi,
     const Primes& primes,
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
    uint64_t i2 = segmentedPi[xpq2];
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
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);

  Status status(x);
  PiTable pi(max(z, max_a_prime));
  SegmentedPiTable segmentedPi(isqrt(x), z, threads);

  // Initialize libdivide vector using primes
  using libdivide_t = libdivide::branchfree_divider<uint64_t>;
  vector<libdivide_t> lprimes(1);
  lprimes.insert(lprimes.end(),
                 primes.begin() + 1,
                 primes.end());

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_x13 = pi[x13];
  int64_t pi_root3_xy = pi[iroot<3>(x / y)];
  int64_t pi_root3_xz = pi[iroot<3>(x / z)];
  int64_t min_b = max(k, pi_root3_xz) + 1;

  int64_t max_b = 0;
  int64_t low = 0;
  int64_t high = 0;
  T xlow = 0;
  T xhigh = 0;

  #pragma omp parallel num_threads(threads)
  {
    // This computes the 1st part of the C formula.
    // Find all special leaves of type:
    // x / (primes[b] * m) <= z.
    // m may be a prime <= y or a square free number <= z
    // who is coprime to the first b primes and whose
    // largest prime factor <= y.
    #pragma omp for schedule(dynamic) reduction(-: sum)
    for (int64_t b = min_b; b <= pi_sqrtz; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t max_m = min(xp / prime, z);
      T min_m128 = max(xp / (prime * prime), z / prime);
      int64_t min_m = min(min_m128, max_m);

      sum -= C1<-1>(xp, b, b, pi_y, 1, min_m, max_m, pi, primes);

      if (is_print())
        status.print(b, pi_x13);
    }

    // This computes A and the 2nd part of the C formula.
    // Find all special leaves of type:
    // x / (primes[b] * primes[i]) <= x^(1/2)
    // with z^(1/2) < primes[b] <= x^(1/3).
    // Since we need to lookup PrimePi[n] values for n <= x^(1/2)
    // we use a segmented PrimePi[n] table of size z (~O(x^1/3))
    // in order to reduce the memory usage.
    while (!segmentedPi.finished())
    {
      #pragma omp single
      {
        // Current segment [low, high[
        low = segmentedPi.low();
        high = segmentedPi.high();
        low = max(low, 1);
        xlow = x / low;
        xhigh = x / high;

        // Lower bounds of C2 formula
        min_b = max3(k, pi_sqrtz, pi_root3_xy);
        min_b = max(min_b, pi[isqrt(low)]);
        min_b = max(min_b, pi[min(xhigh / y, x_star)]);
        min_b += 1;

        // Upper bound of A & C2 formulas:
        // x / (p * q) >= low
        // p * next_prime(p) <= x / low
        // p <= sqrt(x / low)
        int64_t sqrt_xlow = min(isqrt(xlow), x13);
        max_b = pi[sqrt_xlow];
      }

      // C2 formula: pi[sqrt(z)] < b <= pi[x_star]
      // A  formula: pi[x_star] < b <= pi[x13]
      #pragma omp for schedule(dynamic) reduction(+: sum)
      for (int64_t b = min_b; b <= max_b; b++)
      {
        int64_t prime = primes[b];
        T xp = x / prime;

        if (b <= pi_x_star)
        {
          if (xp <= numeric_limits<uint64_t>::max())
            sum += C2(xlow, xhigh, (uint64_t) xp, y, b, prime, pi, lprimes, segmentedPi);
          else
            sum += C2(xlow, xhigh, xp, y, b, pi, primes, segmentedPi);
        }
        else
        {
          if (xp <= numeric_limits<uint64_t>::max())
            sum += A(xlow, xhigh, (uint64_t) xp, y, b, prime, pi, lprimes, segmentedPi);
          else
            sum += A(xlow, xhigh, xp, y, b, pi, primes, segmentedPi);
        }

        if (is_print())
          status.print(b, pi_x13);
      }

      segmentedPi.next();
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
           int threads)
{
  print("");
  print("=== AC(x, y) ===");
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

int128_t AC(int128_t x,
            int64_t y,
            int64_t z,
            int64_t k,
            int threads)
{
  print("");
  print("=== AC(x, y) ===");
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

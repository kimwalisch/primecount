///
/// @file  C4.cpp
/// @brief Implementation of the C(x, y) formula in Xavier Gourdon's
///        prime counting algorithm. In this version the memory usage
///        has been reduced from O(x^(1/2)) to O(z) by segmenting
///        the pi[x] lookup table. In each segment we process the
///        leaves that satisfy: low <= x / (prime1 * prime2) < high.
///
///        In this implementation the easy special leaves have been
///        split up into 2 distinct types. Below sqrt(z) the leaves
///        are composed of a prime and a square free number. But when
///        the prime factors are > sqrt(z) then all leaves are
///        composed of exactly 2 primes.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
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
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <S2Status.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Recursively iterate over the square free numbers coprime
/// to the first b primes. This algorithm is described in
/// section 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March
/// 2015.
///
template <int MU, typename T, typename Primes>
T C1(T xp,
     int64_t b,
     uint64_t i,
     int64_t m,
     int64_t min_m,
     int64_t max_m,
     Primes& primes,
     PiTable& pi)
{
  T sum = 0;

  for (i++; i < primes.size(); i++)
  {
    // Calculate next m
    T m128 = (T) m * primes[i];
    if (m128 > (T) max_m)
      return sum;

    int64_t m64 = (int64_t) m128;
    if (m64 > min_m) {
      int64_t xpm = fast_div64(xp, m64);
      sum += MU * (pi[xpm] - b + 2);
    }

    sum += C1<-MU>(xp, b, i, m64, min_m, max_m, primes, pi);
  }

  return sum;
}

template <typename T, typename Primes>
T C_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t k,
           Primes& primes,
           int threads)
{
  int64_t root3_xy = iroot<3>(x / y);
  int64_t root3_xz = iroot<3>(x / z);
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);

  S2Status status(x);
  PiTable pi(max(root3_xy, z));
  SegmentedPiTable segmentedPi(isqrt(x), z, threads);

  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_root3_xy = pi[root3_xy];
  int64_t pi_root3_xz = pi[root3_xz];
  int64_t min_b = max(k, pi_root3_xz) + 1;
  T sum = 0;

  // This computes the 1st part of the C formula.
  // Find all special leaves of type:
  // x / (primes[b] * m) <= z.
  // m may be a prime <= y or a square free number <= z
  // who is coprime to the first b primes and whose
  // largest prime factor <= y.
  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(-: sum)
  for (int64_t b = min_b; b <= pi_sqrtz; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t max_m = min(xp / prime, z);
    T min_m128 = max(x / ipow<T>(prime, 3), z / prime);
    int64_t min_m = min(min_m128, max_m);

    sum -= C1<-1>(xp, b, b, 1, min_m, max_m, primes, pi);

    if (is_print())
      status.print(b, pi_x_star);
  }

  // Find all special leaves of type:
  // z < x / (primes[b] * primes[i]) <= x^(1/2)
  // with z^(1/2) < primes[b] <= x_star.
  // Since we need to lookup PrimePi[n] values for n <= x^(1/2)
  // we use a segmented PrimePi[n] table of size z (~O(x^1/3))
  // in order to reduce the memory usage.
  for (; !segmentedPi.finished(); segmentedPi.next())
  {
    // Current segment [low, high[
    int64_t low = segmentedPi.low();
    int64_t high = segmentedPi.high();
    low = max(low, 1);
    T xlow = x / low;
    T xhigh = x / high;

    min_b = max3(k, pi_sqrtz, pi_root3_xy);
    min_b = max(min_b, pi[isqrt(low)]);
    min_b = max(min_b, pi[min(xhigh / y, x_star)]);
    min_b += 1;

    #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
    for (int64_t b = min_b; b <= pi_x_star; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t max_m = min3(xlow / prime, xp / prime, y);
      T min_m128 = max3(xhigh / prime, x / ipow<T>(prime, 3), prime);
      int64_t min_m = min(min_m128, max_m);

      int64_t i = pi[max_m];
      int64_t pi_min_m = pi[min_m];
      int64_t min_clustered = (int64_t) isqrt(xp);
      min_clustered = in_between(min_m, min_clustered, max_m);
      int64_t pi_min_clustered = pi[min_clustered];

      // Find all clustered easy leaves where
      // successive leaves are identical.
      // n = primes[b] * primes[i]
      // Which satisfy: n > z && primes[i] <= y
      while (i > pi_min_clustered)
      {
        int64_t xpq = fast_div64(xp, primes[i]);
        int64_t phi_xpq = segmentedPi[xpq] - b + 2;
        int64_t xpq2 = fast_div64(xp, primes[b + phi_xpq - 1]);
        int64_t i2 = segmentedPi[xpq2];
        sum += phi_xpq * (i - i2);
        i = i2;
      }

      // Find all sparse easy leaves where
      // successive leaves are different.
      // n = primes[b] * primes[i]
      // Which satisfy: n > z && primes[i] <= y
      for (; i > pi_min_m; i--)
      {
        int64_t xpq = fast_div64(xp, primes[i]);
        sum += segmentedPi[xpq] - b + 2;
      }

      if (is_print())
        status.print(b, pi_x_star);
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t C(int64_t x,
          int64_t y,
          int64_t z,
          int64_t k,
          int threads)
{
  print("");
  print("=== C(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  auto primes = generate_primes<int32_t>(y);
  int64_t sum = C_OpenMP((intfast64_t) x, y, z, k, primes, threads);

  print("C", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t C(int128_t x,
           int64_t y,
           int64_t z,
           int64_t k,
           int threads)
{
  print("");
  print("=== C(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int128_t sum;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    sum = C_OpenMP((intfast128_t) x, y, z, k, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    sum = C_OpenMP((intfast128_t) x, y, z, k, primes, threads);
  }

  print("C", sum, time);
  return sum;
}

#endif

} // namespace

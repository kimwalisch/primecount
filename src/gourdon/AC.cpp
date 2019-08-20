///
/// @file  AC.cpp
/// @brief Implementation of the A + C formulas in Xavier Gourdon's
///        prime counting algorithm. In this version the memory usage
///        has been reduced from O(x^(1/2)) to O(z) by segmenting
///        the pi[x] lookup table. In each segment we process the
///        leaves that satisfy: low <= x / (prime1 * prime2) < high.
///
///        This is an optimized version of AC(x, y) which uses
///        libdivide. libdivide allows to replace expensive integer
///        divsion instructions by a sequence of shift, add and
///        multiply instructions that will calculate the integer
///        division much faster.
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
T C(T xp,
    uint64_t i,
    uint64_t m,
    uint64_t min_m,
    uint64_t max_m,
    uint64_t b,
    Primes& primes,
    PiTable& pi)
{
  T sum = 0;

  for (i++; i < primes.size(); i++)
  {
    uint64_t next_m = m * primes[i];

    // next_m may be 128-bit
    if ((T) m * primes[i] > max_m)
      return sum;

    if (next_m > min_m) {
      uint64_t xpm = fast_div64(xp, next_m);
      sum += MU * (pi[xpm] - b + 2);
    }

    sum += C<-MU>(xp, i, next_m, min_m, max_m, b, primes, pi);
  }

  return sum;
}

template <typename T, typename Primes>
T AC_OpenMP(T x,
            int64_t y,
            int64_t z,
            int64_t k,
            int64_t x_star,
            int64_t max_a_prime,
            Primes& primes,
            int threads)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t x_div_y = x / y;
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);

  S2Status status(x);
  PiTable pi(z);
  SegmentedPiTable segmentedPi(isqrt(x), z, threads);
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_x13 = pi[x13];

  // This computes the 1st part of the C formula.
  // Find all special leaves of type:
  // x / (primes[b] * m) <= z.
  // m may be a prime or a square free number
  // who is coprime to the first b primes and
  // whose largest prime factor <= y.
  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
  for (int64_t b = k + 1; b <= pi_sqrtz; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t max_m = min(xp / prime, z);
    int64_t min_m = x / ipow<T>(prime, 3);
    min_m = max(min_m, prime, z / prime);

    sum += C<1>(xp, b, 1, min_m, max_m, b, primes, pi);

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
  for (; !segmentedPi.finished(); segmentedPi.next())
  {
    // Current segment [low, high[
    int64_t low = segmentedPi.low();
    int64_t high = segmentedPi.high();
    low = max(low, 1);
    int64_t x_div_low = x / low;
    int64_t x_div_high = x / high;

    int64_t min_prime1 = min(x_div_y / high, primes[pi_x_star]);
    int64_t min_prime2 = min(isqrt(low), primes[pi_x_star]);
    int64_t min_prime3 = min(iroot<3>(x_div_y), primes[pi_x_star]);

    int64_t min_b = max(k, pi_sqrtz);
    min_b = max(min_b, pi[min_prime1]);
    min_b = max(min_b, pi[min_prime2]);
    min_b = max(min_b, pi[min_prime3]);

    // x / (primes[i] * primes[i+1]) >= low
    // primes[i] * primes[i+1] <= x / low
    // primes[i] < sqrt(x / low)
    // primes[i+1] <= || >= sqrt(x / low)
    int64_t sqrt_low = min(isqrt(x_div_low), x13);
    int64_t max_b = pi[sqrt_low];

    if (primes[max_b] < max_a_prime &&
        (T) primes[max_b] * (T) primes[max_b + 1] > (T) x_div_low)
      max_b -= 1;

    min_b = min(min_b, pi_x_star + 1);
    max_b = max(max_b, pi_x_star);

    // 2nd part C formula: pi[sqrt(z)] < b <= pi[x_star]
    // A formula: pi[x_star] < b <= pi[x13]
    #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
    for (int64_t b = min_b + 1; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;

      if (b <= pi_x_star)
      {
        // This computes the 2nd part of the C formula
        int64_t max_m = min(xp / prime, x_div_low / prime, y);
        int64_t min_m = max(x / ipow<T>(prime, 3), x_div_high / prime);
        min_m = max(min_m, prime, z / prime);
        min_m = min(min_m, max_m);

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
      }
      else
      {
        // This computes the A formula
        int64_t min_2nd_prime = min(max_a_prime, x_div_high / prime);
        int64_t i = pi[min_2nd_prime] + 1;
        i = max(i, b + 1);
        int64_t max_2nd_prime = min(isqrt(xp), x_div_low / prime);
        int64_t max_i = pi[max_2nd_prime];

        // x / (p * q) >= y
        for (; i <= max_i; i++)
        {
          int64_t xpq = fast_div64(xp, primes[i]);
          if (xpq < y)
            break;
          sum += segmentedPi[xpq];
        }

        // x / (p * q) < y
        for (; i <= max_i; i++)
        {
          int64_t xpq = fast_div64(xp, primes[i]);
          sum += segmentedPi[xpq] * 2;
        }
      }

      if (is_print())
        status.print(b, pi_x13);
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
  print(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_c_prime = y;
  int64_t max_prime = max(max_a_prime, max_c_prime);
  auto primes = generate_primes<int32_t>(max_prime);
  int64_t ac = AC_OpenMP((intfast64_t) x, y, z, k, x_star, max_a_prime, primes, threads);

  print("A + C", ac, time);
  return ac;
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
  print(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_c_prime = y;
  int64_t max_prime = max(max_a_prime, max_c_prime);
  int128_t ac;

  // uses less memory
  if (max_prime <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(max_prime);
    ac = AC_OpenMP((intfast128_t) x, y, z, k, x_star, max_a_prime, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(max_prime);
    ac = AC_OpenMP((intfast128_t) x, y, z, k, x_star, max_a_prime,  primes, threads);
  }

  print("A + C", ac, time);
  return ac;
}

#endif

} // namespace

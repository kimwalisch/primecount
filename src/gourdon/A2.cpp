///
/// @file  A2.cpp
/// @brief Simple demonstration implementation of the A(x, y) formula
///        in Xavier Gourdon's prime counting algorithm.
///        In this version the memory usage has been reduced from
///        O(x^(1/2)) to O(z) by segmenting the pi[x] lookup table.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
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

template <typename T, typename Primes>
T A_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t x_star,
           Primes& primes,
           int threads)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  int64_t max_prime = primes.back();
  threads = ideal_num_threads(threads, x13, thread_threshold);

  PiTable pi(isqrt(x / x_star));
  SegmentedPiTable segmentedPi(isqrt(x), z);
  int64_t pi_x13 = pi[x13];
  S2Status status(x);

  // while (low <= sqrt(x))
  for (; !segmentedPi.finished(); segmentedPi.next())
  {
    // Current segment [low, high[
    int64_t low = segmentedPi.low();
    int64_t high = segmentedPi.high();
    low = max(low, 1);
    int64_t x_div_low = x / low;
    int64_t x_div_high = x / high;

    // x / (primes[i] * primes[i+1]) >= low
    // primes[i] * primes[i+1] <= x / low
    // primes[i] < sqrt(x / low)
    // primes[i+1] <= || >= sqrt(x / low)
    int64_t sqrt_low = min(isqrt(x_div_low), x13);
    int64_t max_b = pi[sqrt_low];
    if (primes[max_b] < max_prime &&
        (T) primes[max_b] * (T) primes[max_b + 1] > (T) x_div_low)
      max_b -= 1;

    #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
    for (int64_t b = pi[x_star] + 1; b <= max_b; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t min_2nd_prime = min(max_prime, x_div_high / prime);
      int64_t j = pi[min_2nd_prime] + 1;
      j = max(j, b + 1);
      int64_t max_2nd_prime = min(isqrt(xp), x_div_low / prime);
      int64_t max_j = pi[max_2nd_prime];

      // x / (p * q) >= y
      for (; j <= max_j; j++)
      {
        int64_t xpq = fast_div64(xp, primes[j]);
        if (xpq < y)
          break;
        sum += segmentedPi[xpq];
      }

      // x / (p * q) < y
      for (; j <= max_j; j++)
      {
        int64_t xpq = fast_div64(xp, primes[j]);
        sum += segmentedPi[xpq] * 2;
      }

      if (is_print())
        status.print(b, pi_x13);
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t A(int64_t x,
          int64_t y,
          int64_t z,
          int threads)
{
  print("");
  print("=== A(x, y) ===");
  print(x, y, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_prime = (int64_t) isqrt(x / x_star);

  auto primes = generate_primes<int32_t>(max_prime);
  int64_t a = A_OpenMP((intfast64_t) x, y, z, x_star, primes, threads);

  print("A", a, time);
  return a;
}

#ifdef HAVE_INT128_T

int128_t A(int128_t x,
           int64_t y,
           int64_t z,
           int threads)
{
  print("");
  print("=== A(x, y) ===");
  print(x, y, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_prime = (int64_t) isqrt(x / x_star);
  int128_t a;

  // uses less memory
  if (max_prime <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(max_prime);
    a = A_OpenMP((intfast128_t) x, y, z, x_star, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(max_prime);
    a = A_OpenMP((intfast128_t) x, y, z, x_star, primes, threads);
  }

  print("A", a, time);
  return a;
}

#endif

} // namespace

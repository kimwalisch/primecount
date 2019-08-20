///
/// @file  C4_libdivide.cpp
/// @brief Implementation of the C(x, y) formula in Xavier Gourdon's
///        prime counting algorithm. In this version the memory usage
///        has been reduced from O(x^(1/2)) to O(z) by segmenting
///        the pi[x] lookup table. In each segment we process the
///        leaves that satisfy: low <= x / (prime1 * prime2) < high.
///
///        This is an optimized version of C(x, y) which uses
///        libdivide. libdivide allows to replace expensive integer
///        divsion instructions by a sequence of shift, add and
///        multiply instructions that will calculate the integer
///        division much faster.
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
#include <libdivide.h>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <S2Status.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

template <typename T>
bool is_libdivide(T x)
{
  return x <= numeric_limits<uint64_t>::max();
}

using fastdiv_t = libdivide::branchfree_divider<uint64_t>;

template <typename Primes>
vector<fastdiv_t>
libdivide_vector(Primes& primes)
{
  vector<fastdiv_t> fastdiv(1);
  fastdiv.insert(fastdiv.end(), primes.begin() + 1, primes.end());
  return fastdiv;
}

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
T C_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t k,
           Primes& primes,
           int threads)
{
  T sum = 0;
  int64_t x_div_y = x / y;
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);
  auto fastdiv = libdivide_vector(primes);

  S2Status status(x);
  PiTable pi(z);
  SegmentedPiTable segmentedPi(isqrt(x), z, threads);
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_y = pi[y];

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
      status.print(b, pi_x_star);
  }

  // Find all special leaves of type:
  // z < x / (primes[b] * primes[i]) <= x^(1/2)
  // with z^(1/2) < primes[b] < x_star.
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

    int64_t max_b = pi[x_star];
    int64_t min_b = std::max(k, pi_sqrtz);
    min_b = std::max(min_b, pi[min_prime1]);
    min_b = std::max(min_b, pi[min_prime2]);
    min_b = std::max(min_b, pi[min_prime3]);

    #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
    for (int64_t b = min_b + 1; b <= pi_x_star; b++)
    {
      int64_t prime = primes[b];
      T xp = x / prime;
      int64_t max_m = min(xp / prime, x_div_low / prime, y);
      int64_t min_m = max(x / ipow<T>(prime, 3), x_div_high / prime);
      min_m = max(min_m, prime, z / prime);
      min_m = min(min_m, max_m);

      int64_t i = pi[max_m];
      int64_t pi_min_m = pi[min_m];
      int64_t min_clustered = (int64_t) isqrt(xp);
      min_clustered = in_between(min_m, min_clustered, max_m);
      int64_t pi_min_clustered = pi[min_clustered];

      // x / (primes[b] * prime[i]) < low
      // Nothing todo, these leaves were already
      // processed in the previous segment.
      if (pi_min_m + 1 <= pi_y &&
          x_div_low < (T) prime * primes[pi_min_m + 1])
        continue;

      if (is_libdivide(xp))
      {
        // Find all clustered easy leaves where
        // successive leaves are identical.
        // n = primes[b] * primes[i]
        // Which satisfy: n > z && primes[i] <= y
        while (i > pi_min_clustered)
        {
          int64_t xpq = (uint64_t) xp / fastdiv[i];
          int64_t phi_xpq = segmentedPi[xpq] - b + 2;
          int64_t xpq2 = (uint64_t) xp / fastdiv[b + phi_xpq - 1];
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
          int64_t xpq = (uint64_t) xp / fastdiv[i];
          sum += segmentedPi[xpq] - b + 2;
        }
      }
      else
      {      
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
  print(x, y, z, k, threads);

  double time = get_time();
  auto primes = generate_primes<int32_t>(y);
  int64_t c = C_OpenMP((intfast64_t) x, y, z, k, primes, threads);

  print("C", c, time);
  return c;
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
  print(x, y, z, k, threads);

  double time = get_time();
  int128_t c;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    c = C_OpenMP((intfast128_t) x, y, z, k, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    c = C_OpenMP((intfast128_t) x, y, z, k, primes, threads);
  }

  print("C", c, time);
  return c;
}

#endif

} // namespace

///
/// @file  C3.cpp
/// @brief Simple demonstration implementation of the C(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation uses O(x^(1/2)) memory instead of O(x^(1/3))
///        in order to simplify the implementation.
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

    if (next_m > min_m)
    {
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
  T y2 = y * (T) y;
  int64_t x_star = max(iroot<4>(x), x / y2);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);

  PiTable pi(isqrt(x));
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  S2Status status(x);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
  for (int64_t b = k + 1; b <= pi_x_star; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t max_m = min(xp / prime, z);
    int64_t min_m = x / ipow<T>(prime, 3);
    min_m = max(min_m, prime, z / prime);

    if (b <= pi_sqrtz)
      sum += C<1>(xp, b, 1, min_m, max_m, b, primes, pi);
    else
    {
      // Above sqrt(z) m is composed of a single
      // prime and that prime must be <= y
      max_m = min(max_m, y);
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
        int64_t m = primes[i];
        int64_t xpm = fast_div64(xp, m);
        int64_t phi_xpm = pi[xpm] - b + 2;
        int64_t m2 = primes[b + phi_xpm - 1];
        int64_t xpm2 = fast_div64(xp, m2);
        int64_t i2 = pi[xpm2];
        sum += phi_xpm * (i - i2);
        i = i2;
      }

      // Find all sparse easy leaves where
      // successive leaves are different.
      // n = primes[b] * primes[i]
      // Which satisfy: n > z && primes[i] <= y
      for (; i > pi_min_m; i--)
      {
        int64_t m = primes[i];
        int64_t xpm = fast_div64(xp, m);
        sum += pi[xpm] - b + 2;
      }
    }

    if (is_print())
      status.print(b, pi_x_star);
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

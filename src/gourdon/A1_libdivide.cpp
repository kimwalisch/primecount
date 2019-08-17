///
/// @file  A1_libdivide.cpp
/// @brief Simple demonstration implementation of the A(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation uses O(x^(1/2)) memory instead of O(x^(1/3))
///        in order to simplify the implementation.
///
///        This is an optimized version of A(x, y) which uses
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

#include <gourdon.hpp>
#include <PiTable.hpp>
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

template <typename T, typename Primes>
T A_OpenMP(T x,
           int64_t y,
           int64_t x_star,
           Primes& primes,
           int threads)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);
  auto fastdiv = libdivide_vector(primes);

  PiTable pi(isqrt(x));
  int64_t pi_x13 = pi[x13];
  S2Status status(x);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
  for (int64_t b = pi[x_star] + 1; b <= pi_x13; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t j = b + 1;
    int64_t max_j = pi[isqrt(xp)];

    if (is_libdivide(xp))
    {
      // x / (p * q) >= y
      for (; j <= max_j; j++)
      {
        int64_t xpq = (uint64_t) xp / fastdiv[j];
        if (xpq < y)
          break;
        sum += pi[xpq];
      }

      // x / (p * q) < y
      for (; j <= max_j; j++)
      {
        int64_t xpq = (uint64_t) xp / fastdiv[j];
        sum += pi[xpq] * 2;
      }
    }
    else
    {
      // x / (p * q) >= y
      for (; j <= max_j; j++)
      {
        int64_t xpq = fast_div64(xp, primes[j]);
        if (xpq < y)
          break;
        sum += pi[xpq];
      }

      // x / (p * q) < y
      for (; j <= max_j; j++)
      {
        int64_t xpq = fast_div64(xp, primes[j]);
        sum += pi[xpq] * 2;
      }
    }

    if (is_print())
      status.print(b, pi_x13);
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t A(int64_t x,
          int64_t y,
          int threads)
{
  print("");
  print("=== A(x, y) ===");
  print(x, y, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_prime = (int64_t) isqrt(x / x_star);

  auto primes = generate_primes<int32_t>(max_prime);
  int64_t a = A_OpenMP((intfast64_t) x, y, x_star, primes, threads);

  print("A", a, time);
  return a;
}

#ifdef HAVE_INT128_T

int128_t A(int128_t x,
           int64_t y,
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
    a = A_OpenMP((intfast128_t) x, y, x_star, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(max_prime);
    a = A_OpenMP((intfast128_t) x, y, x_star, primes, threads);
  }

  print("A", a, time);
  return a;
}

#endif

} // namespace

///
/// @file  C1.cpp
/// @brief Simple demonstration implementation of the C(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation uses O(x^(1/2)) memory instead of O(x^(1/3))
///        in order to simplify the implementation.
///
///        Currently this implementation is quite slow when compared
///        to Xavier Gourdon's fastpix11.exe binary. This
///        implementation is slow mainly because it iterates over all
///        integers and for each integer checks whether it is coprime
///        to the first b primes. It is possible to iterate only over
///        the square free integers which are coprime to the first b
///        primes which is obviously much faster (see C2.cpp).
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

template <typename T, typename Primes>
T C_OpenMP(T x,
           int64_t y,
           int64_t z,
           int64_t k,
           Primes& primes,
           int threads)
{
  T sum = 0;
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);

  PiTable pi(isqrt(x));
  int64_t pi_x_star = pi[x_star];
  S2Status status(x);

  auto mu = generate_moebius(z);
  auto lpf = generate_lpf(z);
  auto mpf = generate_mpf(z);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(-: sum)
  for (int64_t b = k + 1; b <= pi_x_star; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t m = min(xp / prime, z);
    int64_t min_m = x / ipow<T>(prime, 3);
    min_m = max(min_m, prime, z / prime);

    for (; m > min_m; m--)
    {
      if (lpf[m] > prime && 
          mpf[m] <= y)
      {
        int64_t xpm = fast_div64(xp, m);
        sum -= mu[m] * (pi[xpm] - b + 2);
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

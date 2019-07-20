///
/// @file  C.cpp
/// @brief Simple demonstration implementation of the C(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation uses O(x^(1/2)) memory instead of O(x^(1/3))
///        in order to simplify the implementation.
///
///        Currently this implementation is very slow, it is more than
///        3x slower than Xavier Gourdon's fastpix11.exe binary. There
///        are 2 obvious options to improve this implementation:
///
///        1) Implement clustered easy leaves like in S2_easy.cpp.
///        2) Use FactorTable lookup table instead of mu, lpf, mpf.
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
  T y2 = y * (T) y;
  int64_t x_star = max(iroot<4>(x), x / y2);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);

  PiTable pi(isqrt(x));
  int64_t pi_x_star = pi[x_star];
  S2Status status(x);

  auto mu = generate_moebius(z);
  auto lpf = generate_lpf(z);
  auto mpf = generate_mpf(z);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: sum)
  for (int64_t b = k + 1; b <= pi_x_star; b++)
  {
    int64_t prime = primes[b];
    T x2 = x / prime;
    int64_t min_trivial = min(x2 / prime, z);
    int64_t min_sparse = x / ipow<T>(prime, 3);
    min_sparse = max(prime, min_sparse);
    int64_t m = min_trivial;

    for (; m > min_sparse; m--)
    {
      if (mu[m] != 0 &&
          lpf[m] > prime && 
          mpf[m] <= y)
      {
        int64_t xn = fast_div64(x2, m);
        sum += mu[m] * (pi[xn] - b + 2);
      }
    }

    if (is_print())
      status.print(b, pi_x_star);
  }

  return -sum;
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
  auto primes = generate_primes<int32_t>(z);
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

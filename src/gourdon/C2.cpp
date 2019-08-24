///
/// @file  C2.cpp
/// @brief Simple demonstration implementation of the C(x, y) formula
///        in Xavier Gourdon's prime counting algorithm. This
///        implementation uses O(x^(1/2)) memory instead of O(x^(1/3))
///        in order to simplify the implementation.
///
///        This implementation runs about 5x faster than C1.cpp
///        because it only iterates over the square free integers
///        which are coprime to the first b primes. C1.cpp iterates
///        over all integer and for each integer checks whether it is
///        coprime to the first b primes.
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
    uint64_t b,
    uint64_t i,
    uint64_t m,
    uint64_t min_m,
    uint64_t max_m,
    Primes& primes,
    PiTable& pi)
{
  T sum = 0;

  for (i++; i < primes.size(); i++)
  {
    // next m may be 128-bit
    T m128 = (T) m * primes[i];
    if (m128 > max_m)
      return sum;

    uint64_t m64 = (uint64_t) m128;
    if (m64 > min_m) {
      uint64_t xpm = fast_div64(xp, m64);
      sum += MU * (pi[xpm] - b + 2);
    }

    sum += C<-MU>(xp, b, i, m64, min_m, max_m, primes, pi);
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
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x_star, thread_threshold);

  PiTable pi(isqrt(x));
  int64_t pi_x_star = pi[x_star];
  S2Status status(x);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(-: sum)
  for (int64_t b = k + 1; b <= pi_x_star; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t max_m = min(xp / prime, z);
    T min_m128 = max(x / ipow<T>(prime, 3), z / prime);
    int64_t min_m = min(min_m128, max_m);

    if (min_m < max_m)
      sum -= C<-1>(xp, b, b, 1, min_m, max_m, primes, pi);

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
  print_gourdon(x, y, z, k, threads);

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
  print_gourdon(x, y, z, k, threads);

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

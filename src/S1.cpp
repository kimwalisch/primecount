///
/// @file  S1.cpp
/// @brief Calculate the contribution of the ordinary leaves in the
///        Lagarias-Miller-Odlyzko and Deleglise-Rivat prime counting
///        algorithms. Note that in the original paper from
///        Lagarias-Miller-Odlyzko the ordinary leaves were named S1
///        however in many newer papers the ordinary leaves are named
///        S0. In primecount the ordinary leaves are named S1.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <pod_vector.hpp>
#include <print.hpp>
#include <S.hpp>

#include <stdint.h>

using std::numeric_limits;
using namespace primecount;

namespace {

/// Recursively iterate over the square free numbers coprime
/// to the first b primes and calculate the sum of the
/// ordinary leaves. This algorithm is described in section
/// 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March
/// 2015.
///
template <int MU, typename T, typename vect>
T S1_thread(T x,
            int64_t y,
            uint64_t b,
            int64_t c,
            T square_free,
            const vect& primes)
{
  T s1 = 0;

  for (b++; b < primes.size(); b++)
  {
    T next = square_free * primes[b];
    if (next > y) break;
    s1 += MU * phi_tiny(x / next, c);
    s1 += S1_thread<-MU>(x, y, b, c, next, primes);
  }

  return s1;
}

/// Parallel computation of the ordinary leaves.
/// Run time: O(y * log(log(y)))
/// Memory usage: O(y / log(y))
///
template <typename X, typename Y>
X S1_OpenMP(X x,
            Y y,
            int64_t c,
            int threads)
{
  threads = ideal_num_threads(threads, y);

  auto primes = generate_primes<Y>(y);
  int64_t pi_y = primes.size() - 1;
  X s1 = phi_tiny(x, c);

  #pragma omp parallel for schedule(static, 1) num_threads(threads) reduction (+: s1)
  for (int64_t b = c + 1; b <= pi_y; b++)
  {
    s1 -= phi_tiny(x / primes[b], c);
    s1 += S1_thread<1>(x, y, b, c, (X) primes[b], primes);
  }

  return s1;
}

} // namespace

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int threads,
           bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== S1(x, y) ===");
    print_vars(x, y, c, threads);
  }

  double time = get_time();
  int64_t s1 = S1_OpenMP(x, y, c, threads);

  if (is_print)
    print("S1", s1, time);

  return s1;
}

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int threads,
            bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== S1(x, y) ===");
    print_vars(x, y, c, threads);
  }

  double time = get_time();
  int128_t s1;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
    s1 = S1_OpenMP(x, (uint32_t) y, c, threads);
  else
    s1 = S1_OpenMP(x, y, c, threads);

  if (is_print)
    print("S1", s1, time);

  return s1;
}

#endif

} // namespace

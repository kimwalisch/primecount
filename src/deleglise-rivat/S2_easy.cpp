///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <int128.hpp>
#include <min_max.hpp>
#include <pmath.hpp>
#include <S2Status.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace S2_easy {

/// Calculate the contribution of the clustered easy leaves and
/// sparse easy leaves corresponding to primes[b] i.e.:
/// n = primes[b] * primes[l] with x / n <= y.
/// @param x2 = x / primes[b].
///
template <typename T1, typename T2, typename T3, typename R = T1>
R S2_easy(T1 x2,
          int64_t y,
          int64_t z,
          int64_t b,
          T2& pi,
          vector<T3>& primes)
{
  R s2_easy = 0;
  int64_t prime = primes[b];
  int64_t min_trivial_leaf = min(x2 / prime, y);
  int64_t min_clustered_easy_leaf = (int64_t) isqrt(x2);
  int64_t min_sparse_easy_leaf = z / prime;
  int64_t min_hard_leaf = max(y / prime, prime);

  min_sparse_easy_leaf = max(min_sparse_easy_leaf, min_hard_leaf);
  min_clustered_easy_leaf = max(min_clustered_easy_leaf, min_hard_leaf);
  int64_t l = pi[min_trivial_leaf];

  // Find all clustered easy leaves:
  // n = primes[b] * primes[l]
  // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
  // where phi(x / n, b - 1) = pi(x / n) - b + 2
  while (primes[l] > min_clustered_easy_leaf)
  {
    int64_t xn = (int64_t) (x2 / primes[l]);
    int64_t phi_xn = pi[xn] - b + 2;
    int64_t last_prime = primes[b + phi_xn - 1];
    int64_t xm = max((int64_t) (x2 / last_prime), min_clustered_easy_leaf);
    int64_t l2 = pi[xm];
    s2_easy += phi_xn * (l - l2);
    l = l2;
  }

  // Find all sparse easy leaves:
  // n = primes[b] * primes[l]
  // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
  for (; primes[l] > min_sparse_easy_leaf; l--)
  {
    int64_t xn = (int64_t) (x2 / primes[l]);
    s2_easy += pi[xn] - b + 2;
  }

  return s2_easy;
}

/// Strategy pattern, chooses an S2_easy() implementation.
template <typename T1, typename T2, typename T3>
T1 S2_strategy(T1 x,
               int64_t y,
               int64_t z,
               int64_t b,
               T2& pi,
               vector<T3>& primes)
{
  T1 x2 = x / primes[b];

  // Avoid slow 128-bit divisions
  if (x2 <= std::numeric_limits<intfast64_t>::max())
    return S2_easy<T1>((intfast64_t) x2, y, z, b, pi, primes);

  return S2_easy(x2, y, z, b, pi, primes);
}

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T1, typename T2, typename T3>
T1 S2_easy(T1 x,
           int64_t y,
           int64_t z,
           int64_t c,
           T2& pi,
           vector<T3>& primes,
           int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_easy(x, y) ===" << endl;
    cout << "Computation of the easy special leaves" << endl;
  }

  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[iroot<3>(x)];
  T1 s2_easy = 0;
  double time = get_wtime();
  S2Status status;

  #pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+: s2_easy)
  for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
  {
    s2_easy += S2_strategy(x, y, z, b, pi, primes);

    if (print_status())
      status.print(b, pi_x13);
  }

  if (print_status())
    print_result("S2_easy", s2_easy, time);

  return s2_easy;
}

} // namespace S2_easy
} // namespace

namespace primecount {

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                vector<int32_t>& pi,
                vector<int32_t>& primes,
                int threads)
{
  return S2_easy::S2_easy((intfast64_t) x, y, z, c, pi, primes, threads);
}

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                PiTable& pi,
                vector<int32_t>& primes,
                int threads)
{
  return S2_easy::S2_easy((intfast64_t) x, y, z, c, pi, primes, threads);
}

#ifdef HAVE_INT128_T

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 PiTable& pi,
                 vector<uint32_t>& primes,
                 int threads)
{
  return S2_easy::S2_easy((intfast128_t) x, y, z, c, pi, primes, threads);
}

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 PiTable& pi,
                 vector<int64_t>& primes,
                 int threads)
{
  return S2_easy::S2_easy((intfast128_t) x, y, z, c, pi, primes, threads);
}

#endif

} // namespace primecount

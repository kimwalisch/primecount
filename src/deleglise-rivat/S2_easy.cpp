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
  T1 S2_total = 0;
  double time = get_wtime();
  S2Status status;

  #pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
  {
    T1 prime = primes[b];
    int64_t prime64 =  primes[b];
    int64_t min_trivial_leaf = min(x / (prime * prime), y);
    int64_t min_clustered_easy_leaf = isqrt(x / prime64);
    int64_t min_sparse_easy_leaf = z / prime64;
    int64_t min_hard_leaf = max(y / prime64, prime64);

    min_sparse_easy_leaf = max(min_sparse_easy_leaf, min_hard_leaf);
    min_clustered_easy_leaf = max(min_clustered_easy_leaf, min_hard_leaf);
    int64_t l = pi[min_trivial_leaf];
    T1 S2_result = 0;

    // Find all clustered easy leaves:
    // x / n <= y and phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (primes[l] > min_clustered_easy_leaf)
    {
      T1 n = prime * primes[l];
      int64_t xn = (int64_t) (x / n);
      int64_t phi_xn = pi[xn] - b + 2;
      T1 m = prime * primes[b + phi_xn - 1];
      int64_t xm = max((int64_t) (x / m), min_clustered_easy_leaf);
      int64_t l2 = pi[xm];
      T1 phi_factor = l - l2;
      S2_result += phi_xn * phi_factor;
      l = l2;
    }

    // Find all sparse easy leaves:
    // x / n <= y and phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; primes[l] > min_sparse_easy_leaf; l--)
    {
      T1 n = prime * primes[l];
      int64_t xn = (int64_t) (x / n);
      S2_result += pi[xn] - b + 2;
    }

    if (print_status())
      status.print(b, pi_x13);

    S2_total += S2_result;
  }

  if (print_status())
    print_result("S2_easy", S2_total, time);

  return S2_total;
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
  return S2_easy::S2_easy(x, y, z, c, pi, primes, threads);
}

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                PiTable& pi,
                vector<int32_t>& primes,
                int threads)
{
  return S2_easy::S2_easy(x, y, z, c, pi, primes, threads);
}

#ifdef HAVE_INT128_T

int128_t S2_easy(uint128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 PiTable& pi,
                 vector<uint32_t>& primes,
                 int threads)
{
  return S2_easy::S2_easy(x, y, z, c, pi, primes, threads);
}

int128_t S2_easy(uint128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 PiTable& pi,
                 vector<int64_t>& primes,
                 int threads)
{
  return S2_easy::S2_easy(x, y, z, c, pi, primes, threads);
}

#endif

} // namespace primecount

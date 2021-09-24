///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
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
#include <RelaxedAtomic.hpp>
#include <StatusS2.hpp>
#include <S.hpp>

#include <stdint.h>

using std::numeric_limits;
using namespace primecount;

namespace {

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T, typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 const Primes& primes,
                 int threads,
                 bool is_print)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);

  StatusS2 status(x);
  PiTable pi(y, threads);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  RelaxedAtomic<int64_t> min_b(max(c, pi_sqrty) + 1);

  // for (b = pi[sqrty] + 1; b <= pi_x13; b++)
  #pragma omp parallel num_threads(threads) reduction(+: sum)
  for (int64_t b = min_b++; b <= pi_x13; b = min_b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t min_trivial = min(xp / prime, y);
    int64_t min_clustered = (int64_t) isqrt(xp);
    int64_t min_sparse = z / prime;

    min_clustered = in_between(prime, min_clustered, y);
    min_sparse = in_between(prime, min_sparse, y);

    int64_t l = pi[min_trivial];
    int64_t pi_min_clustered = pi[min_clustered];
    int64_t pi_min_sparse = pi[min_sparse];

    // Find all clustered easy leaves where
    // successive leaves are identical.
    // pq = primes[b] * primes[l]
    // Which satisfy: pq > z && x / pq <= y
    // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
    while (l > pi_min_clustered)
    {
      int64_t xpq = fast_div64(xp, primes[l]);
      int64_t pi_xpq = pi[xpq];
      int64_t phi_xpq = pi_xpq - b + 2;
      int64_t xpq2 = fast_div64(xp, primes[pi_xpq + 1]);
      int64_t lmin = pi[xpq2];
      sum += phi_xpq * (l - lmin);
      l = lmin;
    }

    // Find all sparse easy leaves where
    // successive leaves are different.
    // pq = primes[b] * primes[l]
    // Which satisfy: pq > z && x / pq <= y
    // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
    for (; l > pi_min_sparse; l--)
    {
      int64_t xpq = fast_div64(xp, primes[l]);
      sum += pi[xpq] - b + 2;
    }

    #pragma omp master
    if (is_print)
      status.print(b, pi_x13);
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int threads,
                bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== S2_easy(x, y) ===");
    print_vars(x, y, c, threads);
  }

  double time = get_time();
  auto primes = generate_primes<uint32_t>(y);
  int64_t sum = S2_easy_OpenMP((uint64_t) x, y, z, c, primes, threads, is_print);

  if (is_print)
    print("S2_easy", sum, time);

  return sum;
}

#ifdef HAVE_INT128_T

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads,
                 bool is_print)
{
  if (is_print)
  {
    print("");
    print("=== S2_easy(x, y) ===");
    print_vars(x, y, c, threads);
  }

  double time = get_time();
  int128_t sum;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    sum = S2_easy_OpenMP((uint128_t) x, y, z, c, primes, threads, is_print);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    sum = S2_easy_OpenMP((uint128_t) x, y, z, c, primes, threads, is_print);
  }

  if (is_print)
    print("S2_easy", sum, time);

  return sum;
}

#endif

} // namespace

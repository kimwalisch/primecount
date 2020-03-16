///
/// @file  S2_easy_libdivide.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP.
///        This is an optimized version of S2_easy(x, y) which uses
///        libdivide. libdivide allows to replace expensive integer
///        divsion instructions by a sequence of shift, add and
///        multiply instructions that will calculate the integer
///        division much faster.
///
///        This implementation is based on the paper:
///        Tomás Oliveira e Silva, Computing pi(x): the combinatorial
///        method, Revista do DETUA, vol. 4, no. 6, March 2006,
///        pp. 759-768.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
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
#include <S.hpp>

#include <libdivide.h>
#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

using fastdiv_t = libdivide::branchfree_divider<uint64_t>;

template <typename Primes>
vector<fastdiv_t>
libdivide_primes(const Primes& primes)
{
  vector<fastdiv_t> lprimes(1);
  lprimes.insert(lprimes.end(), primes.begin() + 1, primes.end());
  return lprimes;
}

/// xp < 2^64
template <typename T, typename LibdividePrimes>
T S2_easy_64(T xp128,
             uint64_t prime,
             uint64_t y,
             uint64_t z,
             uint64_t b,
             const PiTable& pi,
             const LibdividePrimes& primes)
{
  uint64_t xp = (uint64_t) xp128;
  uint64_t min_trivial = min(xp / prime, y);
  uint64_t min_clustered = isqrt(xp);
  uint64_t min_sparse = z / prime;
  min_clustered = in_between(prime, min_clustered, y);
  min_sparse = in_between(prime, min_sparse, y);
  uint64_t l = pi[min_trivial];
  uint64_t pi_min_clustered = pi[min_clustered];
  uint64_t pi_min_sparse = pi[min_sparse];

  T s2_easy = 0;

  // Find all clustered easy leaves where
  // successive leaves are identical.
  // pq = primes[b] * primes[l]
  // Which satisfy: pq > z && x / pq <= y
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  while (l > pi_min_clustered)
  {
    uint64_t xpq = xp / primes[l];
    uint64_t phi_xpq = pi[xpq] - b + 2;
    uint64_t xpq2 = xp / primes[b + phi_xpq - 1];
    uint64_t l2 = pi[xpq2];
    s2_easy += phi_xpq * (l - l2);
    l = l2;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // pq = primes[b] * primes[l]
  // Which satisfy: pq > z && x / pq <= y
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  for (; l > pi_min_sparse; l--)
  {
    uint64_t xpq = xp / primes[l];
    s2_easy += pi[xpq] - b + 2;
  }

  return s2_easy;
}

/// xp >= 2^64
template <typename T, typename Primes>
T S2_easy_128(T xp,
              uint64_t prime,
              uint64_t y,
              uint64_t z,
              uint64_t b,
              const PiTable& pi,
              const Primes& primes)
{
  uint64_t min_trivial = min(xp / prime, y);
  uint64_t min_clustered = (uint64_t) isqrt(xp);
  uint64_t min_sparse = z / prime;
  min_clustered = in_between(prime, min_clustered, y);
  min_sparse = in_between(prime, min_sparse, y);
  uint64_t l = pi[min_trivial];
  uint64_t pi_min_clustered = pi[min_clustered];
  uint64_t pi_min_sparse = pi[min_sparse];

  T s2_easy = 0;

  // Find all clustered easy leaves where
  // successive leaves are identical.
  // pq = primes[b] * primes[l]
  // Which satisfy: pq > z && x / pq <= y
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  while (l > pi_min_clustered)
  {
    uint64_t xpq = fast_div64(xp, primes[l]);
    uint64_t phi_xpq = pi[xpq] - b + 2;
    uint64_t xpq2 = fast_div64(xp, primes[b + phi_xpq - 1]);
    uint64_t l2 = pi[xpq2];
    s2_easy += phi_xpq * (l - l2);
    l = l2;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // pq = primes[b] * primes[l]
  // Which satisfy: pq > z && x / pq <= y
  // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
  for (; l > pi_min_sparse; l--)
  {
    uint64_t xpq = fast_div64(xp, primes[l]);
    s2_easy += pi[xpq] - b + 2;
  }

  return s2_easy;
}

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
template <typename T, typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 const Primes& primes,
                 int threads)
{
  T s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  threads = ideal_num_threads(threads, x13, 1000);
  auto lprimes = libdivide_primes(primes);

  PiTable pi(y);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  S2Status status(x);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: s2_easy)
  for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
  {
    // Unsigned integer division is usually slightly
    // faster than signed integer division
    using UT = typename make_unsigned<T>::type;

    int64_t prime = primes[b];
    UT xp = x / prime;

    if (xp <= numeric_limits<uint64_t>::max())
      s2_easy += S2_easy_64(xp, prime, y, z, b, pi, lprimes);
    else
      s2_easy += S2_easy_128(xp, prime, y, z, b, pi, primes);

    if (is_print())
      status.print(b, pi_x13);
  }

  return s2_easy;
}

} // namespace

namespace primecount {

int64_t S2_easy(int64_t x,
                int64_t y,
                int64_t z,
                int64_t c,
                int threads)
{
#ifdef ENABLE_MPI
  if (mpi_num_procs() > 1)
    return S2_easy_mpi(x, y, z, c, threads);
#endif

  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print_vars(x, y, c, threads);

  double time = get_time();
  auto primes = generate_primes<uint32_t>(y);
  int64_t s2_easy = S2_easy_OpenMP(x, y, z, c, primes, threads);

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#ifdef HAVE_INT128_T

int128_t S2_easy(int128_t x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads)
{
#ifdef ENABLE_MPI
  if (mpi_num_procs() > 1)
    return S2_easy_mpi(x, y, z, c, threads);
#endif

  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print_vars(x, y, c, threads);

  double time = get_time();
  int128_t s2_easy;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    s2_easy = S2_easy_OpenMP(x, y, z, c, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    s2_easy = S2_easy_OpenMP(x, y, z, c, primes, threads);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

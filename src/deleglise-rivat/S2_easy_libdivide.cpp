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
#include <S2.hpp>

#include <libdivide.h>
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

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
template <typename T, typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 Primes& primes,
                 int threads)
{
  T s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  threads = ideal_num_threads(threads, x13, 1000);
  auto fastdiv = libdivide_vector(primes);

  PiTable pi(y);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  S2Status status(x);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: s2_easy)
  for (int64_t b = max(c, pi_sqrty) + 1; b <= pi_x13; b++)
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

    if (is_libdivide(xp))
    {
      // Find all clustered easy leaves where
      // successive leaves are identical.
      // pq = primes[b] * primes[l]
      // Which satisfy: pq > z && x / pq <= y
      // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
      while (l > pi_min_clustered)
      {
        int64_t xpq = (uint64_t) xp / fastdiv[l];
        int64_t phi_xpq = pi[xpq] - b + 2;
        int64_t xpq2 = (uint64_t) xp / fastdiv[b + phi_xpq - 1];
        int64_t l2 = pi[xpq2];
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
        int64_t xpq = (uint64_t) xp / fastdiv[l];
        s2_easy += pi[xpq] - b + 2;
      }
    }
    else
    {
      // Find all clustered easy leaves where
      // successive leaves are identical.
      // pq = primes[b] * primes[l]
      // Which satisfy: pq > z && x / pq <= y
      // where phi(x / pq, b - 1) = pi(x / pq) - b + 2
      while (l > pi_min_clustered)
      {
        int64_t xpq = fast_div64(xp, primes[l]);
        int64_t phi_xpq = pi[xpq] - b + 2;
        int64_t xpq2 = fast_div64(xp, primes[b + phi_xpq - 1]);
        int64_t l2 = pi[xpq2];
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
        int64_t xpq = fast_div64(xp, primes[l]);
        s2_easy += pi[xpq] - b + 2;
      }
    }

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
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return S2_easy_mpi(x, y, z, c, threads);
#endif

  print("");
  print("=== S2_easy(x, y) ===");
  print("Computation of the easy special leaves");
  print_vars(x, y, c, threads);

  double time = get_time();
  auto primes = generate_primes<int32_t>(y);
  int64_t s2_easy = S2_easy_OpenMP((intfast64_t) x, y, z, c, primes, threads);

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
#ifdef HAVE_MPI
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
    s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, primes, threads);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

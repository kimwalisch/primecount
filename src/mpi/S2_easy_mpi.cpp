///
/// @file  S2_easy_mpi.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <mpi_reduce_sum.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <S2Status.hpp>
#include <fast_div.hpp>
#include <print.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T, typename Primes>
T S2_easy_mpi_master(T x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     Primes& primes,
                     int threads)
{
  T s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);

  PiTable pi(y);
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t pi_x13 = pi[x13];
  S2Status status(x);

  int proc_id = mpi_proc_id();
  int procs = mpi_num_procs();

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: s2_easy)
  for (int64_t b = max(c, pi_sqrty) + 1 + proc_id; b <= pi_x13; b += procs)
  {
    int64_t prime = primes[b];
    T x2 = x / prime;
    int64_t min_trivial = min(x2 / prime, y);
    int64_t min_clustered = (int64_t) isqrt(x2);
    int64_t min_sparse = z / prime;

    min_clustered = in_between(prime, min_clustered, y);
    min_sparse = in_between(prime, min_sparse, y);

    int64_t l = pi[min_trivial];
    int64_t pi_min_clustered = pi[min_clustered];
    int64_t pi_min_sparse = pi[min_sparse];

    // Find all clustered easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (l > pi_min_clustered)
    {
      int64_t xn = fast_div64(x2, primes[l]);
      int64_t phi_xn = pi[xn] - b + 2;
      int64_t xm = fast_div64(x2, primes[b + phi_xn - 1]);
      int64_t l2 = pi[xm];
      s2_easy += phi_xn * (l - l2);
      l = l2;
    }

    // Find all sparse easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; l > pi_min_sparse; l--)
    {
      int64_t xn = fast_div64(x2, primes[l]);
      s2_easy += pi[xn] - b + 2;
    }

    if (is_print())
      status.print(b, pi_x13);
  }

  s2_easy = mpi_reduce_sum(s2_easy);

  return s2_easy;
}

} // namespace

namespace primecount {

int64_t S2_easy_mpi(int64_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads)
{
  print("");
  print("=== S2_easy_mpi(x, y) ===");
  print("Computation of the easy special leaves");
  print(x, y, c, threads);

  double time = get_time();
  auto primes = generate_primes<int32_t>(y);
  int64_t s2_easy = S2_easy_mpi_master((intfast64_t) x, y, z, c, primes, threads);

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#ifdef HAVE_INT128_T

int128_t S2_easy_mpi(int128_t x,
                     int64_t y,
                     int64_t z,
                     int64_t c,
                     int threads)
{
  print("");
  print("=== S2_easy_mpi(x, y) ===");
  print("Computation of the easy special leaves");
  print(x, y, c, threads);

  double time = get_time();
  int128_t s2_easy;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    s2_easy = S2_easy_mpi_master((intfast128_t) x, y, z, c, primes, threads);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    s2_easy = S2_easy_mpi_master((intfast128_t) x, y, z, c, primes, threads);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

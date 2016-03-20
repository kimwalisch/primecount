///
/// @file  S2_easy_mpi_libdivide.cpp
/// @brief This is an optimized version of S2_easy which uses
///        libdivide. libdivide allows to replace expensive integer
///        divides with comparatively cheap multiplication and
///        bitshifts.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <int128.hpp>
#include <LibdividePrimes.hpp>
#include <min_max.hpp>
#include <mpi_reduce_sum.hpp>
#include <pmath.hpp>
#include <S2Status.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T1, typename T2>
T1 S2_easy_mpi_master(T1 x,
                      int64_t y,
                      int64_t z,
                      int64_t c,
                      LibdividePrimes<T2>& primes,
                      int threads)
{
  T1 s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = validate_threads(threads, x13, thread_threshold);

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
    T1 x2 = x / prime;
    int64_t min_trivial = min(x2 / prime, y);
    int64_t min_clustered = (int64_t) isqrt(x2);
    int64_t min_sparse = z / prime;
    int64_t min_hard = max(y / prime, prime);

    min_clustered = in_between(min_hard, min_clustered, y);
    min_sparse = in_between(min_hard, min_sparse, y);

    int64_t l = pi[min_trivial];
    int64_t pi_min_clustered = pi[min_clustered];
    int64_t pi_min_sparse = pi[min_sparse];

    T1 sum = 0;

    // Find all clustered easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (l > pi_min_clustered)
    {
      int64_t xn = (int64_t) primes.libdivide(x2, l);
      int64_t phi_xn = pi[xn] - b + 2;
      int64_t xm = (int64_t) primes.libdivide(x2, b + phi_xn - 1);
      xm = max(xm, min_clustered);
      int64_t l2 = pi[xm];
      sum += phi_xn * (l - l2);
      l = l2;
    }

    // Find all sparse easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; l > pi_min_sparse; l--)
    {
      int64_t xn = (int64_t) primes.libdivide(x2, l);
      sum += pi[xn] - b + 2;
    }

    if (print_status())
      status.print(b, pi_x13);

    s2_easy += sum;
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

  double time = get_wtime();
  LibdividePrimes<int32_t> primes(y);
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

  double time = get_wtime();
  int128_t s2_easy;

  // uses less memory
  if (y <= std::numeric_limits<uint32_t>::max())
  {
    LibdividePrimes<uint32_t> primes(y);
    s2_easy = S2_easy_mpi_master((intfast128_t) x, y, z, c, primes, threads);
  }
  else
  {
    LibdividePrimes<int64_t> primes(y);
    s2_easy = S2_easy_mpi_master((intfast128_t) x, y, z, c, primes, threads);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace primecount

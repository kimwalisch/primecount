///
/// @file  Phi0_mpi.cpp
/// @brief Implementation of the Phi0 formula (from Xavier Gourdon's
///        algorithm) that has been distributed using MPI (Message
///        Passing Interface) and multi-threaded using OpenMP.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <mpi_reduce_sum.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <print.hpp>

#include <stdint.h>
#include <vector>

using namespace std;
using namespace primecount;

namespace {

/// Recursively iterate over the square free numbers coprime
/// to the first b primes and calculate the sum of the
/// ordinary leaves. This algorithm is described in section
/// 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March
/// 2015.
///
template <int MU, typename T, typename P>
T Phi0_thread(T x,
              int64_t z,
              uint64_t b,
              int64_t k,
              T square_free,
              vector<P>& primes)
{
  T phi0 = 0;

  for (b++; b < primes.size(); b++)
  {
    T next = square_free * primes[b];
    if (next > z) break;
    phi0 += MU * phi_tiny(x / next, k);
    phi0 += Phi0_thread<-MU>(x, z, b, k, next, primes);
  }

  return phi0;
}

/// Parallel computation of the ordinary leaves.
/// Run time: O(z)
/// Memory usage: O(y / log(y))
///
template <typename X, typename Y>
X Phi0_OpenMP(X x,
              Y y,
              int64_t z,
              int64_t k,
              int threads)
{
  threads = ideal_num_threads(threads, y);
  int proc_id = mpi_proc_id();
  int procs = mpi_num_procs();
  auto primes = generate_primes<Y>(y);
  int64_t pi_y = primes.size() - 1;
  X phi0 = 0;

  if (proc_id == 0)
    phi0 += phi_tiny(x, k);

  #pragma omp parallel for schedule(static, 1) num_threads(threads) reduction (+: phi0)
  for (int64_t b = k + 1 + proc_id; b <= pi_y; b += procs)
  {
    phi0 -= phi_tiny(x / primes[b], k);
    phi0 += Phi0_thread<1>(x, z, b, k, (X) primes[b], primes);
  }

  phi0 = mpi_reduce_sum(phi0);

  return phi0;
}

} // namespace

namespace primecount {

int64_t Phi0_mpi(int64_t x,
                 int64_t y,
                 int64_t z,
                 int64_t k,
                 int threads)
{
  print("");
  print("=== Phi0_mpi(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int64_t phi0 = Phi0_OpenMP(x, y, z, k, threads);

  print("Phi0", phi0, time);
  return phi0;
}

#ifdef HAVE_INT128_T

int128_t Phi0_mpi(int128_t x,
                  int64_t y,
                  int64_t z,
                  int64_t k,
                  int threads)
{
  print("");
  print("=== Phi0_mpi(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int128_t phi0;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
    phi0 = Phi0_OpenMP(x, (uint32_t) y, z, k, threads);
  else
    phi0 = Phi0_OpenMP(x, y, z, k, threads);

  print("Phi0", phi0, time);
  return phi0;
}

#endif

} // namespace

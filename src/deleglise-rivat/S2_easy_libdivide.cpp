///
/// @file  S2_easy_libdivide.cpp
/// @brief This is an optimized version of S2_easy which uses
///        libdivide. libdivide allows to replace expensive integer
///        divides with comparatively cheap multiplication and
///        bitshifts.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <S2Status.hpp>
#include <S2.hpp>
#include <json.hpp>

#include <libdivide.h>
#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

template <typename T>
void backup(T x,
            int64_t y,
            int64_t z,
            int64_t c,
            int64_t start,
            int64_t pi_x13,
            T s2_easy,
            double percent,
            double time)
{
  auto j = load_backup();

  j["S2_easy"]["x"] = to_string(x);
  j["S2_easy"]["y"] = y;
  j["S2_easy"]["z"] = z;
  j["S2_easy"]["c"] = c;
  j["S2_easy"]["start"] = start;
  j["S2_easy"]["pi_x13"] = pi_x13;
  j["S2_easy"]["s2_easy"] = to_string(s2_easy);
  j["S2_easy"]["percent"] = percent;
  j["S2_easy"]["seconds"] = get_wtime() - time;

  store_backup(j);
}

template <typename T>
void print_resume(T x,
                  int64_t start,
                  int64_t pi_x13,
                  T s2_easy,
                  double seconds,
                  double percent)
{
  if (!print_variables())
    print_log("");

  print_log("=== Resuming from " + backup_file() + " ===");
  print_log("start", start);
  print_log("pi_x13", pi_x13);
  print_log("s2_easy", s2_easy);
  print_log_seconds(seconds);
  print_status(percent, x);
}

template <typename T>
bool resume(T x,
            int64_t y,
            int64_t z,
            int64_t c,
            int64_t& start,
            int64_t& pi_x13,
            T& s2_easy,
            double& time)
{
  auto j = load_backup();

  if (is_resume(j, "S2_easy", x, y, z))
  {
    double percent = j["S2_easy"]["percent"];
    double seconds = j["S2_easy"]["seconds"];

    start = j["S2_easy"]["start"];
    pi_x13 = j["S2_easy"]["pi_x13"];
    s2_easy = calculator::eval<T>(j["S2_easy"]["s2_easy"]);
    time = get_wtime() - seconds;
    print_resume(x, start, pi_x13, s2_easy, seconds, percent);

    return true;
  }

  return false;
}

using fastdiv_t = libdivide::divider<uint64_t, libdivide::BRANCHFREE>;

vector<fastdiv_t> libdivide_primes(int64_t pi_y)
{
  vector<fastdiv_t> primes;
  primes.reserve(pi_y);
  primes.emplace_back(0);

  primesieve::iterator it;

  for (int64_t i = 0; i <= pi_y; i++)
  {
    uint64_t prime = it.next_prime();
    primes.emplace_back(prime);
  }

  return primes;
}

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
template <typename T>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 int threads,
                 double& time)
{
  T s2_easy = 0;
  int64_t start;
  int64_t pi_x13;
  double backup_time = get_wtime();
  bool is_resume = resume(x, y, z, c, start, pi_x13, s2_easy, time);

  if (is_resume && start >= pi_x13)
    return s2_easy;

  PiTable pi(y);
  S2Status status(x);
  auto primes = libdivide_primes(pi[y]);

  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t x13 = iroot<3>(x);
  int64_t max_dist = 1;
  threads = ideal_num_threads(threads, x13, 1000);

  if (!is_resume)
  {
    start = max(c, pi_sqrty) + 1;
    pi_x13 = pi[x13];
  }

  while (start <= pi_x13)
  {
    int64_t stop = min(start + max_dist * threads, pi_x13);

    #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(+: s2_easy)
    for (int64_t b = start; b <= stop; b++)
    {
      int64_t prime = primes[b].recover_divisor();
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
        int64_t xn = (uint64_t) x2 / primes[l];
        int64_t phi_xn = pi[xn] - b + 2;
        int64_t xm = (uint64_t) x2 / primes[b + phi_xn - 1];
        int64_t l2 = pi[xm];
        s2_easy += phi_xn * (l - l2);
        l = l2;
      }

      // Find all sparse easy leaves:
      // n = primes[b] * primes[l]
      // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
      for (; l > pi_min_sparse; l--)
      {
        int64_t xn = (uint64_t) x2 / primes[l];
        s2_easy += pi[xn] - b + 2;
      }

      if (is_print())
        status.print(b, pi_x13);
    }

    start = stop + 1;

    if (get_wtime() - backup_time < 300)
      max_dist *= 2;
    else
    {
      max_dist = ceil_div(max_dist, 2);
      double percent = status.getPercent(start, pi_x13, start, pi_x13);
      backup(x, y, z, c, start, pi_x13, s2_easy, percent, time);
      backup_time = get_wtime();
    }
  }

  backup(x, y, z, c, pi_x13, pi_x13, s2_easy, 100, time);

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
  print(x, y, c, threads);

  double time = get_wtime();
  int64_t s2_easy = S2_easy_OpenMP((intfast64_t) x, y, z, c, threads, time);

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
  print(x, y, c, threads);

  double time = get_wtime();
  int128_t s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, threads, time);

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

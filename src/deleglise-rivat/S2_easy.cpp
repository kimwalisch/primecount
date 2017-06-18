///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <calculator.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <S2Status.hpp>
#include <S2.hpp>
#include <json.hpp>

#include <stdint.h>
#include <vector>
#include <iostream>
#include <iomanip>

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
  if (is_print())
  {
    if (!print_variables())
      cout << endl;

    cout << "=== Resuming from primecount.backup ===" << endl;
    cout << "start = " << start << endl;
    cout << "pi_x13 = " << pi_x13 << endl;
    cout << "s2_easy = " << s2_easy << endl;
    cout << "Seconds: " << seconds << endl << endl;
    cout << "Status: " << fixed << setprecision(get_status_precision(x)) << percent << '%' << flush;
  }
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

/// Calculate the contribution of the clustered easy leaves
/// and the sparse easy leaves.
/// @param T  either int64_t or uint128_t.
///
template <typename T, typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 Primes& primes,
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

  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t x13 = iroot<3>(x);
  int64_t max_dist = 1;
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);

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
        int64_t xn = (int64_t) fast_div(x2, primes[l]);
        int64_t phi_xn = pi[xn] - b + 2;
        int64_t xm = (int64_t) fast_div(x2, primes[b + phi_xn - 1]);
        int64_t l2 = pi[xm];
        s2_easy += phi_xn * (l - l2);
        l = l2;
      }

      // Find all sparse easy leaves:
      // n = primes[b] * primes[l]
      // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
      for (; l > pi_min_sparse; l--)
      {
        int64_t xn = (int64_t) fast_div(x2, primes[l]);
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
  auto primes = generate_primes<int32_t>(y);
  int64_t s2_easy = S2_easy_OpenMP((intfast64_t) x, y, z, c, primes, threads, time);

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

  int128_t s2_easy;
  double time = get_wtime();

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(y);
    s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, primes, threads, time);
  }
  else
  {
    auto primes = generate_primes<int64_t>(y);
    s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, primes, threads, time);
  }

  print("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

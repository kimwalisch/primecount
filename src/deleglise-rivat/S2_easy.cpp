///
/// @file  S2_easy.cpp
/// @brief Calculate the contribution of the clustered easy leaves
///        and the sparse easy leaves in parallel using OpenMP
///        (Deleglise-Rivat algorithm).
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
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
#include <string>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {

/// backup to file every 60 seconds
bool is_backup(double time)
{
  return get_wtime() - time > 60;
}

/// backup to file
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            int threads,
            double percent,
            double time)
{
  json["S2_easy"]["x"] = to_string(x);
  json["S2_easy"]["y"] = y;
  json["S2_easy"]["z"] = z;
  json["S2_easy"]["threads"] = threads;
  json["S2_easy"]["percent"] = percent;
  json["S2_easy"]["seconds"] = get_wtime() - time;

  store_backup(json);
}

/// backup thread
template <typename T, typename J>
void backup(J& json,
            int64_t start,
            int64_t b,
            int64_t thread_id,
            T s2_easy)
{
  string tid = "thread" + to_string(thread_id);

  json["S2_easy"]["start"] = start;
  json["S2_easy"][tid]["b"] = b;
  json["S2_easy"][tid]["s2_easy"] = to_string(s2_easy);
}

/// backup result
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            T s2_easy,
            double time)
{
  json.erase("S2_easy");

  json["S2_easy"]["x"] = to_string(x);
  json["S2_easy"]["y"] = y;
  json["S2_easy"]["z"] = z;
  json["S2_easy"]["s2_easy"] = to_string(s2_easy);
  json["S2_easy"]["percent"] = 100;
  json["S2_easy"]["seconds"] = get_wtime() - time;

  store_backup(json);
}

template <typename T, typename J>
int get_start(J& json,
              T x,
              int64_t y,
              int64_t z,
              int64_t start)
{
  if (is_resume(json, "S2_easy", x, y, z) &&
      json["S2_easy"].count("start"))
  {
    start = json["S2_easy"]["start"];
  }

  return start;
}

template <typename T, typename J>
int get_threads(J& json,
                T x,
                int64_t y,
                int64_t z,
                int threads)
{
  if (is_resume(json, "S2_easy", x, y, z) &&
      json["S2_easy"].count("threads"))
  {
    int threads2 = json["S2_easy"]["threads"];
    threads = max(threads, threads2);
    if (!print_variables())
      print_log("");
  }

  return threads;
}

template <typename T, typename J>
double get_time(J& json,
                T x,
                int64_t y,
                int64_t z,
                double time)
{
  if (is_resume(json, "S2_easy", x, y, z))
  {
    double seconds = json["S2_easy"]["seconds"];
    time = get_wtime() - seconds;
  }

  return time;
}

template <typename T>
void print_resume(T x,
                  T s2_easy,
                  double seconds)
{
  if (!print_variables())
    print_log("");

  print_log("=== Resuming from " + backup_file() + " ===");
  print_log("s2_easy", s2_easy);
  print_log_seconds(seconds);
}

template <typename T>
void print_resume(int thread_id,
                  int64_t b,
                  T s2_easy)
{
  print_log("\r=== Resuming from " + backup_file() + " ===");
  print_log("thread", thread_id);
  print_log("b", b);
  print_log("s2_easy", s2_easy);
  print_log("");
}

/// resume result
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            T& s2_easy,
            double& time)
{
  if (is_resume(json, "S2_easy", x, y, z) &&
      json["S2_easy"].count("s2_easy"))
  {
    double seconds = json["S2_easy"]["seconds"];
    s2_easy = calculator::eval<T>(json["S2_easy"]["s2_easy"]);
    time = get_wtime() - seconds;
    print_resume(x, s2_easy, seconds);
    return true;
  }

  return false;
}

/// resume thread
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t& b,
            int thread_id,
            T& s2_easy)
{
  bool resumed = false;

  #pragma omp critical (s2_easy_resume)
  if (is_resume(json, "S2_easy", thread_id, x, y, z))
  {
    resumed = true;
    string tid = "thread" + to_string(thread_id);
    b = json["S2_easy"][tid]["b"];
    s2_easy = calculator::eval<T>(json["S2_easy"][tid]["s2_easy"]);
    print_resume(thread_id, b, s2_easy);
  }

  return resumed;
}

template <typename T, typename Primes>
T S2_easy(T x,
          int64_t y,
          int64_t z,
          int64_t b,
          Primes& primes,
          PiTable& pi)
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
  T s2_easy = 0;

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

  return s2_easy;
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
  T s2 = 0;
  auto json = load_backup();
  auto copy = json;

  if (resume(json, x, y, z, s2, time))
    return s2;

  if (!is_resume(json, "S2_easy", x, y, z))
    json.erase("S2_easy");

  PiTable pi(y);
  S2Status status(x);
  double backup_time = get_wtime();
  time = get_time(json, x, y, z, time);

  int64_t x13 = iroot<3>(x);
  int64_t pi_x13 = pi[x13];
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t start = max(c, pi_sqrty) + 1;
  start = get_start(json, x, y, z, start);

  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);
  threads = get_threads(json, x, y, z, threads);

  #pragma omp parallel for reduction(+: s2)
  for (int i = 0; i < threads; i++)
  {
    T s2_easy = 0;
    int64_t b = 0;

    if (resume(copy, x, y, z, b, i, s2_easy))
      s2_easy += S2_easy(x, y, z, b, primes, pi);

    while (true)
    {
      if (is_print())
        status.print(b, pi_x13);

      #pragma omp critical (s2_easy_backup)
      {
        b = start++;

        backup(json, start, b, i, s2_easy);

        if (is_backup(backup_time))
        {
          double percent = status.getPercent(start, pi_x13, start, pi_x13);
          backup(json, x, y, z, threads, percent, time);
          backup_time = get_wtime();
        }
      }

      if (b > pi_x13)
        break;

      s2_easy += S2_easy(x, y, z, b, primes, pi);
    }

    s2 += s2_easy;
  }

  backup(json, x, y, z, s2, time);

  return s2;
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

  print_log("");
  print_log("=== S2_easy(x, y) ===");
  print_log("Computation of the easy special leaves");
  print_log(x, y, c, threads);

  double time = get_wtime();
  auto primes = generate_primes<int32_t>(y);
  int64_t s2_easy = S2_easy_OpenMP((intfast64_t) x, y, z, c, primes, threads, time);

  print_log("S2_easy", s2_easy, time);
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

  print_log("");
  print_log("=== S2_easy(x, y) ===");
  print_log("Computation of the easy special leaves");
  print_log(x, y, c, threads);

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

  print_log("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

///
/// @file  P2.cpp
/// @brief 2nd partial sieve function. P2(x, y) counts the
///        numbers <= x that have exactly 2 prime factors
///        each exceeding the a-th prime.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <calculator.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <json.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;
using namespace nlohmann;

namespace {

/// Count the primes inside [prime, stop]
template <typename T>
int64_t count_primes(primesieve::iterator& it, int64_t& prime, T stop)
{
  int64_t count = 0;

  for (; prime <= stop; count++)
    prime = it.next_prime();

  return count;
}

/// Calculate the thread sieving distance. The idea is to
/// gradually increase the thread_distance in order to
/// keep all CPU cores busy.
///
void balanceLoad(int64_t* thread_distance, 
                 int64_t low,
                 int64_t z,
                 int threads,
                 double start_time)
{
  double seconds = get_wtime() - start_time;

  int64_t min_distance = 1 << 23;
  int64_t max_distance = ceil_div(z - low, threads);

  if (seconds < 60)
    *thread_distance *= 2;
  if (seconds > 60)
    *thread_distance /= 2;

  *thread_distance = in_between(min_distance, *thread_distance, max_distance);
}

template <typename T>
T P2_thread(T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t thread_num,
            int64_t thread_distance,
            int64_t& pix,
            int64_t& pix_count)
{
  pix = 0;
  pix_count = 0;
  low += thread_distance * thread_num;
  z = min(low + thread_distance, z);
  int64_t start = (int64_t) max(x / z, y);
  int64_t stop = (int64_t) min(x / low, isqrt(x));

  primesieve::iterator rit(stop + 1, start);
  primesieve::iterator it(low - 1, z);

  int64_t next = it.next_prime();
  int64_t prime = rit.prev_prime();
  T p2 = 0;

  // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i])
  while (prime > start &&
         x / prime < z)
  {
    pix += count_primes(it, next, x / prime);
    p2 += pix;
    pix_count++;
    prime = rit.prev_prime();
  }

  pix += count_primes(it, next, z - 1);

  return p2;
}

template <typename T>
void backup(T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t thread_distance,
            T pix_total,
            T p2,
            double time)
{
  ifstream ifs("primecount.backup");
  json j;

  if (ifs.is_open())
  {
    ifs >> j;
    ifs.close();
  }

  double percent = get_percent(low, z);

  j["P2"]["x"] = to_string(x);
  j["P2"]["y"] = y;
  j["P2"]["z"] = z;
  j["P2"]["low"] = low;
  j["P2"]["thread_distance"] = thread_distance;
  j["P2"]["pix_total"] = to_string(pix_total);
  j["P2"]["p2"] = to_string(p2);
  j["P2"]["percent"] = percent;
  j["P2"]["seconds"] = get_wtime() - time;

  ofstream ofs("primecount.backup");
  ofs << setw(4) << j << endl;
}

template <typename T>
void resume(T x,
            int64_t y,
            int64_t z,
            int64_t& low,
            int64_t& thread_distance,
            T& pix_total,
            T& p2,
            double& time)
{
  ifstream ifs("primecount.backup");
  json j;

  if (ifs.is_open())
  {
    ifs >> j;
    ifs.close();
  }

  if (j.find("P2") != j.end() &&
      x == calculator::eval<T>(j["P2"]["x"]) &&
      y == j["P2"]["y"] &&
      z == j["P2"]["z"])
  {
    double percent = j["P2"]["percent"];
    double seconds = j["P2"]["seconds"];

    low = j["P2"]["low"];
    thread_distance = j["P2"]["thread_distance"];
    pix_total = calculator::eval<T>(j["P2"]["pix_total"]);
    p2 = calculator::eval<T>(j["P2"]["p2"]);
    time = get_wtime() - seconds;

    if (is_print())
    {
      if (!print_variables())
        cout << endl;

      cout << "=== Resuming from primecount.backup ===" << endl;
      cout << "low = " << low << endl;
      cout << "thread_distance = " << thread_distance << endl;
      cout << "pix_total = " << pix_total << endl;
      cout << "p2 = " << p2 << endl;
      cout << "Seconds: " << fixed << setprecision(3) << seconds << endl << endl;
      cout << "Status: " << fixed << setprecision(get_status_precision(x)) << percent << '%' << flush;
    }
  }
}

/// P2(x, y) counts the numbers <= x that have exactly 2
/// prime factors each exceeding the a-th prime.
/// Run-time: O(z log log z)
///
template <typename T>
T P2_OpenMP(T x, int64_t y, int threads, double& time)
{
  static_assert(prt::is_signed<T>::value,
                "P2(T x, ...): T must be signed integer type");

  if (x < 4)
    return 0;

  T a = pi_legendre(y, threads);
  T b = pi_legendre((int64_t) isqrt(x), threads);

  if (a >= b)
    return 0;

  // \sum_{i=a+1}^{b} -(i - 1)
  T p2 = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;
  T pix_total = 0;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, 1));
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;

  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  resume(x, y, z, low, thread_distance, pix_total, p2, time);

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double t = get_wtime();

    #pragma omp parallel for num_threads(threads) reduction(+: p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_thread(x, y, z, low, i, thread_distance, pix[i], pix_counts[i]);

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, t);

    // add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      p2 += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    backup(x, y, z, low, thread_distance, pix_total, p2, time);

    if (is_print())
    {
      double percent = get_percent(low, z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  return p2;
}

} // namespace

namespace primecount {

int64_t P2(int64_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int64_t p2 = P2_OpenMP(x, y, threads, time);

  print("P2", p2, time);
  return p2;
}

#ifdef HAVE_INT128_T

int128_t P2(int128_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print("");
  print("=== P2(x, y) ===");
  print("Computation of the 2nd partial sieve function");
  print(x, y, threads);

  double time = get_wtime();
  int128_t p2 = P2_OpenMP(x, y, threads, time);

  print("P2", p2, time);
  return p2;
}

#endif

} // namespace

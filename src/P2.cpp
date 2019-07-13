///
/// @file  P2.cpp
/// @brief P2(x, a) is the 2nd partial sieve function.
///        P2(x, a) counts the numbers <= x that have exactly 2 prime
///        factors each exceeding the a-th prime. This implementation
///        uses the primesieve library for quickly iterating over
///        primes using next_prime() and prev_prime() which greatly
///        simplifies the implementation.
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

#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <calculator.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <json.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace primecount;

namespace {

/// backup to file every 60 seconds
bool is_backup(double time)
{
  double seconds = get_time() - time;
  return seconds > 60;
}

/// Count the primes inside [prime, stop]
int64_t count_primes(primesieve::iterator& it, int64_t& prime, int64_t stop)
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
  double seconds = get_time() - start_time;

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
  T p2 = 0;
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

  // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i])
  while (prime > start)
  {
    int64_t xp = (int64_t)(x / prime);
    if (xp >= z) break;
    pix += count_primes(it, next, xp);
    pix_count++;
    p2 += pix;
    prime = rit.prev_prime();
  }

  pix += count_primes(it, next, z - 1);

  return p2;
}

// backup intermediate result
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t thread_distance,
            T pix_total,
            T p2,
            double time)
{
  double percent = get_percent(low, z);

  json["P2"]["x"] = to_string(x);
  json["P2"]["y"] = y;
  json["P2"]["z"] = z;
  json["P2"]["low"] = low;
  json["P2"]["thread_distance"] = thread_distance;
  json["P2"]["pix_total"] = to_string(pix_total);
  json["P2"]["p2"] = to_string(p2);
  json["P2"]["percent"] = percent;
  json["P2"]["seconds"] = get_time() - time;

  store_backup(json);
}

// backup result
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            T p2,
            double time)
{
  if (json.find("P2") != json.end())
    json.erase("P2");

  json["P2"]["x"] = to_string(x);
  json["P2"]["y"] = y;
  json["P2"]["z"] = z;
  json["P2"]["p2"] = to_string(p2);
  json["P2"]["percent"] = 100.0;
  json["P2"]["seconds"] = get_time() - time;

  store_backup(json);
}

// resume computation
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t& low,
            int64_t& thread_distance,
            T& pix_total,
            T& p2,
            double& time)
{
  if (is_resume(json, "P2", x, y, z))
  {
    double seconds = json["P2"]["seconds"];
    low = json["P2"]["low"];
    thread_distance = json["P2"]["thread_distance"];
    pix_total = calculator::eval<T>(json["P2"]["pix_total"]);
    p2 = calculator::eval<T>(json["P2"]["p2"]);
    time = get_time() - seconds;
    return true;
  }

  return false;
}

/// resume result
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            T& p2,
            double& time)
{
  if (is_resume(json, "P2", x, y, z))
  {
    double percent = json["P2"]["percent"];
    double seconds = json["P2"]["seconds"];
    print_resume(percent, x);

    if (!json["P2"].count("low"))
    {
      p2 = calculator::eval<T>(json["P2"]["p2"]);
      time = get_time() - seconds;
      return true;
    }
  }

  return false;
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

  T p2 = 0;
  int64_t z = (int64_t)(x / max(y, 1));
  auto json = load_backup();
  if (resume(json, x, y, z, p2, time))
    return p2;

  double backup_time = get_time();
  T a = pi_legendre(y, threads);
  T b = pi_legendre((int64_t) isqrt(x), threads);

  if (a >= b)
    return 0;

  // \sum_{i=a+1}^{b} -(i - 1)
  p2 = (a - 2) * (a + 1) / 2 - (b - 2) * (b + 1) / 2;
  T pix_total = 0;

  int64_t low = 2;
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;

  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  if (!resume(json, x, y, z, low, thread_distance, pix_total, p2, time))
    if (json.find("P2") != json.end())
      json.erase("P2");

  // \sum_{i=a+1}^{b} pi(x / primes[i])
  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double t = get_time();

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

    if (is_backup(backup_time))
    {
      backup(json, x, y, z, low, thread_distance, pix_total, p2, time);
      backup_time = get_time();
    }

    if (is_print())
    {
      double percent = get_percent(low, z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  backup(json, x, y, z, p2, time);

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

  print_log("");
  print_log("=== P2(x, y) ===");
  print_log("Computation of the 2nd partial sieve function");
  print_log(x, y, threads);

  double time = get_time();
  int64_t p2 = P2_OpenMP(x, y, threads, time);

  print_log("P2", p2, time);
  return p2;
}

#ifdef HAVE_INT128_T

int128_t P2(int128_t x, int64_t y, int threads)
{
#ifdef HAVE_MPI
  if (mpi_num_procs() > 1)
    return P2_mpi(x, y, threads);
#endif

  print_log("");
  print_log("=== P2(x, y) ===");
  print_log("Computation of the 2nd partial sieve function");
  print_log(x, y, threads);

  double time = get_time();
  int128_t p2 = P2_OpenMP(x, y, threads, time);

  print_log("P2", p2, time);
  return p2;
}

#endif

} // namespace

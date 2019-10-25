///
/// @file  B.cpp
/// @brief The B formula is a partial computation of the P2(x, a)
///        formula from the Lagarias-Miller-Odlyzko and Deleglise-
///        Rivat prime counting algorithms. P2(x, a) counts the
///        numbers <= x that have exactly 2 prime factors each
///        exceeding the a-th prime. Both P2 and B have a runtime
///        complexity of O(z log log z) and use O(z^(1/2)) memory,
///        with z = x / y. This implementation is a simplified
///        version of P2.cpp.
///
///        B(x, y) formula:
///        \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <calculator.hpp>
#include <aligned_vector.hpp>
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

/// backup intermediate result
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t thread_distance,
            T pix_total,
            T b,
            double time)
{
  double percent = get_percent(low, z);

  json["B"]["x"] = to_string(x);
  json["B"]["y"] = y;
  json["B"]["low"] = low;
  json["B"]["thread_distance"] = thread_distance;
  json["B"]["sieve_limit"] = z;
  json["B"]["pix_total"] = to_string(pix_total);
  json["B"]["b"] = to_string(b);
  json["B"]["percent"] = percent;
  json["B"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// backup result
template <typename T>
void backup(T x,
            int64_t y,
            int64_t z,
            T b,
            double time)
{
  auto json = load_backup();

  if (json.find("B") != json.end())
    json.erase("B");

  json["B"]["x"] = to_string(x);
  json["B"]["y"] = y;
  json["B"]["b"] = to_string(b);
  json["B"]["sieve_limit"] = z;
  json["B"]["percent"] = 100.0;
  json["B"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// resume computation
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t& low,
            int64_t& thread_distance,
            T& pix_total,
            T& b,
            double& time)
{
  if (is_resume(json, "B", x, y))
  {
    double seconds = json["B"]["seconds"];
    low = json["B"]["low"];
    thread_distance = json["B"]["thread_distance"];
    pix_total = calculator::eval<T>(json["B"]["pix_total"]);
    b = calculator::eval<T>(json["B"]["b"]);
    time = get_time() - seconds;
    return true;
  }

  return false;
}

/// resume result
template <typename T>
bool resume(T x,
            int64_t y,
            T& b,
            double& time)
{
  auto json = load_backup();

  if (is_resume(json, "B", x, y))
  {
    double percent = json["B"]["percent"];
    double seconds = json["B"]["seconds"];
    print_resume(percent, x);

    if (!json["B"].count("low"))
    {
      b = calculator::eval<T>(json["B"]["b"]);
      time = get_time() - seconds;
      return true;
    }
  }

  return false;
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
T B_thread(T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t thread_num,
            int64_t thread_distance,
            int64_t& pix,
            int64_t& pix_count)
{
  T sum = 0;
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
    sum += pix;
    prime = rit.prev_prime();
  }

  pix += count_primes(it, next, z - 1);

  return sum;
}

/// \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
/// Run time: O(z log log z)
/// Memory usage: O(z^(1/2))
///
template <typename T>
T B_OpenMP(T x, 
           int64_t y,
           int64_t z,
           int threads,
           double& time)
{
  if (x < 4)
    return 0;

  T sum = 0;
  T pix_total = 0;

  int64_t low = 2;
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;

  aligned_vector<int64_t> pix(threads);
  aligned_vector<int64_t> pix_counts(threads);

  auto json = load_backup();
  if (!resume(json, x, y, low, thread_distance, pix_total, sum, time))
    if (json.find("B") != json.end())
      json.erase("B");

  double last_backup_time = get_time();

  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double t = get_time();

    #pragma omp parallel for num_threads(threads) reduction(+: sum)
    for (int i = 0; i < threads; i++)
      sum += B_thread(x, y, z, low, i, thread_distance, pix[i], pix_counts[i]);

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, t);

    // add missing sum contributions in order
    for (int i = 0; i < threads; i++)
    {
      sum += pix_total * pix_counts[i];
      pix_total += pix[i];
    }

    if (is_backup(last_backup_time))
    {
      backup(json, x, y, z, low, thread_distance, pix_total, sum, time);
      last_backup_time = get_time();
    }

    if (is_print())
    {
      double percent = get_percent(low, z);
      cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
           << percent << '%' << flush;
    }
  }

  return sum;
}

} // namespace

namespace primecount {

int64_t B(int64_t x, int64_t y, int threads)
{
  print("");
  print("=== B(x, y) ===");
  print_gourdon(x, y, threads);

  double time = get_time();
  int64_t b = 0;

  if (!resume(x, y, b, time))
  {
    int64_t z = (int64_t)(x / max(y, 1));
    b = B_OpenMP((intfast64_t) x, y, z, threads, time);
    backup(x, y, z, b, time);
  }

  print("B", b, time);
  return b;
}

#ifdef HAVE_INT128_T

int128_t B(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== B(x, y) ===");
  print_gourdon(x, y, threads);

  double time = get_time();
  int128_t b = 0;

  if (!resume(x, y, b, time))
  {
    int64_t z = (int64_t)(x / max(y, 1));
    b = B_OpenMP((intfast128_t) x, y, z, threads, time);
    backup(x, y, z, b, time);
  }

  print("B", b, time);
  return b;
}

#endif

} // namespace

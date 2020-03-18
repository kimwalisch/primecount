///
/// @file  B.cpp
/// @brief The B formula is a partial computation of the P2(x, a)
///        formula from the Lagarias-Miller-Odlyzko and Deleglise-Rivat
///        prime counting algorithms. P2(x, a) counts the numbers <= x
///        that have exactly 2 prime factors each exceeding the a-th
///        prime. Both P2 and B have a runtime complexity of
///        O(z log log z) and use O(z^(1/2)) memory, with z = x / y.
///
///        B(x, y) formula:
///        \sum_{i=pi[y]+1}^{pi[x^(1/2)]} pi(x / primes[i])
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <gourdon.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <aligned_vector.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <noinline.hpp>
#include <imath.hpp>
#include <json.hpp>
#include <backup.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <type_traits>

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
template <typename T>
void backup(nlohmann::json& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t low,
            int64_t pi_low_minus_1,
            int64_t thread_distance,
            T sum,
            double time)
{
  using ST = typename make_signed<T>::type;
  double percent = get_percent(low, z);

  auto& B = json["B"];
  B["x"] = to_str(x);
  B["y"] = y;
  B["alpha_y"] = get_alpha_y(x, y);
  B["low"] = low;
  B["pi_low_minus_1"] = pi_low_minus_1;
  B["thread_distance"] = thread_distance;
  B["sieve_limit"] = z;
  B["sum"] = to_str((ST) sum);
  B["percent"] = percent;
  B["seconds"] = get_time() - time;

  store_backup(json);
}

/// backup result
template <typename T>
void backup(T x,
            int64_t y,
            int64_t z,
            T sum,
            double time)
{
  using ST = typename make_signed<T>::type;
  auto json = load_backup();

  if (json.find("B") != json.end())
    json.erase("B");

  auto& B = json["B"];
  B["x"] = to_str(x);
  B["y"] = y;
  B["alpha_y"] = get_alpha_y(x, y);
  B["sum"] = to_str((ST) sum);
  B["sieve_limit"] = z;
  B["percent"] = 100.0;
  B["seconds"] = get_time() - time;

  store_backup(json);
}

/// resume computation
template <typename T>
bool resume(nlohmann::json& json,
            T x,
            int64_t y,
            int64_t& low,
            int64_t& pi_low_minus_1,
            int64_t& thread_distance,
            T& sum,
            double& time)
{
  if (is_resume(json, "B", x, y))
  {
    double seconds = json["B"]["seconds"];
    low = json["B"]["low"];
    pi_low_minus_1 = json["B"]["pi_low_minus_1"];
    thread_distance = json["B"]["thread_distance"];
    sum = (T) to_maxint(json["B"]["sum"]);
    time = get_time() - seconds;
    return true;
  }

  return false;
}

/// resume result
template <typename T>
bool resume(T x,
            int64_t y,
            T& sum,
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
      sum = (T) to_maxint(json["B"]["sum"]);
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
NOINLINE std::tuple<T, int64_t, int64_t>
B_thread(T x,
         int64_t y,
         int64_t z,
         int64_t low,
         int64_t thread_num,
         int64_t thread_distance)
{
  T sum = 0;
  int64_t pix = 0;
  int64_t pix_count = 0;

  // thread sieves [low, z[
  low += thread_distance * thread_num;
  z = min(low + thread_distance, z);
  int64_t start = (int64_t) max(x / z, y);
  int64_t stop = (int64_t) min(x / low, isqrt(x));

  primesieve::iterator rit(stop + 1, start);
  primesieve::iterator it(low - 1, z);

  int64_t next = it.next_prime();
  int64_t prime = rit.prev_prime();

  // \sum_{i = pi[start]+1}^{pi[stop]} pi(x / primes[i]) - pi(low - 1)
  while (prime > start)
  {
    int64_t xp = (int64_t)(x / prime);
    if (xp >= z) break;
    pix += count_primes(it, next, xp);
    pix_count++;
    sum += pix;
    prime = rit.prev_prime();
  }

  // prime count [low, z[
  pix += count_primes(it, next, z - 1);
  auto res = make_tuple(sum, pix, pix_count);

  return res;
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
  int64_t low = 2;
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;
  int64_t pi_low_minus_1 = 0;

  // prevents CPU false sharing
  using res_t = std::tuple<T, int64_t, int64_t>;
  aligned_vector<res_t> res(threads);

  auto json = load_backup();
  if (!resume(json, x, y, low, pi_low_minus_1, thread_distance, sum, time))
    if (json.find("B") != json.end())
      json.erase("B");

  double last_backup_time = get_time();

  while (low < z)
  {
    int64_t max_threads = ceil_div(z - low, thread_distance);
    threads = in_between(1, threads, max_threads);
    double t = get_time();

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++)
      res[i] = B_thread(x, y, z, low, i, thread_distance);

    // The threads above have computed the sum of:
    // PrimePi(n) - PrimePi(thread_low - 1)
    // for many different values of n. However we actually want to
    // compute the sum of PrimePi(n). In order to get the complete
    // sum we now have to calculate the missing sum contributions in
    // sequential order as each thread depends on values from the
    // previous thread. The missing sum contribution for each thread
    // can be calculated using pi_low_minus_1 * thread_count.
    for (int i = 0; i < threads; i++)
    {
      auto thread_sum = std::get<0>(res[i]);
      auto thread_pix = std::get<1>(res[i]);
      auto thread_count = std::get<2>(res[i]);
      thread_sum += (T) pi_low_minus_1 * thread_count;
      sum += thread_sum;
      pi_low_minus_1 += thread_pix;
    }

    low += thread_distance * threads;
    balanceLoad(&thread_distance, low, z, threads, t);

    if (is_backup(last_backup_time))
    {
      backup(json, x, y, z, low, pi_low_minus_1, thread_distance, sum, time);
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
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int64_t sum = 0;

  if (!resume(x, y, sum, time))
  {
    int64_t z = (int64_t)(x / max(y, 1));
    sum = B_OpenMP((uint64_t) x, y, z, threads, time);
    backup(x, y, z, sum, time);
  }

  print("B", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t B(int128_t x, int64_t y, int threads)
{
  print("");
  print("=== B(x, y) ===");
  print_gourdon_vars(x, y, threads);

  double time = get_time();
  int128_t sum = 0;

  if (!resume(x, y, sum, time))
  {
    int64_t z = (int64_t)(x / max(y, 1));
    sum = B_OpenMP((uint128_t) x, y, z, threads, time);
    backup(x, y, z, sum, time);
  }

  print("B", sum, time);
  return sum;
}

#endif

} // namespace

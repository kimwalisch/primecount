///
/// @file  S2_easy_libdivide.cpp
/// @brief This is an optimized version of S2_easy which uses
///        libdivide. libdivide allows to replace expensive integer
///        divides with comparatively cheap multiplication and
///        bitshifts.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <S2Status.hpp>
#include <S2.hpp>
#include <json.hpp>

#include <libdivide.h>
#include <stdint.h>
#include <vector>
#include <string>

using namespace std;
using namespace primecount;

namespace {

/// backup to file every 60 seconds
bool is_backup(double time)
{
  double seconds = get_time() - time;
  return seconds > 60;
}

/// backup to file
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            double percent,
            double time)
{
  json["S2_easy"]["x"] = to_string(x);
  json["S2_easy"]["y"] = y;
  json["S2_easy"]["z"] = z;
  json["S2_easy"]["percent"] = percent;
  json["S2_easy"]["seconds"] = get_time() - time;

  store_backup(json);
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
  json["S2_easy"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// update json
template <typename T, typename J>
void update(J& json,
            int64_t start,
            int64_t b,
            int64_t thread_id,
            T s2_easy)
{
  string tid = "thread" + to_string(thread_id);

  json["S2_easy"][tid]["b"] = b;
  json["S2_easy"]["start"] = start;
  json["S2_easy"]["s2_easy"] = to_string(s2_easy);
}

/// resume thread
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t& b,
            int thread_id)
{
  if (is_resume(json, "S2_easy", thread_id, x, y, z))
  {
    string tid = "thread" + to_string(thread_id);
    b = json["S2_easy"][tid]["b"];
    return true;
  }

  return false;
}

/// resume vars
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            T& s2_easy,
            int64_t& start,
            double& time)
{
  if (is_resume(json, "S2_easy", x, y, z))
  {
    double seconds = json["S2_easy"]["seconds"];
    s2_easy = calculator::eval<T>(json["S2_easy"]["s2_easy"]);
    start = json["S2_easy"]["start"];
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
            T& s2_easy,
            double& time)
{
  if (is_resume(json, "S2_easy", x, y, z))
  {
    double percent = json["S2_easy"]["percent"];
    double seconds = json["S2_easy"]["seconds"];
    print_resume(percent, x);

    if (!json["S2_easy"].count("start"))
    {
      s2_easy = calculator::eval<T>(json["S2_easy"]["s2_easy"]);
      time = get_time() - seconds;
      return true;
    }
  }

  return false;
}

template <typename T>
bool is_libdivide(T x)
{
  return x <= numeric_limits<uint64_t>::max();
}

using fastdiv_t = libdivide::branchfree_divider<uint64_t>;

template <typename Primes>
vector<fastdiv_t>
libdivide_vector(Primes& primes)
{
  vector<fastdiv_t> fastdiv(1);
  fastdiv.insert(fastdiv.end(), primes.begin() + 1, primes.end());
  return fastdiv;
}

template <typename T, typename Primes, typename Fastdiv>
T S2_easy(T x,
          int64_t y,
          int64_t z,
          int64_t b,
          Primes& primes,
          Fastdiv& fastdiv,
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

  if (is_libdivide(x2))
  {
    // Find all clustered easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (l > pi_min_clustered)
    {
      int64_t xn = (uint64_t) x2 / fastdiv[l];
      int64_t phi_xn = pi[xn] - b + 2;
      int64_t xm = (uint64_t) x2 / fastdiv[b + phi_xn - 1];
      int64_t l2 = pi[xm];
      s2_easy += phi_xn * (l - l2);
      l = l2;
    }

    // Find all sparse easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; l > pi_min_sparse; l--)
    {
      int64_t xn = (uint64_t) x2 / fastdiv[l];
      s2_easy += pi[xn] - b + 2;
    }
  }
  else
  {
    // Find all clustered easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) == phi(x / m, b - 1)
    // where phi(x / n, b - 1) = pi(x / n) - b + 2
    while (l > pi_min_clustered)
    {
      int64_t xn = (int64_t) (x2 / primes[l]);
      int64_t phi_xn = pi[xn] - b + 2;
      int64_t xm = (int64_t) (x2 / primes[b + phi_xn - 1]);
      int64_t l2 = pi[xm];
      s2_easy += phi_xn * (l - l2);
      l = l2;
    }

    // Find all sparse easy leaves:
    // n = primes[b] * primes[l]
    // x / n <= y && phi(x / n, b - 1) = pi(x / n) - b + 2
    for (; l > pi_min_sparse; l--)
    {
      int64_t xn = (int64_t) (x2 / primes[l]);
      s2_easy += pi[xn] - b + 2;
    }
  }

  return s2_easy;
}

/// Calculate the contribution of the clustered easy
/// leaves and the sparse easy leaves.
///
template <typename T, typename Primes>
T S2_easy_OpenMP(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c,
                 Primes& primes,
                 int threads,
                 double& time,
                 nlohmann::json& json)
{
  PiTable pi(y);
  S2Status status(x);
  auto fastdiv = libdivide_vector(primes);
  double backup_time = get_time();
  nlohmann::json copy = json;

  T s2_easy = 0;
  int64_t x13 = iroot<3>(x);
  int64_t pi_x13 = pi[x13];
  int64_t pi_sqrty = pi[isqrt(y)];
  int64_t start = max(c, pi_sqrty) + 1;
  threads = ideal_num_threads(threads, x13, 1000);

  if (!resume(json, x, y, z, s2_easy, start, time))
    json.erase("S2_easy");

  int resume_threads = calculate_resume_threads(json, "S2_easy");

  #pragma omp parallel for num_threads(threads)
  for (int i = 0; i < threads; i++)
  {
    int64_t b = 0;

    // 1st resume computations from backup file
    for (int j = i; j < resume_threads; j += threads)
    {
      if (resume(copy, x, y, z, b, j))
      {
        T s2_thread = S2_easy(x, y, z, b, primes, fastdiv, pi);

        #pragma omp critical (s2_easy)
        {
          s2_easy += s2_thread;
          string tid = "thread" + to_string(j);

          json["S2_easy"]["s2_easy"] = to_string(s2_easy);
          json["S2_easy"].erase(tid);
        }
      }
    }

    T s2_thread = 0;

    // 2nd, run new computations
    while (true)
    {
      if (is_print())
        status.print(b, pi_x13);

      #pragma omp critical (s2_easy)
      {
        s2_easy += s2_thread;
        b = start++;

        update(json, start, b, i, s2_easy);

        if (is_backup(backup_time))
        {
          double percent = status.getPercent(start, pi_x13, start, pi_x13);
          backup(json, x, y, z, percent, time);
          backup_time = get_time();
        }
      }

      if (b > pi_x13)
        break;

      s2_thread = S2_easy(x, y, z, b, primes, fastdiv, pi);
    }
  }

  backup(json, x, y, z, s2_easy, time);

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

  print_log("");
  print_log("=== S2_easy(x, y) ===");
  print_log("Computation of the easy special leaves");
  print_log(x, y, c, threads);

  auto json = load_backup();
  double time = get_time();
  int64_t s2_easy = 0;

  if (!resume(json, x, y, z, s2_easy, time))
  {
    auto primes = generate_primes<int32_t>(y);
    s2_easy = S2_easy_OpenMP((intfast64_t) x, y, z, c, primes, threads, time, json);
  }

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

  auto json = load_backup();
  double time = get_time();
  int128_t s2_easy = 0;

  if (!resume(json, x, y, z, s2_easy, time))
  {
    // uses less memory
    if (y <= numeric_limits<uint32_t>::max())
    {
      auto primes = generate_primes<uint32_t>(y);
      s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, primes, threads, time, json);
    }
    else
    {
      auto primes = generate_primes<int64_t>(y);
      s2_easy = S2_easy_OpenMP((intfast128_t) x, y, z, c, primes, threads, time, json);
    }
  }

  print_log("S2_easy", s2_easy, time);
  return s2_easy;
}

#endif

} // namespace

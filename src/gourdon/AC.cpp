///
/// @file  AC.cpp
/// @brief Implementation of the A + C formulas in Xavier Gourdon's
///        prime counting algorithm. In this implementation the memory
///        usage of the pi[x] lookup table has been reduced from
///        O(x^(1/2)) to O(x^(1/3)) by using a segmented pi[x] lookup
///        table. In each segment we process the leaves that satisfy:
///        low <= x / (prime * m) < high.
///
///        The A & C formulas roughly correspond to the easy special
///        leaves in the Deleglise-Rivat algorithm. Since both
///        formulas use a very similar segmented algorithm that goes
///        up to x^(1/2) it makes sense to merge the A & C formulas
///        hence reducing the runtime complexity by a factor of
///        O(x^(1/2) * ln ln x^(1/2)) and avoiding initializing some
///        data structures twice. Merging the A & C formulas also
///        improves scaling on systems with many CPU cores.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <gourdon.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <StatusAC.hpp>
#include <json.hpp>
#include <backup.hpp>

#include <stdint.h>
#include <vector>
#include <string>
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
            int64_t k,
            int64_t x_star,
            double percent,
            double time)
{
  if (json["AC"].count("percent") > 0)
  {
    double percent2 = json["AC"]["percent"];
    percent = std::max(percent, percent2);
  }

  json["AC"]["x"] = to_str(x);
  json["AC"]["y"] = y;
  json["AC"]["z"] = z;
  json["AC"]["k"] = k;
  json["AC"]["x_star"] = x_star;
  json["AC"]["alpha_y"] = get_alpha_y(x, y);
  json["AC"]["alpha_z"] = get_alpha_z(y, z);
  json["AC"]["percent"] = percent;
  json["AC"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// backup result
template <typename T>
void backup_result(nlohmann::json& json,
                   T x,
                   int64_t y,
                   int64_t z,
                   int64_t k,
                   int64_t x_star,
                   T sum,
                   double time)
{
  if (json.find("AC") != json.end())
    json.erase("AC");

  using ST = typename make_signed<T>::type;
  json["AC"]["x"] = to_str(x);
  json["AC"]["y"] = y;
  json["AC"]["z"] = z;
  json["AC"]["k"] = k;
  json["AC"]["x_star"] = x_star;
  json["AC"]["alpha_y"] = get_alpha_y(x, y);
  json["AC"]["alpha_z"] = get_alpha_z(y, z);
  json["AC"]["sum"] = to_str((ST) sum);
  json["AC"]["percent"] = 100.0;
  json["AC"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// Update from resume thread.
/// Update backup without storing to disk.
template <typename T>
void update(nlohmann::json& json,
            const std::string& tid,
            const std::string& formula,
            T sum)
{
  using ST = typename make_signed<T>::type;
  auto& AC = json["AC"];
  AC[formula] = to_str((ST) sum);
  AC.erase(tid);
}

/// Update from regular thread.
/// Update backup without storing to disk.
template <typename T>
void update(nlohmann::json& json,
            const std::string& tid,
            const std::string& formula,
            int64_t b,
            int64_t next_b,
            int64_t max_b,
            T sum)
{
  using ST = typename make_signed<T>::type;
  auto& AC = json["AC"];
  AC["next_b"] = next_b;
  AC[formula] = to_str((ST) sum);

  if (b <= max_b)
    AC[tid]["b"] = b;
  else
  {
    // finished
    if (AC.find(tid) != AC.end())
      AC.erase(tid);
  }
}

/// resume thread
template <typename T>
bool resume(nlohmann::json& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t k,
            int64_t& b,
            int thread_id)
{
  if (is_resume(json, "AC", thread_id, x, y, z, k))
  {
    string tid = "thread" + to_str(thread_id);
    b = json["AC"][tid]["b"];
    return true;
  }

  return false;
}

/// Resume C1 algorithm
template <typename T>
bool resume_c1(nlohmann::json& json,
               T x,
               int64_t y,
               int64_t z,
               int64_t k,
               int64_t& next_b,
               T& sum,
               double& time)
{
  if (is_resume(json, "AC", x, y, z, k) &&
      json["AC"].count("sum_c1") > 0)
  {
    double seconds = json["AC"]["seconds"];
    sum = (T) to_maxint(json["AC"]["sum_c1"]);
    next_b = json["AC"]["next_b"];
    time = get_time() - seconds;
    return true;
  }

  return false;
}

/// Resume the A and C2 algorithms which
/// make use of SegmentedPi.
///
template <typename T>
bool resume_a_c2(nlohmann::json& json,
                 T x,
                 int64_t y,
                 int64_t z,
                 int64_t k,
                 int64_t& next_b,
                 int64_t& low,
                 T& sum,
                 double& time)
{
  if (is_resume(json, "AC", x, y, z, k) &&
      json["AC"].count("sum_ac") > 0)
  {
    double seconds = json["AC"]["seconds"];
    sum = (T) to_maxint(json["AC"]["sum_ac"]);
    next_b = json["AC"]["next_b"];
    low = json["AC"]["low"];
    time = get_time() - seconds;
    return true;
  }

  return false;
}

/// resume result
template <typename T>
bool resume(nlohmann::json& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t k,
            T& sum,
            double& time)
{
  if (is_resume(json, "AC", x, y, z, k))
  {
    double percent = json["AC"]["percent"];
    double seconds = json["AC"]["seconds"];
    print_resume(percent, x);

    if (json["AC"].count("sum") > 0)
    {
      sum = (T) to_maxint(json["AC"]["sum"]);
      time = get_time() - seconds;
      return true;
    }
  }

  return false;
}

/// Compute the A formula.
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename Primes>
T A(T x,
    T xlow,
    T xhigh,
    uint64_t y,
    uint64_t b,
    const Primes& primes,
    const PiTable& pi,
    const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t prime = primes[b];
  T xp = x / prime;
  uint64_t sqrt_xp = (uint64_t) isqrt(xp);
  uint64_t min_2nd_prime = min(xhigh / prime, sqrt_xp);
  uint64_t max_2nd_prime = min(xlow / prime, sqrt_xp);
  uint64_t i = pi[max(prime, min_2nd_prime)] + 1;
  uint64_t max_i1 = pi[min(xp / y, max_2nd_prime)];
  uint64_t max_i2 = pi[max_2nd_prime];

  // x / (p * q) >= y
  for (; i <= max_i1; i++)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq];
  }

  // x / (p * q) < y
  for (; i <= max_i2; i++)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq] * 2;
  }

  return sum;
}

/// Compute the 1st part of the C formula.
/// pi[(x/z)^(1/3)] < b <= pi[sqrt(z)]
/// x / (primes[b] * m) <= z
///
/// m may be a prime <= y or a square free number <= z which is
/// coprime to the first b primes and whose largest prime factor <= y.
/// This algorithm recursively iterates over the square free numbers
/// coprime to the first b primes. This algorithm is described in
/// section 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March 2015.
///
template <int MU, 
          typename T, 
          typename Primes>
T C1(T xp,
     uint64_t b,
     uint64_t i,
     uint64_t pi_y,
     uint64_t m,
     uint64_t min_m,
     uint64_t max_m,
     const Primes& primes,
     const PiTable& pi)
{
  T sum = 0;

  for (i++; i <= pi_y; i++)
  {
    // Calculate next m
    T m128 = (T) m * primes[i];
    if (m128 > max_m)
      return sum;

    uint64_t m64 = (uint64_t) m128;

    if (m64 > min_m) {
      uint64_t xpm = fast_div64(xp, m64);
      T phi_xpm = pi[xpm] - b + 2;
      sum += phi_xpm * MU;
    }

    sum += C1<-MU>(xp, b, i, pi_y, m64, min_m, max_m, primes, pi);
  }

  return sum;
}

template <typename T,
          typename Primes>
T C1(T x,
     uint64_t z,
     uint64_t b,
     uint64_t pi_y,
     const Primes& primes,
     const PiTable& pi)
{
  uint64_t prime = primes[b];
  T xp = x / prime;
  uint64_t max_m = min(xp / prime, z);
  T min_m128 = max(xp / (prime * prime), z / prime);
  uint64_t min_m = min(min_m128, max_m);

  return C1<-1>(xp, b, b, pi_y, 1, min_m, max_m, primes, pi);
}

/// Compute the 2nd part of the C formula.
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) < x^(1/2)
///
template <typename T,
          typename Primes>
T C2(T x,
     T xlow,
     T xhigh,
     uint64_t y,
     uint64_t b,
     const Primes& primes,
     const PiTable& pi,
     const SegmentedPiTable& segmentedPi)
{
  T sum = 0;

  uint64_t prime = primes[b];
  T xp = x / prime;
  uint64_t max_m = min3(xlow / prime, xp / prime, y);
  T min_m128 = max3(xhigh / prime, xp / (prime * prime), prime);
  uint64_t min_m = min(min_m128, max_m);
  uint64_t i = pi[max_m];
  uint64_t pi_min_m = pi[min_m];
  uint64_t min_clustered = (uint64_t) isqrt(xp);
  min_clustered = in_between(min_m, min_clustered, max_m);
  uint64_t pi_min_clustered = pi[min_clustered];

  // Find all clustered easy leaves where
  // successive leaves are identical.
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
  while (i > pi_min_clustered)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    uint64_t phi_xpq = segmentedPi[xpq] - b + 2;
    uint64_t xpq2 = fast_div64(xp, primes[b + phi_xpq - 1]);
    uint64_t i2 = pi[max(xpq2, min_clustered)];
    sum += phi_xpq * (i - i2);
    i = i2;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
  for (; i > pi_min_m; i--)
  {
    uint64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq] - b + 2;
  }

  return sum;
}

/// Compute A + C
template <typename T,
          typename Primes>
T AC_OpenMP(T x,
            int64_t y,
            int64_t z,
            int64_t k,
            int64_t x_star,
            int64_t max_a_prime,
            const Primes& primes,
            int threads,
            double& time)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);
  StatusAC status(x);

  // PiTable's size = z because of the C1 formula.
  // PiTable is accessed much less frequently than
  // SegmentedPiTable, hence it is OK that PiTable's size
  // is fairly large and does not fit into the CPU's cache.
  PiTable pi(max(z, max_a_prime), threads);

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_x13 = pi[x13];
  int64_t pi_root3_xy = pi[iroot<3>(x / y)];
  int64_t pi_root3_xz = pi[iroot<3>(x / z)];
  int64_t min_b = max(k, pi_root3_xz) + 1;

  auto backup_time = get_time();
  auto json = load_backup();
  auto copy = json;
  int resume_threads = calculate_resume_threads(json, "AC");
  int64_t next_b = 0;
  int64_t low = 0;

  if (resume_a_c2(json, x, y, z, k, next_b, low, sum, time))
  {
    // Currently our implementation cannot accurately estimate
    // the completion percentage near the end. In order to fix
    // this we would need to implement backwards sieving in
    // SegmentedPi.cpp which would add some more complexity. I
    // don't think this is worth doing.
    if (low > 0)
    {
      if (json["AC"].count("percent") > 0)
      {
        double percent = json["AC"]["percent"];
        status.setPercent(percent);
      }
    }
  }
  else
  {
    if (!resume_c1(json, x, y, z, k, next_b, sum, time))
      if (json.find("AC") != json.end())
        json.erase("AC");

    next_b = max(next_b, min_b);

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++)
    {
      int64_t b = 0;

      // 1st resume computations from backup file
      for (int j = i; j < resume_threads; j += threads)
      {
        if (resume(copy, x, y, z, k, b, j))
        {
          T sum_thread = C1(x, z, b, pi_y, primes, pi);
          string thread_id = "thread" + to_str(j);

          #pragma omp critical (ac)
          {
            sum -= sum_thread;
            update(json, thread_id, "sum_c1", sum);
          }
        }
      }

      T sum_thread = 0;
      string thread_id = "thread" + to_str(i);

      // 2nd, run new computations
      while (true)
      {
        if (is_print())
          status.print(b, pi_x13);

        #pragma omp critical (ac)
        {
          sum -= sum_thread;
          b = next_b++;

          update(json, thread_id, "sum_c1", b, next_b, pi_sqrtz, sum);

          if (is_backup(backup_time))
          {
            double percent = status.getPercent(next_b, pi_x13);
            backup(json, x, y, z, k, x_star, percent, time);
            backup_time = get_time();
          }
        }

        if (b > pi_sqrtz)
          break;

        sum_thread = C1(x, z, b, pi_y, primes, pi);
      }
    }

    // Resume finished, reset vars so there is
    // no reset attempt the next iteration.
    resume_threads = 0;
    next_b = 0;
    copy.clear();
  }

  // SegmentedPiTable is accessed very frequently.
  // In order to get good performance it is important that
  // SegmentedPiTable fits into the CPU's cache.
  // Hence we use a small size of x^(1/3).
  SegmentedPiTable segmentedPi(low, sqrtx, x13, threads);

  // This computes A and the 2nd part of the C formula.
  // Find all special leaves of type:
  // x / (primes[b] * primes[i]) <= x^(1/2)
  // with z^(1/2) < primes[b] <= x^(1/3).
  // Since we need to lookup PrimePi[n] values for n <= x^(1/2)
  // we use a segmented PrimePi[n] table of size z (~O(x^1/3))
  // in order to reduce the memory usage.
  for (; !segmentedPi.finished(); segmentedPi.next())
  {
    // Current segment [low, high[
    low = segmentedPi.low();
    int64_t high = segmentedPi.high();
    json["AC"]["low"] = segmentedPi.low();

    low = max(low, 1);
    T xlow = x / low;
    T xhigh = x / high;

    // Lower bounds of C2 formula
    min_b = max(k, pi_root3_xy);
    min_b = max(min_b, pi_sqrtz);
    min_b = max(min_b, pi[isqrt(low)]);
    min_b = max(min_b, pi[min(xhigh / y, x_star)]);
    min_b += 1;
    next_b = max(next_b, min_b);

    int64_t min_a = min(xhigh / high, x13);
    min_a = pi[max(x_star, min_a)] + 1;
    if (next_b == pi_x_star + 1) next_b = min_a;

    // Upper bound of A & C2 formulas:
    // x / (p * q) >= low
    // p * next_prime(p) <= x / low
    // p <= sqrt(x / low)
    T sqrt_xlow = isqrt(xlow);
    int64_t max_b = pi[min(sqrt_xlow, x13)];

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < threads; i++)
    {
      int64_t b = 0;

      // 1st resume computations from backup file
      for (int j = i; j < resume_threads; j += threads)
      {
        if (resume(copy, x, y, z, k, b, j))
        {
          T sum_thread = 0;
          string thread_id = "thread" + to_str(j);

          if (b <= pi_x_star)
            sum_thread = C2(x, xlow, xhigh, y, b, primes, pi, segmentedPi);
          else
            sum_thread = A(x, xlow, xhigh, y, b, primes, pi, segmentedPi);

          #pragma omp critical (ac)
          {
            sum += sum_thread;
            update(json, thread_id, "sum_ac", sum);
          }
        }
      }

      T sum_thread = 0;
      string thread_id = "thread" + to_str(i);

      // 2nd, run new computations
      while (true)
      {
        if (is_print())
          status.print(b, max_b);

        #pragma omp critical (ac)
        {
          sum += sum_thread;
          if (next_b == pi_x_star + 1) next_b = min_a;
          b = next_b++;

          update(json, thread_id, "sum_ac", b, next_b, max_b, sum);

          if (is_backup(backup_time))
          {
            double percent = status.getPercent(next_b, max_b);
            backup(json, x, y, z, k, x_star, percent, time);
            backup_time = get_time();
          }
        }

        if (b > max_b)
          break;

        if (b <= pi_x_star)
          sum_thread = C2(x, xlow, xhigh, y, b, primes, pi, segmentedPi);
        else
          sum_thread = A(x, xlow, xhigh, y, b, primes, pi, segmentedPi);
      }
    }

    // Resume finished, reset vars so there is
    // no reset attempt the next iteration.
    resume_threads = 0;
    next_b = 0;
    copy.clear();
    status.next();
  }

  backup_result(json, x, y, z, k, x_star, sum, time);

  return sum;
}

} // namespace

namespace primecount {

int64_t AC(int64_t x,
           int64_t y,
           int64_t z,
           int64_t k,
           int threads)
{
#ifdef ENABLE_MPI
  if (mpi_num_procs() > 1)
    return AC_mpi(x, y, z, k, threads);
#endif

  print("");
  print("=== AC(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  int64_t sum = 0;

  auto json = load_backup();

  if (!resume(json, x, y, z, k, sum, time))
  {
    auto primes = generate_primes<uint32_t>(max_prime);
    sum = AC_OpenMP((uint64_t) x, y, z, k, x_star, max_a_prime, primes, threads, time);
  }

  print("A + C", sum, time);
  return sum;
}

#ifdef HAVE_INT128_T

int128_t AC(int128_t x,
            int64_t y,
            int64_t z,
            int64_t k,
            int threads)
{
#ifdef ENABLE_MPI
  if (mpi_num_procs() > 1)
    return AC_mpi(x, y, z, k, threads);
#endif

  print("");
  print("=== AC(x, y) ===");
  print_gourdon_vars(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  int128_t sum = 0;

  auto json = load_backup();

  if (!resume(json, x, y, z, k, sum, time))
  {
    // uses less memory
    if (max_prime <= numeric_limits<uint32_t>::max())
    {
      auto primes = generate_primes<uint32_t>(max_prime);
      sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, max_a_prime, primes, threads, time);
    }
    else
    {
      auto primes = generate_primes<uint64_t>(max_prime);
      sum = AC_OpenMP((uint128_t) x, y, z, k, x_star, max_a_prime, primes, threads, time);
    }
  }

  print("A + C", sum, time);
  return sum;
}

#endif

} // namespace

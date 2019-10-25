///
/// @file  AC_libdivide.cpp
/// @brief Implementation of the A + C formulas in Xavier Gourdon's
///        prime counting algorithm. In this version the memory usage
///        has been reduced from O(x^(1/2)) to O(z) by segmenting
///        the pi[x] lookup table. In each segment we process the
///        leaves that satisfy: low <= x / (prime * m) < high.
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
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <SegmentedPiTable.hpp>
#include <primecount-internal.hpp>
#include <fast_div.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <min.hpp>
#include <imath.hpp>
#include <print.hpp>
#include <S2Status.hpp>
#include <calculator.hpp>
#include <json.hpp>

#include <stdint.h>
#include <vector>
#include <string>
#include <iostream>

using namespace std;
using namespace primecount;

namespace {

/// backup to file every 60 seconds
bool is_backup(double time)
{
  double seconds = get_time() - time;
  return seconds > 1;
}

/// backup intermediate result
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t k,
            double percent,
            double time)
{
  json["AC"]["x"] = to_string(x);
  json["AC"]["y"] = y;
  json["AC"]["z"] = z;
  json["AC"]["k"] = k;
  json["AC"]["percent"] = percent;
  json["AC"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// backup result
template <typename T, typename J>
void backup(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t k,
            T ac,
            double time)
{
  if (json.find("AC") != json.end())
    json.erase("AC");

  json["AC"]["x"] = to_string(x);
  json["AC"]["y"] = y;
  json["AC"]["z"] = z;
  json["AC"]["k"] = k;
  json["AC"]["ac"] = to_string(ac);
  json["AC"]["percent"] = 100.0;
  json["AC"]["seconds"] = get_time() - time;

  store_backup(json);
}

/// update backup (without storing to disk)
template <typename T, typename J>
void update(J& json,
            int64_t b,
            int64_t min_b,
            int64_t max_b,
            int64_t thread_id,
            T ac)
{
  string tid = "thread" + to_string(thread_id);

  json["AC"]["min_b"] = min_b;
  json["AC"]["ac"] = to_string(ac);

  if (b <= max_b)
    json["AC"][tid]["b"] = b;
  else
  {
    // finished
    if (json["AC"].find(tid) != json["AC"].end())
      json["AC"].erase(tid);
  }
}

/// resume thread
template <typename T, typename J>
bool resume(J& json,
            T x,
            int64_t y,
            int64_t z,
            int64_t k,
            int64_t& b,
            int thread_id)
{
  if (is_resume(json, "AC", thread_id, x, y, z, k))
  {
    string tid = "thread" + to_string(thread_id);
    b = json["AC"][tid]["b"];
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
            int64_t k,
            int64_t& min_b,
            int64_t& low,
            int64_t& high,
            T& ac,
            double& time)
{
  if (is_resume(json, "AC", x, y, z, k))
  {
    double seconds = json["AC"]["seconds"];
    ac = calculator::eval<T>(json["AC"]["ac"]);
    min_b = json["AC"]["min_b"];
    low = json["AC"]["low"];
    high = json["AC"]["high"];
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
            int64_t k,
            T& ac,
            double& time)
{
  if (is_resume(json, "AC", x, y, z, k))
  {
    double percent = json["AC"]["percent"];
    double seconds = json["AC"]["seconds"];
    print_resume(percent, x);

    if (!json["AC"].count("min_b"))
    {
      ac = calculator::eval<T>(json["AC"]["ac"]);
      time = get_time() - seconds;
      return true;
    }
  }

  return false;
}

/// Compute the A formula.
/// pi[x_star] < b <= pi[x^(1/3)]
/// x / (primes[b] * primes[i]) <= x^(1/2)
///
template <typename T, typename Primes>
T A(T x,
    int64_t y,
    int64_t b,
    T x_div_low,
    T x_div_high,
    Primes& primes,
    PiTable& pi,
    SegmentedPiTable& segmentedPi)
{
  int64_t prime = primes[b];
  T xp = x / prime;
  T sum = 0;

  int64_t sqrt_xp = isqrt(xp);
  int64_t min_2nd_prime = min(x_div_high / prime, sqrt_xp);
  int64_t i = pi[min_2nd_prime];
  i = max(i, b) + 1;
  int64_t max_2nd_prime = min(x_div_low / prime, sqrt_xp);
  int64_t max_i = pi[max_2nd_prime];

  // x / (p * q) >= y
  for (; i <= max_i; i++)
  {
    int64_t xpq = fast_div64(xp, primes[i]);
    if (xpq < y)
      break;
    sum += segmentedPi[xpq];
  }

  // x / (p * q) < y
  for (; i <= max_i; i++)
  {
    int64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq] * 2;
  }

  return sum;
}

/// Compute the 1st part of the C formula.
/// k < b <= pi[sqrt(z)]
/// x / (primes[b] * m) <= z
/// 
/// Recursively iterate over the square free numbers coprime
/// to the first b primes. This algorithm is described in
/// section 2.2 of the paper: Douglas Staple, "The Combinatorial
/// Algorithm For Computing pi(x)", arXiv:1503.01839, 6 March
/// 2015.
///
template <int MU, typename T, typename Primes>
T C1(T xp,
     int64_t b,
     int64_t i,
     int64_t pi_y,
     int64_t m,
     int64_t min_m,
     int64_t max_m,
     Primes& primes,
     PiTable& pi)
{
  T sum = 0;

  for (i++; i <= pi_y; i++)
  {
    // Calculate next m
    T m128 = (T) m * primes[i];
    if (m128 > (T) max_m)
      return sum;

    int64_t m64 = (int64_t) m128;
    if (m64 > min_m) {
      int64_t xpm = fast_div64(xp, m64);
      sum += MU * (pi[xpm] - b + 2);
    }

    sum += C1<-MU>(xp, b, i, pi_y, m64, min_m, max_m, primes, pi);
  }

  return sum;
}

/// Compute the 2nd part of the C formula.
/// pi[sqrt(z)] < b <= pi[x_star]
/// x / (primes[b] * primes[i]) <= x^(1/2)
///
template <typename T, typename Primes>
T C2(T x,
     int64_t y,
     int64_t b,
     T x_div_low,
     T x_div_high,
     Primes& primes,
     PiTable& pi,
     SegmentedPiTable& segmentedPi)
{
  int64_t prime = primes[b];
  T xp = x / prime;
  T sum = 0;

  int64_t max_m = min3(x_div_low / prime, xp / prime, y);
  T min_m128 = max3(x_div_high / prime, x / ipow<T>(prime, 3), prime);
  int64_t min_m = min(min_m128, max_m);

  int64_t i = pi[max_m];
  int64_t pi_min_m = pi[min_m];
  int64_t min_clustered = (int64_t) isqrt(xp);
  min_clustered = in_between(min_m, min_clustered, max_m);
  int64_t pi_min_clustered = pi[min_clustered];

  // Find all clustered easy leaves where
  // successive leaves are identical.
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
  while (i > pi_min_clustered)
  {
    int64_t xpq = fast_div64(xp, primes[i]);
    int64_t phi_xpq = segmentedPi[xpq] - b + 2;
    int64_t xpq2 = fast_div64(xp, primes[b + phi_xpq - 1]);
    int64_t i2 = segmentedPi[xpq2];
    sum += phi_xpq * (i - i2);
    i = i2;
  }

  // Find all sparse easy leaves where
  // successive leaves are different.
  // n = primes[b] * primes[i]
  // Which satisfy: n > z && primes[i] <= y
  for (; i > pi_min_m; i--)
  {
    int64_t xpq = fast_div64(xp, primes[i]);
    sum += segmentedPi[xpq] - b + 2;
  }

  return sum;
}

/// Compute A + C
template <typename T, typename Primes>
T AC_OpenMP(T x,
            int64_t y,
            int64_t z,
            int64_t k,
            int64_t x_star,
            int64_t max_a_prime,
            Primes& primes,
            int threads,
            double& time)
{
  T sum = 0;
  int64_t x13 = iroot<3>(x);
  int64_t thread_threshold = 1000;
  threads = ideal_num_threads(threads, x13, thread_threshold);

  S2Status status(x);
  PiTable pi(max(z, max_a_prime));
  SegmentedPiTable segmentedPi(isqrt(x), z, threads);

  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[isqrt(z)];
  int64_t pi_x_star = pi[x_star];
  int64_t pi_x13 = pi[x13];
  int64_t pi_root3_xy = pi[iroot<3>(x / y)];
  int64_t pi_root3_xz = pi[iroot<3>(x / z)];
  int64_t min_b = max(k, pi_root3_xz) + 1;

  // This computes the 1st part of the C formula.
  // Find all special leaves of type:
  // x / (primes[b] * m) <= z.
  // m may be a prime <= y or a square free number <= z
  // who is coprime to the first b primes and whose
  // largest prime factor <= y.
  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction(-: sum)
  for (int64_t b = min_b; b <= pi_sqrtz; b++)
  {
    int64_t prime = primes[b];
    T xp = x / prime;
    int64_t max_m = min(xp / prime, z);
    T min_m128 = max(x / ipow<T>(prime, 3), z / prime);
    int64_t min_m = min(min_m128, max_m);

    sum -= C1<-1>(xp, b, b, pi_y, 1, min_m, max_m, primes, pi);

    if (is_print())
      status.print(b, pi_x13);
  }

  auto json = load_backup();
  auto copy = json;
  double backup_time = get_time();
  int64_t resume_min_b = 0;
  int64_t low = 0;
  int64_t high = 0;

  if (!resume(json, x, y, z, k, resume_min_b, low, high, sum, time))
    if (json.find("AC") != json.end())
      json.erase("AC");

  int64_t resume_threads = calculate_resume_threads(json, "AC");
  for (; !segmentedPi.finished() && segmentedPi.low() < low; segmentedPi.next());

  std::cout << "low: " << segmentedPi.low() << std::endl;
  std::cout << "high: " << segmentedPi.high() << std::endl;

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
    high = segmentedPi.high();
    json["AC"]["low"] = segmentedPi.low();
    json["AC"]["high"] = segmentedPi.high();

    low = max(low, 1);
    T x_div_low = x / low;
    T x_div_high = x / high;

    if (resume_min_b != 0)
    {
      min_b = resume_min_b;
      resume_min_b = 0; 
    }
    else
    {
      min_b = max3(k, pi_sqrtz, pi_root3_xy);
      min_b = max(min_b, pi[isqrt(low)]);
      min_b = max(min_b, pi[min(x_div_high / y, x_star)]);
      min_b = min(min_b, pi_x_star) + 1;
    }

    // x / (primes[i] * primes[i+1]) >= low
    // primes[i] * primes[i+1] <= x / low
    // primes[i] <= floor(sqrt(x / low))
    int64_t sqrt_low = min(isqrt(x_div_low), x13);
    int64_t max_b = pi[sqrt_low];
    max_b = max(max_b, pi_x_star);

    #pragma omp parallel for num_threads(threads)
    for (int64_t i = 0; i < threads; i++)
    {
      int64_t b = 0;

      // 1st resume computations from backup file
      for (int64_t j = i; j < resume_threads; j += threads)
      {
        if (resume(copy, x, y, z, k, b, j))
        {
          T sum_thread = 0;

          if (b <= pi_x_star)
            sum_thread = C2(x, y, b, x_div_low, x_div_high, primes, pi, segmentedPi);
          else
            sum_thread = A(x, y, b, x_div_low, x_div_high, primes, pi, segmentedPi);

          #pragma omp critical (ac)
          {
            sum += sum_thread;
            string tid = "thread" + to_string(j);
            json["AC"]["ac"] = to_string(sum);
            json["AC"].erase(tid);
          }
        }
      }

      T sum_thread = 0;

      // 2nd, run new computations
      while (true)
      {
        if (is_print())
          status.print(b, pi_x13);

        #pragma omp critical (ac)
        {
          sum += sum_thread;
          b = min_b++;

          update(json, b, min_b, max_b, i, sum);

          if (is_backup(backup_time))
          {
            double percent = status.getPercent(min_b, pi_x13, min_b, pi_x13);
            backup(json, x, y, z, k, percent, time);
            backup_time = get_time();
          }
        }

        if (b > max_b)
          break;

        if (b <= pi_x_star)
          sum_thread = C2(x, y, b, x_div_low, x_div_high, primes, pi, segmentedPi);
        else
          sum_thread = A(x, y, b, x_div_low, x_div_high, primes, pi, segmentedPi);
      }
    }
  }

  json["AC"]["low"] = segmentedPi.low();
  json["AC"]["high"] = segmentedPi.high();
  backup(json, x, y, z, k, 100.0, time);

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
  print("");
  print("=== AC(x, y) ===");
  print_gourdon(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  auto primes = generate_primes<int32_t>(max_prime);

  int64_t ac = AC_OpenMP((intfast64_t) x, y, z, k, x_star, max_a_prime, primes, threads, time);

  print("A + C", ac, time);
  return ac;
}

#ifdef HAVE_INT128_T

int128_t AC(int128_t x,
            int64_t y,
            int64_t z,
            int64_t k,
            int threads)
{
  print("");
  print("=== AC(x, y) ===");
  print_gourdon(x, y, z, k, threads);

  double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t) isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  int128_t ac;

  // uses less memory
  if (max_prime <= numeric_limits<uint32_t>::max())
  {
    auto primes = generate_primes<uint32_t>(max_prime);
    ac = AC_OpenMP((intfast128_t) x, y, z, k, x_star, max_a_prime, primes, threads, time);
  }
  else
  {
    auto primes = generate_primes<int64_t>(max_prime);
    ac = AC_OpenMP((intfast128_t) x, y, z, k, x_star, max_a_prime, primes, threads, time);
  }

  print("A + C", ac, time);
  return ac;
}

#endif

} // namespace

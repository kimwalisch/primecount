///
/// @file  S2_trivial.cpp
/// @brief Calculate the contribution of the trivial special leaves.
///        Since this can be calculated very quickly using only
///        about O(alpha * n^(1/3)) time, there is no need to use
///        multi-threading.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <calculator.hpp>
#include <generate.hpp>
#include <int128_t.hpp>
#include <json.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;
using namespace primecount;

namespace {

template <typename T>
void backup(T x,
            int64_t y,
            int64_t z,
            T s2_trivial,
            double time)
{
  auto json = load_backup();

  json["S2_trivial"]["x"] = to_string(x);
  json["S2_trivial"]["y"] = y;
  json["S2_trivial"]["z"] = z;
  json["S2_trivial"]["s2_trivial"] = to_string(s2_trivial);
  json["S2_trivial"]["percent"] = 100.0;
  json["S2_trivial"]["seconds"] = get_time() - time;

  store_backup(json);
}

template <typename T>
bool resume(T x,
            int64_t y,
            int64_t z,
            T& s2_trivial,
            double& time)
{
  auto json = load_backup();

  if (is_resume(json, "S2_trivial", x, y, z))
  {
    double percent = json["S2_trivial"]["percent"];
    double seconds = json["S2_trivial"]["seconds"];

    s2_trivial = calculator::eval<T>(json["S2_trivial"]["s2_trivial"]);
    time = get_time() - seconds;
    print_resume(percent, x);
    return true;
  }

  return false;
}

/// Find all trivial leaves: n = primes[b] * primes[l]
/// which satisfy phi(x / n), b - 1) = 1.
/// Hence we only need to calculate their number!
///
template <typename T>
T get_S2_trivial(T x,
                 int64_t y,
                 int64_t z,
                 int64_t c)
{
  PiTable pi(y);
  int64_t pi_y = pi[y];
  int64_t sqrtz = isqrt(z);
  int64_t prime_c = nth_prime(c);
  int64_t start = max(prime_c, sqrtz) + 1;
  primesieve::iterator it(start - 1, y);

  T s2_trivial = 0;
  int64_t prime;

  // Goes up to ~ x^(1/3)
  while ((prime = it.next_prime()) < y)
  {
    T n = (T) prime * prime;
    int64_t xn = (int64_t)(x / n);
    if (xn <= prime) break;
    s2_trivial += pi_y - pi[xn];
  }

  if (prime < y)
  {
    // \sum_{i = pi[prime]}^{pi[y-1]} pi[y] - i
    // Formula above can be calculated using:
    // https://en.wikipedia.org/wiki/Arithmetic_progression
    // sum = n * (a1 + a2) / 2
    T n = (pi[y-1] - pi[prime]) + 1;
    T a1 = pi[y] - pi[y-1];
    T a2 = pi[y] - pi[prime];
    s2_trivial += n * (a1 + a2) / 2;
  }

  return s2_trivial;
}

} // namespace

namespace primecount {

int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c)
{
  print_log("");
  print_log("=== S2_trivial(x, y) ===");
  print_log("Computation of the trivial special leaves");
  print_log(x, y, c, /* threads */ 1);

  double time = get_time();
  int64_t s2_trivial;

  if (!resume(x, y, z, s2_trivial, time))
  {
    s2_trivial = get_S2_trivial(x, y, z, c);
    backup(x, y, z, s2_trivial, time);
  }

  print_log("S2_trivial", s2_trivial, time);
  return s2_trivial;
}

#ifdef HAVE_INT128_T

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c)
{
  print_log("");
  print_log("=== S2_trivial(x, y) ===");
  print_log("Computation of the trivial special leaves");
  print_log(x, y, c, /* threads */ 1);

  double time = get_time();
  int128_t s2_trivial;

  if (!resume(x, y, z, s2_trivial, time))
  {
    s2_trivial = get_S2_trivial(x, y, z, c);
    backup(x, y, z, s2_trivial, time);
  }

  print_log("S2_trivial", s2_trivial, time);
  return s2_trivial;
}

#endif

} // namespace

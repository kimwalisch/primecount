///
/// @file  S2_trivial.cpp
/// @brief Calculate the contribution of the trivial special leaves.
///        Since this can be calculated very quickly using only
///        about O(alpha * n^(1/3)) time, there is no need to use
///        multi-threading.
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

#include <PiTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <int128_t.hpp>
#include <imath.hpp>
#include <print.hpp>

#include <stdint.h>
#include <algorithm>

using namespace std;
using namespace primecount;

namespace {

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

  // For all primes[b] > z^(1/2) && < x^(1/3):
  // (primes[b] < x / primes[b]^2 < y)
  // s2_trivial += pi[y] - pi[x / primes[b]^2]
  while ((prime = it.next_prime()) < y)
  {
    T pp = (T) prime * prime;
    int64_t xpp = (int64_t)(x / pp);
    if (xpp <= prime) break;
    s2_trivial += pi_y - pi[xpp];
  }

  // For all primes[b] >= x^(1/3) && < y:
  // (x / primes[b]^2 <= primes[b])
  // s2_trivial += pi[y] - b
  //
  // \sum_{b = pi[prime]}^{pi[y-1]} (pi[y] - b)
  // Formula above can be calculated using:
  // https://en.wikipedia.org/wiki/Arithmetic_progression
  // sum = n * (a1 + a2) / 2
  if (prime < y)
  {
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
  print("");
  print("=== S2_trivial(x, y) ===");
  print("Computation of the trivial special leaves");
  print(x, y, c, /* threads */ 1);

  double time = get_time();
  int64_t s2_trivial = get_S2_trivial(x, y, z, c);

  print("S2_trivial", s2_trivial, time);
  return s2_trivial;
}

#ifdef HAVE_INT128_T

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c)
{
  print("");
  print("=== S2_trivial(x, y) ===");
  print("Computation of the trivial special leaves");
  print(x, y, c, /* threads */ 1);

  double time = get_time();
  int128_t s2_trivial = get_S2_trivial(x, y, z, c);

  print("S2_trivial", s2_trivial, time);
  return s2_trivial;
}

#endif

} // namespace

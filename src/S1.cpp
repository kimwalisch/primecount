///
/// @file  S1.cpp
/// @brief Functions to calculate the contribution of the ordinary
///        leaves in the Lagarias-Miller-Odlyzko and Deleglise-Rivat
///        prime counting algorithms.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S1.hpp>
#include <primecount-internal.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace S1 {

/// Recursively iterate over the square free numbers coprime to the
/// first b primes and calculate the sum of the ordinary leaves.
/// This algorithm is based on section 2.2 of the paper:
/// Douglas Staple, "The Combinatorial Algorithm For Computing pi(x)",
/// arXiv:1503.01839, 6 March 2015.
///
template <int MU, typename T, typename P>
T S1(T x,
     int64_t y,
     int64_t b,
     int64_t c,
     T square_free,
     vector<P>& primes)
{
  T s1 = 0;

  for (b += 1; b < (int64_t) primes.size(); b++)
  {
    T next = square_free * primes[b];
    if (next > y) break;
    s1 += MU * phi_tiny(x / next, c);
    s1 += S1<-MU>(x, y, b, c, next, primes);
  }

  return s1;
}

/// Calculate the contribution of the ordinary leaves in parallel.
/// Run time: O(y * log(log(y))) operations.
/// Space complexity: O(y / log(y)).
///
template <typename X, typename Y>
X S1(X x,
     Y y,
     int64_t c,
     int threads)
{
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  vector<Y> primes = generate_primes<Y>(y);
  X s1 = phi_tiny(x, c);

  #pragma omp parallel for schedule(dynamic) num_threads(threads) reduction (+: s1)
  for (int64_t b = c + 1; b < (int64_t) primes.size(); b++)
  {
    s1 += -1 * phi_tiny(x / primes[b], c);
    s1 += S1<1>(x, y, b, c, (X) primes[b], primes);
  }

  return s1;
}

} // namespace S1
} // namespace

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int threads)
{
  print("");
  print("=== S1(x, y) ===");
  print("Computation of the ordinary leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  int64_t s1 = S1::S1(x, y, c, threads);

  print("S1", s1, time);
  return s1;
}

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int threads)
{
  print("");
  print("=== S1(x, y) ===");
  print("Computation of the ordinary leaves");
  print(x, y, c, threads);

  double time = get_wtime();
  int128_t s1;

  // uses less memory
  if (y <= numeric_limits<uint32_t>::max())
    s1 = S1::S1(x, (uint32_t) y, c, threads);
  else
    s1 = S1::S1(x, y, c, threads);

  print("S1", s1, time);
  return s1;
}

#endif

} // namespace primecount

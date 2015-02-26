///
/// @file  S2_trivial.cpp
/// @brief Calculate the contribution of the trivial special leaves
///        in parallel using OpenMP.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <primesieve.hpp>
#include <generate.hpp>
#include <int128.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace S2_trivial {

/// primes[1] = 2, primes[2] = 3, ...
const int64_t primes[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19 };

template <typename T>
T S2_trivial(T x,
             int64_t y,
             int64_t z,
             int64_t c,
             int threads)
{
  int64_t thread_threshold = 1000000;
  threads = validate_threads(threads, y, thread_threshold);

  PiTable pi(y);
  int64_t pi_y = pi[y];
  int64_t sqrtz = isqrt(z);

  T S2_total = 0;
  double time = get_wtime();

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
  for (int64_t i = 0; i < threads; i++)
  {
    int64_t start = max(primes[c], sqrtz + 1);
    int64_t thread_interval = ceil_div(y - start, threads);
    start += thread_interval * i;
    int64_t stop = min(start + thread_interval, y);
    primesieve::iterator iter(start - 1, stop);
    T prime;

    while ((prime = iter.next_prime()) < stop)
    {
      int64_t xn = (int64_t) max(x / (prime * prime), prime);
      S2_total += pi_y - pi[xn];
    }
  }

  if (print_status())
    print_result("S2_trivial", S2_total, time);

  return S2_total;
}

} // namespace S2_trivial
} // namespace

namespace primecount {

int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_trivial(x, y) ===" << endl;
    cout << "Computation of the trivial special leaves" << endl;
  }

  return S2_trivial::S2_trivial(x, y, z, c, threads);
}

#ifdef HAVE_INT128_T

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_trivial(x, y) ===" << endl;
    cout << "Computation of the trivial special leaves" << endl;
  }

  return S2_trivial::S2_trivial(x, y, z, c, threads);
}

#endif

} // namespace primecount

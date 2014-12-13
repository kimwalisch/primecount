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
#include <primecount-internal.hpp>
#include <int128.hpp>

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace S2_trivial {

template <typename T1, typename T2, typename T3>
T1 S2_trivial(T1 x,
              int64_t y,
              int64_t z,
              int64_t c,
              T2& pi,
              vector<T3>& primes,
              int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S2_trivial(x, y) ===" << endl;
    cout << "Computation of the trivial special leaves" << endl;
  }

  T1 S2_total = 0;
  int64_t pi_y = pi[y];
  int64_t pi_sqrtz = pi[min(isqrt(z), y)];
  double time = get_wtime();

  // Find all trivial leaves: n = primes[b] * primes[l]
  // which satisfy phi(x / n), b - 1) = 1
  #pragma omp parallel for num_threads(threads) reduction(+: S2_total)
  for (int64_t b = max(c, pi_sqrtz + 1); b < pi_y; b++)
  {
    T1 prime = primes[b];
    uint64_t xn = (uint64_t) max(x / (prime * prime), prime);
    S2_total += pi_y - pi[xn];
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
                   vector<int32_t>& pi,
                   vector<int32_t>& primes,
                   int threads)
{
  assert(x >= 0);
  return S2_trivial::S2_trivial((uint64_t) x, y, z, c, pi, primes, threads);
}

int64_t S2_trivial(int64_t x,
                   int64_t y,
                   int64_t z,
                   int64_t c,
                   PiTable& pi,
                   vector<int32_t>& primes,
                   int threads)
{
  assert(x >= 0);
  return S2_trivial::S2_trivial((uint64_t) x, y, z, c, pi, primes, threads);
}

#ifdef HAVE_INT128_T

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    PiTable& pi,
                    vector<uint32_t>& primes,
                    int threads)
{
  assert(x >= 0);
  return S2_trivial::S2_trivial((uint128_t) x, y, z, c, pi, primes, threads);
}

int128_t S2_trivial(int128_t x,
                    int64_t y,
                    int64_t z,
                    int64_t c,
                    PiTable& pi,
                    vector<int64_t>& primes,
                    int threads)
{
  assert(x >= 0);
  return S2_trivial::S2_trivial((uint128_t) x, y, z, c, pi, primes, threads);
}

#endif

} // namespace primecount

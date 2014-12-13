///
/// @file  S1.cpp
/// @brief Functions to calculate the contribution of the ordinary
///        leaves in the Lagarias-Miller-Odlyzko and Deleglise-Rivat
///        prime counting algorithms.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <S1.hpp>
#include <primecount-internal.hpp>
#include <FactorTable.hpp>
#include <PhiTiny.hpp>
#include <generate.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <iostream>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace primecount;

namespace {
namespace S1 {

/// Run time: O(y) operations, O(y) space.
/// @pre is_phi_tiny(c) == true
///
template <typename T>
T S1(T x, int64_t y, int64_t c, int64_t prime_c, vector<int32_t>& lpf, vector<int32_t>& mu, int threads)
{
  if (print_status())
  {
    cout << endl;
    cout << "=== S1(x, y) ===" << endl;
    cout << "Computation of the ordinary leaves" << endl;
  }

  T sum = 0;
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  double time = get_wtime();

  #pragma omp parallel for num_threads(threads) reduction (+: sum)
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > prime_c)
      sum += mu[n] * phi_tiny(x / n, c);

  if (print_status())
    print_result("S1", sum, time);

  return sum;
}

/// This version uses less memory (due to FactorTable compression).
/// Run time: O(y) operations, O(y) space.
/// @pre is_phi_tiny(c) == true
///
template <typename T, typename F>
T S1(T x, int64_t y, int64_t c, int64_t prime_c, FactorTable<F>& factors, int threads)
{
  // the factors lookup table contains only numbers
  // which are coprime to 2, 3, 5 and 7
  if (prime_c <= 7)
  {
    vector<int32_t> mu = generate_moebius(y);
    vector<int32_t> lpf = generate_least_prime_factors(y);
    return S1(x, y, c, prime_c, lpf, mu, threads);
  }

  if (print_status())
  {
    cout << endl;
    cout << "=== S1(x, y) ===" << endl;
    cout << "Computation of the ordinary leaves" << endl;
  }

  T sum = 0;
  int64_t limit = factors.get_index(y);
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  double time = get_wtime();

  #pragma omp parallel for num_threads(threads) reduction (+: sum)
  for (int64_t i = factors.get_index(1); i <= limit; i++)
    if (factors.lpf(i) > prime_c)
      sum += factors.mu(i) * phi_tiny(x / factors.get_number(i), c);

  if (print_status())
    print_result("S1", sum, time);

  return sum;
}

} // namespace S1
} // namespace

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int64_t prime_c,
           vector<int32_t>& lpf,
           vector<int32_t>& mu,
           int threads)
{
  return S1::S1((intfast64_t) x, y, c, prime_c, lpf, mu, threads);
}

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int64_t prime_c,
           FactorTable<uint16_t>& factors,
           int threads)
{
  return S1::S1((intfast64_t) x, y, c, prime_c, factors, threads);
}

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int64_t prime_c,
            vector<int32_t>& lpf,
            vector<int32_t>& mu,
            int threads)
{
  return S1::S1((intfast128_t) x, y, c, prime_c, lpf, mu, threads);
}

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            uint32_t prime_c,
            FactorTable<uint16_t>& factors,
            int threads)
{
  return S1::S1((intfast128_t) x, y, c, prime_c, factors, threads);
}

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int64_t prime_c,
            FactorTable<uint32_t>& factors,
            int threads)
{
  return S1::S1((intfast128_t) x, y, c, prime_c, factors, threads);
}

#endif

} // namespace primecount

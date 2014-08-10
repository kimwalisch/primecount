///
/// @file  S1.hpp
/// @brief Functions to calculate the contribution of the ordinary
///        leaves in the Lagarias-Miller-Odlyzko algorithm.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>
#include <FactorTable.hpp>
#include <PhiTiny.hpp>
#include <pmath.hpp>

#include <stdint.h>
#include <iostream>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace primecount {

/// Run time: O(y) operations, O(y) space.
/// @pre is_phi_tiny(c) == true
///
template <typename T, typename V>
T S1(T x, int64_t y, int64_t c, int64_t prime_c, V& lpf, V& mu, int threads)
{
  T sum = 0;
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  double time1 = get_wtime();

  if (print_status())
  {
    std::cout << std::endl;
    std::cout << "=== S1(x, y) ===" << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "threads = " << threads << std::endl;
  }

  #pragma omp parallel for num_threads(threads) reduction (+: sum)
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > prime_c)
      sum += mu[n] * phi_tiny(x / n, c);

  if (print_status())
  {
    std::cout << "S1 = " << sum << std::endl;
    std::cout << "Seconds: " << get_wtime() - time1 << std::endl;
  }

  return sum;
}

/// This version uses less memory (due to factors lookup table).
/// Run time: O(y) operations, O(y) space.
/// @pre is_phi_tiny(c) == true
///
template <typename T, typename F>
T S1(T x, int64_t y, int64_t c, int64_t prime_c, F& factors, int threads)
{
  // the factors lookup table contains only numbers
  // which are coprime to 2, 3, 5 and 7
  if (prime_c <= 7)
  {
    std::vector<int32_t> mu = generate_moebius(y);
    std::vector<int32_t> lpf = generate_least_prime_factors(y);
    return S1(x, y, c, prime_c, lpf, mu, threads);
  }

  T sum = 0;
  int64_t limit = factors.get_index(y);
  int64_t thread_threshold = ipow(10, 6);
  threads = validate_threads(threads, y, thread_threshold);
  double time1 = get_wtime();

  if (print_status())
  {
    std::cout << std::endl;
    std::cout << "=== S1(x, y) ===" << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "threads = " << threads << std::endl;
  }

  #pragma omp parallel for num_threads(threads) reduction (+: sum)
  for (int64_t i = factors.get_index(1); i <= limit; i++)
    if (factors.lpf(i) > prime_c)
      sum += factors.mu(i) * phi_tiny(x / factors.get_number(i), c);

  if (print_status())
  {
    std::cout << "S1 = " << sum << std::endl;
    std::cout << "Seconds: " << get_wtime() - time1 << std::endl;
  }

  return sum;
}

} // namespace primecount

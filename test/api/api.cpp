///
/// @file   api.cpp
/// @brief  Test primecount's C++ API.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  // Test 64-bit pi(-x)
  int64_t n = -1;
  int64_t res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 0);

  n = -9223372036854775807ll;
  res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 0);

  n = (int64_t) 1e10;
  res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 455052511);

  pc_int128_t n128;
  n128.lo = (uint64_t) 1e9;
  n128.hi = 0;
  pc_int128_t res128 = pi(n128);
  std::cout << "pi(" << n128.lo << ") = " << res128.lo;
  check(res128.lo == 50847534 && res128.hi == 0);

  // Check x >= primecount max x of 10^31.
  // primecount must detect issue and throw an exception.
  try {
    std::cout << "pi(2^114) throws primecount_error: ";
    n128.lo = 0;
    n128.hi = 1ull << 50;
    n128 = pi(n128);
    std::cout << "  ERROR" << std::endl;
    if (n128.lo != 0 || n128.hi != 0)
      std::cout << "ERROR: this code path should be unreachable!" << std::endl;
    return 1;
  }
  catch(primecount_error&)
  {
    std::cout << "  OK" << std::endl;
  }

  n = 455052511;
  res = nth_prime(n);
  std::cout << "nth_prime(" << n << ") = " << res;
  check(res == 9999999967);

  n128.lo = (uint64_t) 1e9;
  n128.hi = 0;
  res128 = nth_prime(n128);
  std::cout << "nth_prime(" << n128.lo << ") = " << res128.lo;
  check(res128.lo == 22801763489 && res128.hi == 0);

  // Check n >= primecount max n of ~ 10^29.
  // primecount must detect issue and throw an exception.
  try {
    std::cout << "nth_prime(2^114) throws primecount_error: ";
    n128.lo = 0;
    n128.hi = 1ull << 50;
    n128 = nth_prime(n128);
    std::cout << "  ERROR" << std::endl;
    if (n128.lo != 0 || n128.hi != 0)
      std::cout << "ERROR: this code path should be unreachable!" << std::endl;
    return 1;
  }
  catch(primecount_error&)
  {
    std::cout << "  OK" << std::endl;
  }

  n = (int64_t) 1e12;
  int64_t a = 78498;
  res = phi(n, a);
  std::cout << "phi(" << n << ", " << a << ") = " << res;
  check(res == 37607833521);

  std::cout << "primecount version: " << primecount_version();
  check(primecount_version() == std::string(PRIMECOUNT_VERSION));

  std::cout << "threads: " << get_num_threads() << std::endl;

  // If multi-threading is disabled setting the
  // number of threads must have no effect.
  if (get_num_threads() <= 1)
  {
    set_num_threads(2);
    std::cout << "new threads: " << get_num_threads();
    check(get_num_threads() == 1);
  }
  else
  {
    set_num_threads(2);
    std::cout << "new threads: " << get_num_threads();
    check(get_num_threads() == 2);
  }

  std::pair<double, double> alphas1 = get_alpha_gourdon(1ull << 50);

  set_double_check(true);
  std::pair<double, double> alphas2 = get_alpha_gourdon(1ull << 50);
  std::cout << "set_double_check(true) alpha_y: " << alphas2.first;
  check(alphas2.first != alphas1.first);
  std::cout << "set_double_check(true) alpha_z: " << alphas2.second;
  check(alphas2.second != alphas1.second);

  set_double_check(false);
  std::pair<double, double> alphas3 = get_alpha_gourdon(1ull << 50);
  std::cout << "set_double_check(false) alpha_y: " << alphas3.first;
  check(alphas3.first == alphas1.first);
  std::cout << "set_double_check(false) alpha_z: " << alphas3.second;
  check(alphas3.second == alphas1.second);

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

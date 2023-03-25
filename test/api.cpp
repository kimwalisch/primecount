///
/// @file   api.cpp
/// @brief  Test primecount's C++ API.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>

#include <stdint.h>
#include <iostream>
#include <string>
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
  std::cout << "primecount version: " << primecount_version() << std::endl;
  std::cout << "threads: " << get_num_threads() << std::endl;
  
  set_num_threads(3);
  std::cout << "new threads: " << get_num_threads() << std::endl;

  // Test 64-bit pi(-1)
  int64_t n = -1;
  int64_t res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 0);

  // Test 128-bit pi(-1)
  std::string in = "-1";
  std::string out = pi(in);
  std::cout << "pi(" << in << ") = " << out;
  check(out == "0");

  n = (int64_t) 1e10;
  res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 455052511);

  n = 455052511;
  res = nth_prime(n);
  std::cout << "nth_prime(" << n << ") = " << res;
  check(res == 9999999967);

  n = (int64_t) 1e12;
  int64_t a = 78498;
  res = phi(n, a);
  std::cout << "phi(" << n << ", " << a << ") = " << res;
  check(res == 37607833521);

  in = "1000000000000";
  out = pi(in);
  std::cout << "pi(" << in << ") = " << out;
  check(out == "37607912018");

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

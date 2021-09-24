///
/// @file   api_c.cpp
/// @brief  Test primecount's C API.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>

#include <stdint.h>
#include <iostream>
#include <string>
#include <cstdlib>

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::cout << "primecount version: " << primecount_version() << std::endl;
  std::cout << "threads: " << primecount_get_num_threads() << std::endl;
  
  primecount_set_num_threads(3);
  std::cout << "new threads: " << primecount_get_num_threads() << std::endl;

  int64_t n = (int64_t) 1e10;
  int64_t res = primecount_pi(n);
  std::cout << "primecount_pi(" << n << ") = " << res;
  check(res == 455052511);

  n = 455052511;
  res = primecount_nth_prime(n);
  std::cout << "primecount_nth_prime(" << n << ") = " << res;
  check(res == 9999999967);

  n = (int64_t) 1e12;
  int64_t a = 78498;
  res = primecount_phi(n, a);
  std::cout << "primecount_phi(" << n << ", " << a << ") = " << res;
  check(res == 37607833521);

  const char* in = "1000000000000";
  char out[32];
  primecount_pi_str(in, out, 32);
  std::cout << "primecount_pi_str(" << in << ") = " << out;
  check(std::string(out) == "37607912018");

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

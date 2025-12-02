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
#include <int128_t.hpp>

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

  // Test 64-bit pi(-x)
  int64_t n = -1;
  int64_t res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 0);

  n = -9223372036854775807ll;
  res = pi(n);
  std::cout << "pi(" << n << ") = " << res;
  check(res == 0);

  std::string in = "-1";
  std::string out = pi(in);
  std::cout << "pi(" << in << ") = " << out;
  check(out == "0");

#ifdef HAVE_INT128_T
  // Test 128-bit pi(-x)
  in = "-1208925819614629174696176";
  out = pi(in);
  std::cout << "pi(" << in << ") = " << out;
  check(out == "0");

  // Test using INT128_MIN+1
  in = "-170141183460469231731687303715884105727";
  out = pi(in);
  std::cout << "pi(" << in << ") = " << out;
  check(out == "0");
#endif

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

  n = 455052511;
  res = nth_prime(n);
  std::cout << "nth_prime(" << n << ") = " << res;
  check(res == 9999999967);

  n128.lo = (uint64_t) 1e9;
  n128.hi = 0;
  res128 = nth_prime(n128);
  std::cout << "nth_prime(" << n128.lo << ") = " << res128.lo;
  check(res128.lo == 22801763489 && res128.hi == 0);

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

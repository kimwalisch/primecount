///
/// @file   api.cpp
/// @brief  Test primecount's C++ API.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>

#include <stdint.h>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;
using namespace primecount;

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  cout << "primecount version: " << primecount_version() << endl;
  cout << "threads: " << get_num_threads() << endl;
  
  set_num_threads(3);
  cout << "new threads: " << get_num_threads() << endl;

  int64_t n = (int64_t) 1e10;
  int64_t res = pi(n);
  cout << "pi(" << n << ") = " << res;
  check(res == 455052511);

  n = 455052511;
  res = nth_prime(n);
  cout << "nth_prime(" << n << ") = " << res;
  check(res == 9999999967);

  n = (int64_t) 1e12;
  int64_t a = 78498;
  res = phi(n, a);
  cout << "phi(" << n << ", " << a << ") = " << res;
  check(res == 37607833521);

  string in("1000000000000");
  string out = pi(in);
  cout << "pi(" << in << ") = " << out;
  check(out == "37607912018");

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

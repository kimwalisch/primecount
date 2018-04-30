///
/// @file  fast_div.cpp
/// @brief Test fast_div(x, y) function
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <fast_div.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

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
  random_device rd;
  mt19937 gen(rd());

  uniform_int_distribution<int32_t> dist_i32(1, std::numeric_limits<int32_t>::max());
  uniform_int_distribution<uint64_t> dist_u64(0, std::numeric_limits<uint64_t>::max());

  static_assert(is_same<fastdiv<int32_t>::type, uint32_t>::value,
                "fastdiv<int32_t>::type != uint32_t");

  static_assert(is_same<fastdiv<uint64_t>::type, uint32_t>::value,
                "fastdiv<uint64_t>::type != uint32_t");

  for (int i = 0; i < 10000; i++)
  {
    uint64_t x = dist_i32(gen);
     int32_t y = dist_i32(gen);
    uint64_t res = fast_div(x, y);

    cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_u64(gen);
    y = dist_i32(gen);
    res = fast_div(x, y);

    cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#ifdef HAVE_INT128_T

  uniform_int_distribution<int128_t> dist_i128(0, prt::numeric_limits<int128_t>::max());

  static_assert(is_same<fastdiv<int128_t>::type, uint64_t>::value,
                "fastdiv<int128_t>::type != uint64_t");

  for (int i = 0; i < 10000; i++)
  {
    int128_t x = dist_u64(gen);
     int32_t y = dist_i32(gen);
    int128_t res = fast_div(x, y);

    cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_i128(gen);
    y = dist_i32(gen);
    res = fast_div(x, y);

    cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#endif

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

///
/// @file  fast_div.cpp
/// @brief Test fast_div(x, y) function
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
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

using namespace primecount;

static_assert(std::is_same<make_smaller<int32_t>::type, uint32_t>::value,
              "make_smaller<int32_t>::type != uint32_t");

#if defined(ENABLE_DIV32)

static_assert(std::is_same<make_smaller<uint64_t>::type, uint32_t>::value,
              "make_smaller<uint64_t>::type != uint32_t");

#else

static_assert(std::is_same<make_smaller<uint64_t>::type, uint64_t>::value,
              "make_smaller<uint64_t>::type != uint64_t");

#endif

#ifdef HAVE_INT128_T

static_assert(std::is_same<make_smaller<int128_t>::type, uint64_t>::value,
              "make_smaller<int128_t>::type != uint64_t");

#endif

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<int32_t> dist_i32(1, std::numeric_limits<int32_t>::max());
  std::uniform_int_distribution<uint64_t> dist_u64(0, std::numeric_limits<uint64_t>::max());

  for (int i = 0; i < 10000; i++)
  {
    uint64_t x = dist_i32(gen);
     int32_t y = dist_i32(gen);
    uint64_t res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_u64(gen);
    y = dist_i32(gen);
    res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#ifdef HAVE_INT128_T

  std::uniform_int_distribution<int128_t> dist_i128(0, std::numeric_limits<int128_t>::max());

  for (int i = 0; i < 10000; i++)
  {
    int128_t x = dist_u64(gen);
     int32_t y = dist_i32(gen);
    int128_t res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_i128(gen);
    y = dist_i32(gen);
    res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

///
/// @file  fast_div.cpp
/// @brief Test fast_div(x, y) function
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <fast_div.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

using namespace primecount;

static_assert(std::is_same<make_smaller<int32_t>::type, uint32_t>::value,
              "make_smaller<int32_t>::type != uint32_t");


static_assert(std::is_same<make_smaller<uint64_t>::type, uint64_t>::value,
              "make_smaller<uint64_t>::type != uint64_t");

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

// This test uses double floating-point division if
// ENABLE_DOUBLE_INTEGER_DIVISION is defined.
NOINLINE uint64_t double_integer_divide(uint64_t x, uint64_t y)
{
  return fast_div(x, y);
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

  for (int i = 0; i < 10000; i++)
  {
    uint64_t x = dist_i32(gen);
     int32_t y = dist_i32(gen);
    uint64_t res = fast_div64(x, y);

    std::cout << "fast_div64(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_u64(gen);
    y = dist_i32(gen);
    res = fast_div64(x, y);

    std::cout << "fast_div64(" << x << ", " << y << ") = " << res;
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

  for (int i = 0; i < 10000; i++)
  {
    int128_t x = dist_u64(gen);
     int32_t y = dist_i32(gen);
     int64_t res = fast_div64(x, y);

    std::cout << "fast_div64(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

#endif

  const std::array<std::array<uint64_t, 3>, 9> doubleDivTests =
  {{
    { 1ull << 53, 1 },
    { 1ull << 53, 2 },
    { 1ull << 53, 1ull << 53 },
    { 1ull << 52, 1ull << 53 },
    { (1ull << 53) - 1, 2 },
    { (1ull << 53) - 3, 7 },
    { (1ull << 53), (1ull << 40) - 1 },
    { (1ull << 53) - 1, (1ull << 40) - 1 },
    { (1ull << 53) - 10000, (1ull << 40) - 1 }
  }};

  // This test uses double floating-point division if
  // ENABLE_DOUBLE_INTEGER_DIVISION is defined.
  // It tests double integer division bounds near 2^53.
  for (const auto& v : doubleDivTests)
  {
    uint64_t x = v[0];
    uint64_t y = v[1];
    uint64_t res = double_integer_divide(x, y);
    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

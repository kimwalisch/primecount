///
/// @file   min.cpp
/// @brief  Test min() and max() functions.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>
#include <macros.hpp>
#include <min.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <string>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  static_assert(!is_comparable_int<int64_t, double>::value,
                "is_comparable_int<int64_t, double> is broken");

  static_assert(is_comparable_int<int64_t, int64_t>::value,
                "is_comparable_int<int64_t, int64_t> is broken");

  static_assert(is_comparable_int<uint64_t, int64_t>::value,
                "is_comparable_int<uint64_t, int64_t> is broken");

  static_assert(is_comparable_int<int64_t, int8_t>::value,
                "is_comparable_int<int64_t, int8_t> is broken");

  static_assert(is_comparable_int<uint64_t, uint8_t>::value,
                "is_comparable_int<uint64_t, uint8_t> is broken");

  static_assert(!is_comparable_int<int8_t, int64_t>::value,
                "is_comparable_int<int8_t, int64_t> is broken");

  uint8_t u8_99 = 99;
  uint8_t u8_100 = 100;
  int64_t i64_99 = 99;
  int64_t i64_100 = 100;
  uint64_t u64_99 = 99;
  uint64_t u64_100 = 100;

  std::cout << "min(100, 99) = " << (int) min(u64_100, u8_99);
  check(min(u64_100, u8_99) == std::min(u64_100, u64_99));

  std::cout << "min(99, 100) = " << (int) min(u64_99, u8_100);
  check(min(u64_99, u8_100) == std::min(u64_99, u64_100));

  std::cout << "min(100, 99) = " << min(u64_100, i64_99);
  check(min(u64_100, i64_99) == std::min(i64_100, i64_99));

  std::cout << "min(99, 100) = " << min(u64_99, i64_100);
  check(min(u64_99, i64_100) == std::min(i64_99, i64_100));

  std::cout << "max(100, 99) = " << max(u64_100, u8_99);
  check(max(u64_100, u8_99) == std::max(u64_100, u64_99));

  std::cout << "max(99, 100) = " << max(u64_99, u8_100);
  check(max(u64_99, u8_100) == std::max(u64_99, u64_100));

  std::cout << "max(100, 99) = " << max(u64_100, i64_99);
  check(max(u64_100, i64_99) == std::max(u64_100, u64_99));

  std::cout << "max(99, 100) = " << max(u64_99, i64_100);
  check(max(u64_99, i64_100) == std::max(u64_99, u64_100));

#if defined(HAVE_INT128_T)

  static_assert(!is_comparable_int<uint128_t, double>::value,
                "is_comparable_int<uint128_t, double> is broken");

  static_assert(is_comparable_int<int128_t, int128_t>::value,
                "is_comparable_int<int128_t, int128_t> is broken");

  static_assert(is_comparable_int<uint128_t, int128_t>::value,
                "is_comparable_int<uint128_t, int128_t> is broken");

  static_assert(is_comparable_int<int128_t, int8_t>::value,
                "is_comparable_int<int128_t, int8_t> is broken");

  static_assert(is_comparable_int<uint128_t, uint8_t>::value,
                "is_comparable_int<uint128_t, uint8_t> is broken");

  static_assert(!is_comparable_int<int8_t, int128_t>::value,
                "is_comparable_int<int8_t, int128_t> is broken");

  int128_t i128_99 = 99;
  int128_t i128_100 = 100;
  uint128_t u128_99 = 99;
  uint128_t u128_100 = 100;

  std::cout << "min(100, 99) = " << (int) min(u128_100, u8_99);
  check(min(u128_100, u8_99) == std::min(u128_100, u128_99));

  std::cout << "min(99, 100) = " << (int) min(u128_99, u8_100);
  check(min(u128_99, u8_100) == std::min(u128_99, u128_100));

  std::cout << "min(100, 99) = " << min(u128_100, i128_99);
  check(min(u128_100, i128_99) == std::min(i128_100, i128_99));

  std::cout << "min(99, 100) = " << min(u128_99, i128_100);
  check(min(u128_99, i128_100) == std::min(i128_99, i128_100));

  std::cout << "max(100, 99) = " << max(u128_100, u8_99);
  check(max(u128_100, u8_99) == std::max(u128_100, u128_99));

  std::cout << "max(99, 100) = " << max(u128_99, u8_100);
  check(max(u128_99, u8_100) == std::max(u128_99, u128_100));

  std::cout << "max(100, 99) = " << max(u128_100, i128_99);
  check(max(u128_100, i128_99) == std::max(u128_100, u128_99));

  std::cout << "max(99, 100) = " << max(u128_99, i128_100);
  check(max(u128_99, i128_100) == std::max(u128_99, u128_100));

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

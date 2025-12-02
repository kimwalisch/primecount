///
/// @file   int128.cpp
/// @brief  Test int128_t and uint128_t types.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>
#include <macros.hpp>

#include <stdint.h>
#include <iostream>

using namespace primecount;

#if defined(HAVE_INT128_T)

#include <cstdlib>
#include <sstream>
#include <string>

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

#endif

int main()
{
#if defined(HAVE_INT128_T)

  static_assert(pstd::numeric_limits<uint128_t>::max() == ~((uint128_t) 0),
                "pstd::numeric_limits<uint128_t>::max() is broken");

  static_assert(pstd::is_integral<int128_t>::value,
                "is_integral<int128_t> != true");

  static_assert(pstd::is_integral<uint128_t>::value,
                "is_integral<uint128_t> != true");

  static_assert(pstd::is_signed<int128_t>::value,
                "is_signed<int128_t> != true");

  static_assert(!pstd::is_signed<uint128_t>::value,
                "is_signed<uint128_t> != false");

  static_assert(!pstd::is_unsigned<int128_t>::value,
                "is_unsigned<int128_t> != false");

  static_assert(pstd::is_unsigned<uint128_t>::value,
                "is_unsigned<uint128_t> != true");

  {
    std::ostringstream s;
    s << (((int128_t) 1 ) << 100);
    std::cout << "2^100 = " << s.str();
    check(s.str() == "1267650600228229401496703205376");
  }

  {
    std::ostringstream s;
    s << pstd::numeric_limits<int128_t>::min() + 1;
    std::cout << "-2^127+1 = " << s.str();
    check(s.str() == "-170141183460469231731687303715884105727");
  }

  {
    std::ostringstream s;
    s << pstd::numeric_limits<int128_t>::max();
    std::cout << "2^127-1 = " << s.str();
    check(s.str() == "170141183460469231731687303715884105727");
  }

  {
    std::ostringstream s;
    s << pstd::numeric_limits<uint128_t>::max();
    std::cout << "2^128-1 = " << s.str();
    check(s.str() == "340282366920938463463374607431768211455");
  }

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

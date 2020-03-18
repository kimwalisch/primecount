///
/// @file   int128.cpp
/// @brief  Test int128_t and uint128_t types.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <int128_t.hpp>

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <limits>
#include <type_traits>

using namespace std;
using namespace primecount;

#if defined(HAVE_INT128_T)

#include <cstdlib>
#include <sstream>
#include <string>

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

#endif

int main()
{
#if defined(HAVE_INT128_T)

  static_assert(numeric_limits<uint128_t>::max() == ~((uint128_t) 0), 
                "numeric_limits<uint128_t>::max() is broken");

  static_assert(is_integral<int128_t>::value, 
                "is_integral<int128_t> != true");

  static_assert(is_integral<uint128_t>::value, 
                "is_integral<uint128_t> != true");

  static_assert(is_signed<int128_t>::value, 
                "is_signed<int128_t> != true");

  static_assert(!is_signed<uint128_t>::value, 
                "is_signed<uint128_t> != false");

  static_assert(!is_unsigned<int128_t>::value, 
                "is_unsigned<int128_t> != false");

  static_assert(is_unsigned<uint128_t>::value, 
                "is_unsigned<uint128_t> != true");

  {
    ostringstream s;
    s << (((int128_t) 1 ) << 100);
    cout << "2^100 = " << s.str();
    check(s.str() == "1267650600228229401496703205376");
  }

  {
    ostringstream s;
    s << numeric_limits<int128_t>::min() + 1;
    cout << "-2^127+1 = " << s.str();
    check(s.str() == "-170141183460469231731687303715884105727");
  }

  {
    ostringstream s;
    s << numeric_limits<int128_t>::max();
    cout << "2^127-1 = " << s.str();
    check(s.str() == "170141183460469231731687303715884105727");
  }

  {
    ostringstream s;
    s << numeric_limits<uint128_t>::max();
    cout << "2^128-1 = " << s.str();
    check(s.str() == "340282366920938463463374607431768211455");
  }

#endif

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

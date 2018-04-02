///
/// @file   number_of_bits.cpp
/// @brief  Test number_of_bits(x) function.
///         Note that number_of_bits(x) does not count the sign
///         bit e.g. number_of_bits(int64_t) = 63.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <imath.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>

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
  int8_t i8 = 0;
  cout << "number_of_bits(int8_t) = " << (int) number_of_bits(i8);
  check(number_of_bits(i8) == 7);

  uint8_t ui8 = 0;
  cout << "number_of_bits(uint8_t) = " << (int) number_of_bits(ui8);
  check(number_of_bits(ui8) == 8);

  int16_t i16 = 0;
  cout << "number_of_bits(int16_t) = " << number_of_bits(i16);
  check(number_of_bits(i16) == 15);

  uint16_t ui16 = 0;
  cout << "number_of_bits(uint16_t) = " << number_of_bits(ui16);
  check(number_of_bits(ui16) == 16);

  int32_t i32 = 0;
  cout << "number_of_bits(int32_t) = " << number_of_bits(i32);
  check(number_of_bits(i32) == 31);

  uint32_t ui32 = 0;
  cout << "number_of_bits(uint32_t) = " << number_of_bits(ui32);
  check(number_of_bits(ui32) == 32);

  int64_t i64 = 0;
  cout << "number_of_bits(int64_t) = " << number_of_bits(i64);
  check(number_of_bits(i64) == 63);

  uint64_t ui64 = 0;
  cout << "number_of_bits(uint64_t) = " << number_of_bits(ui64);
  check(number_of_bits(ui64) == 64);

#ifdef HAVE_INT128_T

  int128_t i128 = 0;
  cout << "number_of_bits(int128_t) = " << number_of_bits(i128);
  check(number_of_bits(i128) == 127);

  uint128_t ui128 = 0;
  cout << "number_of_bits(uint128_t) = " << number_of_bits(ui128);
  check(number_of_bits(ui128) == 128);

#endif

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

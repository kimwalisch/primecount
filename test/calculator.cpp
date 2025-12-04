///
/// @file   calculator.cpp
/// @brief  test program for calculator.hpp
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <calculator.hpp>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>

// If stdint.h defines the INT128_MAX macro we assume the int128_t
// and uint128_t are well supported by the C++ library.
#include <stdint.h>

// The compiler supports __int128_t but the C++ standard library
// does not, at least for GCC <= 15 and Clang <= 21.
// Hence we need to write some __int128_t helper functions.
#if defined(__SIZEOF_INT128__) && \
   !defined(INT128_MAX) && \
   !defined(_MSC_VER)

using int128_t = __int128_t;
using uint128_t = __uint128_t;

std::string to_string(uint128_t n)
{
  std::string str;

  while (n > 0)
  {
    str += '0' + n % 10;
    n /= 10;
  }

  if (str.empty())
    str = "0";

  std::reverse(str.begin(), str.end());

  return str;
}

std::string to_string(int128_t n)
{
  if (n >= 0)
    return to_string((uint128_t) n);
  else
  {
    // -n causes undefined behavior for n = INT128_MIN.
    // Hence we use the defined two's complement negation: ~n + 1.
    // Casting ~n to unsigned ensures the result of the addition
    // (2^127 for INT128_MIN) is safely stored in a uint128_t
    // without signed overflow.
    uint128_t abs_n = uint128_t(~n) + 1;
    return "-" + to_string(abs_n);
  }
}

std::ostream& operator<<(std::ostream& stream, int128_t n)
{
  stream << to_string(n);
  return stream;
}

std::ostream& operator<<(std::ostream& stream, uint128_t n)
{
  stream << to_string(n);
  return stream;
}

#endif

template <typename T>
void compare(T result, const std::string& str)
{
  T r = calculator::eval<T>(str);
  std::cout << (r == result ? "Correct: " : "Error: ");
  std::cout << std::setw(50) << str << " = " << std::setw(10) << r;
  if (r != result)
  {
    std::cout << " != " << result << std::endl;
    std::exit(1);
  }
  std::cout << std::endl;
}

void signed_integer_tests()
{
  std::cout << std::endl;
  std::cout << "=== Signed integer tests ===" << std::endl;
  std::cout << std::endl;

  #define STR1(s) #s
  #define TOSTRING(s) STR1(s)

  /// Test expressions
  #define EXPR1 45345 + 0 + 0xdf234 - 1000 % 7
  #define EXPR2 (0 + 0xdf234 - 1000) * 3 / 2 % 999
  #define EXPR3 1 << 16
  #define EXPR4 (0 + ~(0xdf234 & 1000) * 3) / -2
  #define EXPR5 ((1 << 16) + (1 << 16)) >> 0X5
  #define EXPR6 1+(((2+(3+(4+(5+6)* -7)/8))&127)<<1) *-3
  #define EXPR7 100000000 + (1 << 16) + (1 << 16)
  #define EXPR8 1-~1
  #define EXPR9 1- ~1*0xfFa/( ((((8+(6|(4 *(2*(1)*3)*5)|7)+9)))))
  #define EXPRa ((12|13)<<8)>>((1|127) %10&(31+7))
  #define EXPRb ((((((((((5))))))  ))))- ((((((((( 6)))))))))

  compare(EXPR1, TOSTRING(EXPR1));
  compare(EXPR2, TOSTRING(EXPR2));
  compare(EXPR3, TOSTRING(EXPR3));
  compare(EXPR4, TOSTRING(EXPR4));
  compare(EXPR5, TOSTRING(EXPR5));
  compare(EXPR6, TOSTRING(EXPR6));
  compare(EXPR7, TOSTRING(EXPR7));
  compare(EXPR8, TOSTRING(EXPR8));
  compare(EXPR9, TOSTRING(EXPR9));
  compare(EXPRa, TOSTRING(EXPRa));
  compare(EXPRb, TOSTRING(EXPRb));

  std::cout << std::endl;

  compare(calculator::eval<int64_t>("300+(-200)"), "100");
  compare(calculator::eval<int64_t>("300-(-200)"), "500");
  compare(calculator::eval<int64_t>("1e18"), "1000000000000000000");
  compare(calculator::eval<int64_t>("3e18"), "3000000000000000000");
  compare(calculator::eval<int64_t>("10^0"), "1");
  compare(calculator::eval<int64_t>("10^1"), "10");
  compare(calculator::eval<int64_t>("37^2"), "1369");
  compare(calculator::eval<int64_t>("101^3"), "1030301");
  compare(calculator::eval<int64_t>("3^30"), "205891132094649");
  compare(calculator::eval<int64_t>("2^62-1"), "4611686018427387903");
  compare(calculator::eval<int64_t>("2^62-1+2^62"), "9223372036854775807");
  compare(calculator::eval<int64_t>("-(2^62)-(2^62)"), "-9223372036854775807-1");

  std::cout << std::endl;

  try {
    compare(calculator::eval<int64_t>("0xfffffffffffffffffff"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("1000000000000000000000000000"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("10^20"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("123456789012345*1234567890"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("9223372036854775700+200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("-9223372036854775700+(-200)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("-9223372036854775700-200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("9223372036854775700-(-200)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int64_t>("-(-9223372036854775807-1)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

#if defined(__SIZEOF_INT128__) && \
   !defined(_MSC_VER)

  std::cout << std::endl;

  compare(calculator::eval<int128_t>("1e25"), "10000000000000000000000000");
  compare(calculator::eval<int128_t>("3e25"), "30000000000000000000000000");
  compare(calculator::eval<int128_t>("5^50"), "88817841970012523233890533447265625");
  compare(calculator::eval<int128_t>("2^120-1"), "1329227995784915872903807060280344575");
  compare(calculator::eval<int128_t>("2^126-1+2^126"), "170141183460469231731687303715884105727");
  compare(calculator::eval<int128_t>("-(2^126)-(2^126)"), "-170141183460469231731687303715884105727-1");

  std::cout << std::endl;

  try {
    compare(calculator::eval<int128_t>("0xfffffffffffffffffffffffffffffffff"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("10000000000000000000000000000000000000000"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("10^40"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("170141183460469231731687303715884105700*2"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("170141183460469231731687303715884105700+200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("-170141183460469231731687303715884105700+(-200)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("-170141183460469231731687303715884105700-200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("170141183460469231731687303715884105700-(-200)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<int128_t>("-(-170141183460469231731687303715884105727-1)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

#endif
}

void unsigned_integer_tests()
{
  std::cout << std::endl;
  std::cout << "=== Unsigned integer tests ===" << std::endl;
  std::cout << std::endl;

  compare(calculator::eval<uint64_t>("300-200"), "100");
  compare(calculator::eval<uint64_t>("1e19"), "10000000000000000000");
  compare(calculator::eval<uint64_t>("11e18"), "11000000000000000000");
  compare(calculator::eval<uint64_t>("10^0"), "1");
  compare(calculator::eval<uint64_t>("10^1"), "10");
  compare(calculator::eval<uint64_t>("37^2"), "1369");
  compare(calculator::eval<uint64_t>("101^3"), "1030301");
  compare(calculator::eval<uint64_t>("3^30"), "205891132094649");
  compare(calculator::eval<uint64_t>("2^63-1"), "9223372036854775807");
  compare(calculator::eval<uint64_t>("2^63-1+2^63"), "18446744073709551615");
  compare(calculator::eval<uint64_t>("0"), "100-50-50");

  std::cout << std::endl;

  try {
    compare(calculator::eval<uint64_t>("0xfffffffffffffffffff"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint64_t>("1000000000000000000000000000"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint64_t>("10^20"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint64_t>("123456789012345*1234567890"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint64_t>("18446744073709551516+200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint64_t>("2-3"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint64_t>("-100+200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

#if defined(__SIZEOF_INT128__) && \
   !defined(_MSC_VER)

  std::cout << std::endl;

  compare(calculator::eval<uint128_t>("1e25"), "10000000000000000000000000");
  compare(calculator::eval<uint128_t>("3e25"), "30000000000000000000000000");
  compare(calculator::eval<uint128_t>("5^50"), "88817841970012523233890533447265625");
  compare(calculator::eval<uint128_t>("2^120-1"), "1329227995784915872903807060280344575");
  compare(calculator::eval<uint128_t>("2^127-1+2^127"), "340282366920938463463374607431768211455");

  std::cout << std::endl;

  try {
    compare(calculator::eval<uint128_t>("0xfffffffffffffffffffffffffffffffff"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint128_t>("10000000000000000000000000000000000000000"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint128_t>("10^40"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint128_t>("340282366920938463463374607431768211356*2"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint128_t>("340282366920938463463374607431768211356+200"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint128_t>("340282366920938463463374607431768211356-340282366920938463463374607431768211357"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

  try {
    compare(calculator::eval<uint128_t>("100-(-100)"), "0");
  }
  catch (calculator::error& e) {
    std::string errorMsg(e.what());

    if (errorMsg.rfind("Error: ", 0) == 0)
      errorMsg = errorMsg.substr(7);

    std::cout << "Correct: " << errorMsg << std::endl;
  }

#endif
}

int main()
{
  std::cout.setf(std::ios::left);

  signed_integer_tests();
  unsigned_integer_tests();

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

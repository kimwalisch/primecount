///
/// @file   isqrt.cpp
/// @brief  Test integer square root function.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <isqrt.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <calculator.hpp>

#include <stdint.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

using std::numeric_limits;
using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  uint64_t n;
  uint64_t res1;
  double res2;

  for (n = 0; n < 100000; n++)
  {
    res1 = isqrt(n);
    res2 = std::sqrt((double) n);
    std::cout << "isqrt(" << n << ") = " << res1;
    check(res1 == (uint64_t) res2);
  }

  n = (1ull << 32) - 1;
  res1 = isqrt(n);
  res2 = std::sqrt((double) n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == (uint64_t) res2);

  n = 1ull << 32;
  res1 = isqrt(n);
  res2 = std::sqrt((double) n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == (uint64_t) res2);

  n = 1000000000000000000ull - 1;
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 999999999);

  n = 1000000000000000000ull;
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 1000000000);

  n = numeric_limits<int8_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 11);

  n = numeric_limits<uint8_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 15);

  n = numeric_limits<int16_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 181);

  n = numeric_limits<uint16_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 255);

  n = numeric_limits<int32_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 46340);

  n = numeric_limits<uint32_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 65535);

  n = numeric_limits<int64_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 3037000499ll);

  n = numeric_limits<uint64_t>::max();
  res1 = isqrt(n);
  std::cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 4294967295ull);

#ifdef HAVE_INT128_T

  for (n = 0; n < 100000; n++)
  {
    res1 = isqrt((int128_t) n);
    res2 = std::sqrt((double) n);
    std::cout << "isqrt(" << n << ") = " << res1;
    check(res1 == (uint64_t) res2);
  }

  int128_t x = ((int128_t) 1) << 100;
  int128_t res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 1ll << 50);

  x -= 1;
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 1125899906842623ll);

  x = ipow<31>((int128_t) 10);
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 3162277660168379ll);

  x = ipow<30>((int128_t) 10);
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 1000000000000000ll);

  x -= 1;
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 999999999999999ll);

  // In my tests the first occurrences where std::sqrt((double) x)
  // is off by more than 1 happened above 10^32. If std::sqrt(x)
  // is off by more than 1 our isqrt(x) function corrects the
  // result using a while loop. Since primecount can only compute
  // pi(x) for x <= 10^31 our isqrt(x) function is guaranteed to
  // execute in O(1) instructions.

  // here std::sqrt((double) x) is 1 too small
  x = calculator::eval<int128_t>("443075998594972078030832658571409090");
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271553ll);

  // here std::sqrt((double) x) is 1 too large
  x = calculator::eval<int128_t>("443075998594972075382716071791084150");
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271551ll);

  // here std::sqrt((double) x) is 38 too small
  x = calculator::eval<int128_t>("443075998594971958032420320541208365");
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271462ll);

  // here std::sqrt((double) x) is 81 too large
  x = calculator::eval<int128_t>("443075998594971969939937761777907585");
  res3 = isqrt(x);
  std::cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271471ll);

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

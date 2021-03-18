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
  uint64_t n;
  uint64_t res1;
  double res2;

  for (n = 0; n < 100000; n++)
  {
    res1 = isqrt(n);
    res2 = sqrt((double) n);
    cout << "isqrt(" << n << ") = " << res1;
    check(res1 == (uint64_t) res2);
  }

  n = (1ull << 32) - 1;
  res1 = isqrt(n);
  res2 = sqrt((double) n);
  cout << "isqrt(" << n << ") = " << res1;
  check(res1 == (uint64_t) res2);

  n = 1ull << 32;
  res1 = isqrt(n);
  res2 = sqrt((double) n);
  cout << "isqrt(" << n << ") = " << res1;
  check(res1 == (uint64_t) res2);

  n = 1000000000000000000ull - 1;
  res1 = isqrt(n);
  cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 999999999);

  n = 1000000000000000000ull;
  res1 = isqrt(n);
  cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 1000000000);

  n = 18446744073709551615ull;
  res1 = isqrt(n);
  cout << "isqrt(" << n << ") = " << res1;
  check(res1 == 4294967295ull);

#ifdef HAVE_INT128_T

  for (n = 0; n < 100000; n++)
  {
    res1 = isqrt((int128_t) n);
    res2 = sqrt((double) n);
    cout << "isqrt(" << n << ") = " << res1;
    check(res1 == (uint64_t) res2);
  }

  int128_t x = ((int128_t) 1) << 100;
  int128_t res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 1ull << 50);

  x -= 1;
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 1125899906842623ull);

  x = ipow((int128_t) 10, 31);
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 3162277660168379ull);

  x = ipow((int128_t) 10, 30);
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 1000000000000000ull);

  x -= 1;
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 999999999999999ull);

  // In my tests the first occurrences where std::sqrt((double) x)
  // is off by more than 1 happened above 10^32. If std::sqrt(x)
  // is off by more than 1 our isqrt(x) function corrects the
  // result using a while loop. Since primecount can only compute
  // pi(x) for x <= 10^31 our isqrt(x) function is guaranteed to
  // execute in O(1) instructions.

  // here std::sqrt((double) x) is 1 too small
  x = calculator::eval<__int128_t>("443075998594972078030832658571409090");
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271553ull);

  // here std::sqrt((double) x) is 1 too large
  x = calculator::eval<__int128_t>("443075998594972075382716071791084150");
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271551ull);

  // here std::sqrt((double) x) is 38 too small
  x = calculator::eval<__int128_t>("443075998594971958032420320541208365");
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271462ull);

  // here std::sqrt((double) x) is 81 too large
  x = calculator::eval<__int128_t>("443075998594971969939937761777907585");
  res3 = isqrt(x);
  cout << "isqrt(" << x << ") = " << res3;
  check(res3 == 665639541039271471ull);

#endif

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

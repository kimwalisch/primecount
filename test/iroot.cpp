///
/// @file   iroot.cpp
/// @brief  Test integer nth root function.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <imath.hpp>
#include <int128_t.hpp>

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

  if (sizeof(long double) > 8)
  {
    for (n = 0; n < 100000; n++)
    {
      res1 = iroot<2>(n);
      res2 = sqrt((long double) n);
      cout << "iroot<2>(" << n << ") = " << res1;
      check(res1 == (uint64_t) res2);
    }

    for (n = 0; n < 100000; n++)
    {
      res1 = iroot<3>(n);
      res2 = cbrt((long double) n);
      cout << "iroot<3>(" << n << ") = " << res1;
      check(res1 == (uint64_t) res2);
    }

    for (n = 0; n < 100000; n++)
    {
      res1 = iroot<4>(n);
      res2 = pow((long double) n, 1.0L / 4);
      cout << "iroot<4>(" << n << ") = " << res1;
      check(res1 == (uint64_t) res2);
    }

    for (n = 0; n < 100000; n++)
    {
      res1 = iroot<6>(n);
      res2 = pow((long double) n, 1.0L / 6);
      cout << "iroot<6>(" << n << ") = " << res1;
      check(res1 == (uint64_t) res2);
    }
  }

  n = 18446744073709551615ull;
  res1 = iroot<2>(n);
  cout << "iroot<2>(" << n << ") = " << res1;
  check(res1 == 4294967295ull);

  n = 18446744073709551615ull;
  res1 = iroot<3>(n);
  cout << "iroot<3>(" << n << ") = " << res1;
  check(res1 == 2642245);

  for (uint64_t i = 2000000; i <= 2100000; i++)
  {
    n = ipow(i, 3);
    res1 = iroot<3>(n);
    cout << "iroot<3>(" << n << ") = " << res1;
    check(res1 == i);

    res1 = iroot<3>(n - 1);
    cout << "iroot<3>(" << n - 1 << ") = " << res1;
    check(res1 == i - 1);
  }

  n = 18446744073709551615ull;
  res1 = iroot<4>(n);
  cout << "iroot<4>(" << n << ") = " << res1;
  check(res1 == 65535);

  n = 18446744073709551615ull;
  res1 = iroot<6>(n);
  cout << "iroot<6>(" << n << ") = " << res1;
  check(res1 == 1625);

  for (uint64_t i = 1; i <= 1625; i++)
  {
    n = ipow(i, 6);
    res1 = iroot<6>(n);
    cout << "iroot<6>(" << n << ") = " << res1;
    check(res1 == i);

    res1 = iroot<6>(n - 1);
    cout << "iroot<6>(" << n - 1 << ") = " << res1;
    check(res1 == i - 1);
  }

#ifdef HAVE_INT128_T

  int128_t m;
  int128_t res;

  if (sizeof(long double) > 8)
  {
    for (m = 0; m < 100000; m++)
    {
      res = iroot<2>(m);
      res2 = sqrt((long double) m);
      cout << "iroot<2>(" << m << ") = " << res;
      check(res == (int128_t) res2);
    }

    for (m = 0; m < 100000; m++)
    {
      res = iroot<3>(m);
      res2 = cbrt((long double) m);
      cout << "iroot<3>(" << m << ") = " << res;
      check(res == (int128_t) res2);
    }

    for (m = 0; m < 100000; m++)
    {
      res = iroot<4>(m);
      res2 = pow((long double) m, 1.0L / 4);
      cout << "iroot<4>(" << m << ") = " << res;
      check(res == (int128_t) res2);
    }

    for (m = 0; m < 100000; m++)
    {
      res = iroot<6>(m);
      res2 = pow((long double) m, 1.0L / 6);
      cout << "iroot<6>(" << m << ") = " << res;
      check(res == (int128_t) res2);
    }
  }

  m = ipow((int128_t) 10, 30);
  res = iroot<2>(m);
  cout << "iroot<2>(" << m << ") = " << res;
  check(res == ipow(10ll, 15));

  res = iroot<2>(m - 1);
  cout << "iroot<2>(" << m - 1 << ") = " << res;
  check(res == ipow(10ll, 15) - 1);

  res = iroot<3>(m);
  cout << "iroot<3>(" << m << ") = " << res;
  check(res == ipow(10ll, 10));

  res = iroot<3>(m - 1);
  cout << "iroot<3>(" << m - 1 << ") = " << res;
  check(res == ipow(10ll, 10) - 1);

  res = iroot<6>(m);
  cout << "iroot<6>(" << m << ") = " << res;
  check(res == ipow(10, 5));

  res = iroot<6>(m - 1);
  cout << "iroot<6>(" << m - 1 << ") = " << res;
  check(res == ipow(10, 5) - 1);

  m = ipow((int128_t) 10, 28);
  res = iroot<4>(m);
  cout << "iroot<4>(" << m << ") = " << res;
  check(res == ipow(10ll, 7));

  res = iroot<4>(m - 1);
  cout << "iroot<4>(" << m - 1 << ") = " << res;
  check(res == ipow(10ll, 7) - 1);

#endif

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

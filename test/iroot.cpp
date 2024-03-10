///
/// @file   iroot.cpp
/// @brief  Test integer nth root function.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <imath.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>

using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  {
    uint64_t n;
    uint64_t res;
    uint64_t root;

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<2>(n);
      if ((root+1) * (root+1) == n) root++;
      std::cout << "iroot<2>(" << n << ") = " << res;
      check(res == root);
    }

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<3>(n);
      if ((root+1) * (root+1) * (root+1) == n) root++;
      std::cout << "iroot<3>(" << n << ") = " << res;
      check(res == root);
    }

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<4>(n);
      if ((root+1) * (root+1) * (root+1) * (root+1) == n) root++;
      std::cout << "iroot<4>(" << n << ") = " << res;
      check(res == root);
    }

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<6>(n);
      if ((root+1) * (root+1) * (root+1) * (root+1) * (root+1) * (root+1) == n) root++;
      std::cout << "iroot<6>(" << n << ") = " << res;
      check(res == root);
    }

    n = 18446744073709551615ull;
    res = iroot<2>(n);
    std::cout << "iroot<2>(" << n << ") = " << res;
    check(res == 4294967295ull);

    n = 18446744073709551615ull;
    res = iroot<3>(n);
    std::cout << "iroot<3>(" << n << ") = " << res;
    check(res == 2642245);

    for (uint64_t i = 2000000; i <= 2050000; i++)
    {
      n = ipow<3>(i);
      res = iroot<3>(n);
      std::cout << "iroot<3>(" << n << ") = " << res;
      check(res == i);

      res = iroot<3>(n - 1);
      std::cout << "iroot<3>(" << n - 1 << ") = " << res;
      check(res == i - 1);
    }

    n = 18446744073709551615ull;
    res = iroot<4>(n);
    std::cout << "iroot<4>(" << n << ") = " << res;
    check(res == 65535);

    n = 18446744073709551615ull;
    res = iroot<6>(n);
    std::cout << "iroot<6>(" << n << ") = " << res;
    check(res == 1625);

    for (uint64_t i = 1; i <= 1625; i++)
    {
      n = ipow<6>(i);
      res = iroot<6>(n);
      std::cout << "iroot<6>(" << n << ") = " << res;
      check(res == i);

      res = iroot<6>(n - 1);
      std::cout << "iroot<6>(" << n - 1 << ") = " << res;
      check(res == i - 1);
    }
  }

#ifdef HAVE_INT128_T
  {
    int128_t n;
    int128_t res;
    int128_t root;

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<2>(n);
      if ((root+1) * (root+1) == n) root++;
      std::cout << "iroot<2>(" << n << ") = " << res;
      check(res == root);
    }

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<3>(n);
      if ((root+1) * (root+1) * (root+1) == n) root++;
      std::cout << "iroot<3>(" << n << ") = " << res;
      check(res == root);
    }

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<4>(n);
      if ((root+1) * (root+1) * (root+1) * (root+1) == n) root++;
      std::cout << "iroot<4>(" << n << ") = " << res;
      check(res == root);
    }

    for (root = 0, n = 0; n < 50000; n++)
    {
      res = iroot<6>(n);
      if ((root+1) * (root+1) * (root+1) * (root+1) * (root+1) * (root+1) == n) root++;
      std::cout << "iroot<6>(" << n << ") = " << res;
      check(res == root);
    }

    n = ipow<30>((int128_t) 10);
    res = iroot<2>(n);
    std::cout << "iroot<2>(" << n << ") = " << res;
    check(res == ipow<15>(10ll));

    res = iroot<2>(n - 1);
    std::cout << "iroot<2>(" << n - 1 << ") = " << res;
    check(res == ipow<15>(10ll) - 1);

    res = iroot<3>(n);
    std::cout << "iroot<3>(" << n << ") = " << res;
    check(res == ipow<10>(10ll));

    res = iroot<3>(n - 1);
    std::cout << "iroot<3>(" << n - 1 << ") = " << res;
    check(res == ipow<10>(10ll) - 1);

    res = iroot<6>(n);
    std::cout << "iroot<6>(" << n << ") = " << res;
    check(res == ipow<5>(10));

    res = iroot<6>(n - 1);
    std::cout << "iroot<6>(" << n - 1 << ") = " << res;
    check(res == ipow<5>(10) - 1);

    n = ipow<28>((int128_t) 10);
    res = iroot<4>(n);
    std::cout << "iroot<4>(" << n << ") = " << res;
    check(res == ipow<7>(10ll));

    res = iroot<4>(n - 1);
    std::cout << "iroot<4>(" << n - 1 << ") = " << res;
    check(res == ipow<7>(10ll) - 1);
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

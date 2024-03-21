///
/// @file   ipow.cpp
/// @brief  Test integer raise to power function.
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
  for (uint64_t n = 1; n < 10; n++)
  {
    uint64_t res = ipow<0>(n);
    std::cout << "pow<0>(" << n << ") = " << res;
    check(res == 1);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<1>(n);
    std::cout << "pow<1>(" << n << ") = " << res;
    check(res == n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<2>(n);
    std::cout << "pow<2>(" << n << ") = " << res;
    check(res == n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<3>(n);
    std::cout << "pow<3>(" << n << ") = " << res;
    check(res == n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<4>(n);
    std::cout << "pow<4>(" << n << ") = " << res;
    check(res == n * n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<5>(n);
    std::cout << "pow<5>(" << n << ") = " << res;
    check(res == n * n * n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<6>(n);
    std::cout << "pow<6>(" << n << ") = " << res;
    check(res == n * n * n * n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<7>(n);
    std::cout << "pow<7>(" << n << ") = " << res;
    check(res == n * n * n * n * n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<8>(n);
    std::cout << "pow<8>(" << n << ") = " << res;
    check(res == n * n * n * n * n * n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<9>(n);
    std::cout << "pow<9>(" << n << ") = " << res;
    check(res == n * n * n * n * n * n * n * n * n);
  }

  for (uint64_t n = 0; n < 10; n++)
  {
    uint64_t res = ipow<10>(n);
    std::cout << "pow<10>(" << n << ") = " << res;
    check(res == n * n * n * n * n * n * n * n * n * n);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

///
/// @file   next_power_of_2.cpp
/// @brief  Round up to nearest power of 2.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>

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
  uint64_t res2;

  for (uint64_t i = 0; i < 64; i++)
  {
    n = 1ull << i;
    res1 = next_power_of_2(n);
    res2 = n;
    std::cout << "next_power_of_2(" << n << ") = " << res1;
    check(res1 == res2);
  }

  for (uint64_t i = 0; i < 63; i++)
  {
    n = (1ull << i) + 1;
    res1 = next_power_of_2(n);
    res2 = 1ull << (i + 1);
    std::cout << "next_power_of_2(" << n << ") = " << res1;
    check(res1 == res2);
  }

  for (uint64_t i = 2; i < 64; i++)
  {
    n = (1ull << i) - 1;
    res1 = next_power_of_2(n);
    res2 = 1ull << i;
    std::cout << "next_power_of_2(" << n << ") = " << res1;
    check(res1 == res2);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

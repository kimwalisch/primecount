///
/// @file   pi_cache.cpp
/// @brief  Test the pi_cache(x) function.
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount-internal.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace primecount;

std::vector<int64_t> pix =
{
    0, 0, 1, 2, 2, 3, 3, 4, 4, 4,
    4, 5, 5, 6, 6, 6, 6, 7, 7, 8,
    8, 8, 8, 9, 9, 9, 9, 9, 9, 10,
    10, 11, 11, 11, 11, 11, 11, 12, 12, 12,
    12, 13, 13, 14, 14, 14, 14, 15, 15, 15,
    15, 15, 15, 16, 16, 16, 16, 16, 16, 17,
    17, 18, 18, 18, 18, 18, 18, 19, 19, 19,
    19, 20, 20, 21, 21, 21, 21, 21, 21
};

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  {
    int64_t x = -1;
    int64_t res = pi_cache(x);
    std::cout << "pi_cache(" << x << ") = " << res;
    check(res == 0);
  }

  for (int64_t x = 0; x < (int64_t) pix.size(); x++)
  {
    int64_t res = pi_cache(x);
    std::cout << "pi_cache(" << x << ") = " << res;
    check(res == pix[x]);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

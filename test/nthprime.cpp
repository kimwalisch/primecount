///
/// @file   nthprime.cpp
/// @brief  Test the nth_prime(n) function
///
/// Copyright (C) 2023 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.hpp>
#include <primesieve.hpp>
#include <PiTable.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <random>

using namespace primecount;

void check_equal(int64_t x,
                 int64_t res1,
                 int64_t res2)
{
  bool OK = (res1 == res2);
  std::cout << "nth_prime(" << x << ") = " << res1 << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  int64_t n = 1;
  int64_t prime = 0;
  int64_t limit_small = PiTable::max_cached() + 100;
  primesieve::iterator iter;

  for (; n < limit_small; n++)
  {
    prime = iter.next_prime();
    check_equal(n, nth_prime(n), prime);
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(1, 100000);
  n = limit_small - 1;

  // Test random increment
  for (int64_t i = 0; i < 100; i++)
  {
    int64_t next = dist(gen);
    n += next;
    prime = primesieve::nth_prime(next, prime);
    check_equal(n, nth_prime(n), prime);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

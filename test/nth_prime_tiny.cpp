///
/// @file   nth_prime_tiny.cpp
/// @brief  Test the nth_prime(n) function for tiny values of n.
///         If n <= pi(PiTable::max_cached()) then nth_prime(n) uses
///         a lookup table under the hood and does not use any of
///         the advanced prime counting function implementations.
///         Large nth_prime(n) computations are tested in
///         test/api/nth_prime.cpp.
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

int main()
{
  int64_t n;

  // nth_prime(-1) must throw an exception
  try {
    n = -1;
    int64_t res = nth_prime(n);
    std::cout << "nth_prime(" << n << ") = " << res << "   ERROR\n";
    std::exit(1);
  }
  catch (const primecount_error& e) {
    std::cout << "nth_prime(" << n << ") = OK, caught exception: " << e.what() << "\n";
  }

  // nth_prime(0) must throw an exception
  try {
    n = 0;
    int64_t res = nth_prime(n);
    std::cout << "nth_prime(" << n << ") = " << res << "   ERROR\n";
    std::exit(1);
  }
  catch (const primecount_error& e) {
    std::cout << "nth_prime(" << n << ") = OK, caught exception: " << e.what() << "\n";
  }

  primesieve::iterator iter;
  int64_t max_n_tiny = PiTable::pi_cache(PiTable::max_cached());

  for (n = 1; n <= max_n_tiny; n++)
  {
    int64_t res1 = nth_prime(n);
    int64_t res2 = iter.next_prime();
    bool OK = (res1 == res2);
    std::cout << "nth_prime(" << n << ") = " << res1 << "   " << (OK ? "OK" : "ERROR") << "\n";
    if (!OK)
      std::exit(1);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

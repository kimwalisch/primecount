///
/// @file   nth_prime.cpp
/// @brief  Test the nth_prime(n) function for large values of n.
///         For large computations nth_prime(n) uses either
///         pi_legendre(x), pi_meissel(x) or pi_gourdon(x) under the
///         hood. This test has been moved to the test/api directory
///         so that it is executed after the Legendre, Meissel and
///         Gourdon algorithm tests.
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

void check_equal(int64_t n,
                 int64_t res1,
                 int64_t res2)
{
  bool OK = (res1 == res2);
  std::cout << "nth_prime(" << n << ") = " << res1 << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  primesieve::iterator iter(PiTable::max_cached() + 1);
  int64_t n = PiTable::pi_cache(PiTable::max_cached()) + 1;
  int64_t limit_small = n + 100;
  int64_t prime = iter.next_prime();

  // Test first few n > pi(PiTable::max_cached())
  for (; n < limit_small; n++)
  {
    check_equal(n, nth_prime(n), prime);
    prime = iter.next_prime();
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> dist(1, 10000);

  // Test random increment, goes up to ~ 5*10^6
  for (int64_t i = 0; i < 1000; i++)
  {
    int64_t next_n = n + dist(gen);
    for (; n < next_n; n++)
      prime = iter.next_prime();

    check_equal(n, nth_prime(n), prime);
  }

  // nth_prime(1e7)
  n = 10000000ll;
  check_equal(n, nth_prime(n), 179424673ll);

  // nth_prime(1e8)
  n = 100000000ll;
  check_equal(n, nth_prime(n), 2038074743ll);

  // nth_prime(1e9)
  n = 1000000000ll;
  check_equal(n, nth_prime(n), 22801763489ll);

  // nth_prime(1e10)
  n = 10000000000ll;
  check_equal(n, nth_prime(n), 252097800623ll);

  // nth_prime(1e11)
  n = 100000000000ll;
  check_equal(n, nth_prime(n), 2760727302517ll);

  // nth_prime(1e12)
  n = 1000000000000ll;
  check_equal(n, nth_prime(n), 29996224275833ll);

  // nth_prime(1e13)
  n = 10000000000000ll;
  check_equal(n, nth_prime(n), 323780508946331ll);

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

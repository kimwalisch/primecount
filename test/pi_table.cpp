///
/// @file   pi_table.cpp
/// @brief  Test the PiTable class
/// @link   https://en.wikipedia.org/wiki/Prime-counting_function
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <primesieve.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using std::size_t;
using namespace primecount;

std::vector<int> pix =
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
  // Test PiTable::pi_cache(x)
  {
    primesieve::iterator it;
    int64_t prime = it.next_prime();
    int64_t count = 0;

    for (int64_t i = 0; i <= PiTable::max_cached(); i++)
    {
      if (i == prime)
      {
        count++;
        prime = it.next_prime();
      }
      std::cout << "pi_cache(" << i << ") = " << PiTable::pi_cache(i);
      check(PiTable::pi_cache(i) == count);
    }
  }

  // Test PiTable::pi(x)
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(1000000, 2000000);

    int threads = 1;
    PiTable pi(dist(gen), threads);

    for (size_t i = 0; i < pix.size(); i++)
    {
      std::cout << "pi(" << i << ") = " << pi[i];
      check(pi[i] == pix[i]);
    }

    primesieve::iterator it;
    uint64_t prime = it.next_prime();
    int count = 1;

    while (prime < pi.size())
    {
      std::cout << "pi(" << prime << ") = " << pi[prime];
      check(pi[prime] == count);
      prime = it.next_prime();
      count++;
    }

    for (int i = 0; i < 10000; i++)
    {
      int n = dist(gen) % pi.size();
      std::cout << "pi(" << n << ") = " << pi[n];
      check(pi[n] == (int64_t) primesieve::count_primes(0, n));
    }
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

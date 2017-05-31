///
/// @file   pi_table.cpp
/// @brief  Test the PiTable class
/// @link   https://en.wikipedia.org/wiki/Prime-counting_function
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
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

using namespace std;
using namespace primecount;

vector<int> pix =
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
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(1000000, 2000000);
  PiTable pi(dist(gen));

  for (size_t i = 0; i < pix.size(); i++)
  {
    cout << "pi(" << i << ") = " << pi[i];
    check(pi[i] == pix[i]);
  }

  primesieve::iterator it;
  int64_t prime = it.next_prime();
  int64_t count = 1;

  while (prime < pi.size())
  {
    cout << "pi(" << prime << ") = " << pi[prime];
    check(pi[prime] == count);
    prime = it.next_prime();
    count++;
  }

  for (int i = 0; i < 10000; i++)
  {
    int n = dist(gen) % pi.size();
    cout << "pi(" << n << ") = " << pi[n];
    check(pi[n] == (int64_t) primesieve::count_primes(0, n));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

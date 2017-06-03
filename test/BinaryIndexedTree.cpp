///
/// @file   BinaryIndexedTree.cpp
/// @brief  Test the BinaryIndexedTree class which counts
///         the number of unsieved elements in the sieve
///         array using only O(log n) operations.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BinaryIndexedTree.hpp>
#include <generate.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

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
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(1000000, 2000000);

  int pre_sieve = 13;
  int low = 1;
  int size = dist(gen);
  size = next_power_of_2(size);

  auto primes = generate_primes<int>(isqrt(size));
  vector<int> sieve(size, 1);
  BinaryIndexedTree tree;

  for (size_t i = 1; i < primes.size(); i++)
  {
    for (int j = primes[i] - low; j < size; j += primes[i])
    {
      if (sieve[j] && primes[i] > pre_sieve)
        tree.update(j);
      sieve[j] = 0;
    }

    if (primes[i] <= pre_sieve)
      tree.init(sieve);

    int rand = dist(gen) % size;
    int count = 0;

    for (int j = 0; j <= rand; j++)
      count += sieve[j];

    cout << "tree.count(" << rand << ") = " << tree.count(0, rand);
    check(count == tree.count(0, rand));
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

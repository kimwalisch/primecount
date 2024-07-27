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
#include <generate_primes.hpp>
#include <imath.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

using std::size_t;
using namespace primecount;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(1000000, 2000000);

  int pre_sieve = 13;
  int low = 1;
  int size = dist(gen);
  size = next_power_of_2(size);

  auto primes = generate_primes<int32_t>(isqrt(size));
  std::vector<int> sieve(size, 1);
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

    std::cout << "tree.count(" << rand << ") = " << tree.count(0, rand);
    check(count == tree.count(0, rand));
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

///
/// @file   BinaryIndexedTree.cpp
/// @brief  Test the BinaryIndexedTree class which counts
///         the number of unsieved elements in the sieve
///         array using only O(log n) operations.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
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
  std::uniform_int_distribution<int64_t> dist(1000000, 2000000);

  int64_t pre_sieve = 13;
  int64_t low = 1;
  int64_t size = dist(gen);
  size = next_power_of_2(size);

  auto primes = generate_primes<uint32_t>(isqrt(size));
  std::vector<bool> sieve(size, true);
  BinaryIndexedTree tree;

  for (size_t i = 1; i < primes.size(); i++)
  {
    for (int64_t j = primes[i] - low; j < size; j += primes[i])
    {
      if (sieve[j] && primes[i] > pre_sieve)
        tree.update(j);
      sieve[j] = false;
    }

    if (primes[i] <= pre_sieve)
      tree.init(sieve);

    int64_t rand = dist(gen) % size;
    int64_t count = 0;

    for (int64_t j = 0; j <= rand; j++)
      count += sieve[j];

    std::cout << "tree.count(" << rand << ") = " << tree.count(0, rand);
    check(count == tree.count(0, rand));
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

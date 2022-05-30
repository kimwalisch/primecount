///
/// @file   pod_vector.cpp
/// @brief  Plain old data vector, like std::vector but does not 
///         default initialize memory.
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesieve/pod_vector.hpp>

#include <cstdlib>
#include <iostream>
#include <random>
#include <numeric>
#include <utility>

using std::size_t;
using primesieve::pod_vector;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  // The pod_vector class uses std::vector internally. For
  // performance reasons we want vector::resize() not to
  // free memory when resizing to a smaller size. The C++
  // standard seems to indirectly guarantee this behavior,
  // but it is not 100% clear. So this tests verifies this
  // behavior.

  // Allocate from 1 KiB to 128 MiB
  for (size_t i = 10; i <= 27; i++)
  {
    pod_vector<char> vect;
    vect.resize(size_t(1) << i);
    auto capacity1 = vect.capacity();
    vect.resize(100);
    auto capacity2 = vect.capacity();

    std::cout << "vect.resize(100).capacity = " << capacity1;
    check(capacity1 == capacity2);
  }

  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> dist(10000, 20000);

    std::size_t n = dist(gen);
    pod_vector<int> vect;
    vect.reserve(n);

    std::cout << "Vect size after reserve: " << vect.size();
    check(vect.size() == 0);
    std::cout << "Vect empty after reserve: " << vect.empty();
    check(vect.empty() == true);
    std::cout << "Vect capacity after reserve: " << vect.capacity();
    check(vect.capacity() == n);
  }

  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(10000, 20000);

    int size = dist(gen);
    pod_vector<int> vect(size);
    std::fill_n(&vect[0], size, 123);

    // Test if resize does not default initilize
    vect.resize(0);
    vect.resize(size);
    int sum = std::accumulate(&vect[0], &vect[0] + vect.size(), 0);
    std::cout << "Vect sum after resize: " << sum;
    check(sum == 123 * size);
    std::cout << "Vect.end(): " << vect.end();
    check(vect.end() == vect.begin() + vect.size());

    // Test reallocation (old content must be copied into new vector)
    vect.resize(vect.size() * 2);
    sum = std::accumulate(&vect[0], &vect[0] + size, 0);
    std::cout << "Vect sum after reallocation: " << sum;
    check(sum == 123 * size);
  }

  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(10000, 20000);

    int size = dist(gen);
    pod_vector<int> vect(size);
    std::fill_n(&vect[0], size, 123);

    pod_vector<int> vect2 = std::move(vect);
    std::cout << "Vect1 empty after std::move: " << vect.empty();
    check(vect.empty() == true);
    int sum = std::accumulate(vect2.begin(), vect2.end(), 0);
    std::cout << "Vect2 sum after std::move: " << sum;
    check(sum == 123 * size);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

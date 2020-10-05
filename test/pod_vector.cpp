///
/// @file   vector_resize.cpp
/// @brief  The primecount::pod_vector class uses pod_vector
///         internally. For performance reasons we want
///         vector::resize() not to free memory when resizing to
///         a smaller size. The C++ standard seems to indirectly
///         guarantee this behavior, but it is not 100% clear.
///         So this tests verifies this behavior.
///
///         If this test ever fails update pod_vector::resize()
///         in primecount.hpp to (and delete this test):
///
///         void resize(std::size_t size)
///         {
///           if (size > size_)
///             vect_.resize(size);
///           size_ = size;
///         }
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <pod_vector.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
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
  // Allocate from 1 KiB to 128 MiB
  for (size_t i = 10; i <= 27; i++)
  {
    pod_vector<char> vect;
    vect.resize(size_t(1) << i);
    auto capacity1 = vect.capacity();
    vect.resize(100);
    auto capacity2 = vect.capacity();

    cout << "vect.resize(100).capacity = " << capacity1;
    check(capacity1 == capacity2);
  }

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<size_t> dist(10000, 20000);

  size_t size = dist(gen);
  pod_vector<size_t> vect(size);
  fill_n(&vect[0], size, 123);
  
  // Test if resize does not default initilize
  vect.resize(0);
  vect.resize(size);
  size_t sum = accumulate(&vect[0], &vect[0] + size, 0);
  cout << "Vect sum after resize: " << sum;
  check(sum == 123 * size);

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

///
/// @file   vector_resize.cpp
/// @brief  The primecount::pod_vector class uses std::vector
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

#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

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
    std::vector<char> vect;
    vect.resize(1 << i);
    auto capacity1 = vect.capacity();
    vect.resize(100);
    auto capacity2 = vect.capacity();

    cout << "vect.resize(100).capacity = " << capacity1;
    check(capacity1 == capacity2);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}

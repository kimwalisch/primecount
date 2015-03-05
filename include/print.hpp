///
/// @file  print.hpp
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRINT_HPP
#define PRINT_HPP

#include <primecount-internal.hpp>
#include <pmath.hpp>
#include <stdint.h>

#include <iostream>
#include <iomanip>
#include <string>

namespace primecount {

void print(const std::string& str)
{
  if (print_status())
    std::cout << str << std::endl;
}

template <typename T>
void print(T x, int64_t y, int64_t c, int threads)
{
  if (print_variables())
  {
    T z = x / y;
    double alpha = (double) y / (double) iroot<3>(x);

    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl;
    std::cout << "c = " << PhiTiny::max_a() << std::endl;
    std::cout << "alpha = " << std::fixed << std::setprecision(3) << alpha << std::endl;
    std::cout << "threads = " << validate_threads(threads) << std::endl;
    std::cout << std::endl;
  }
}

} // namespace

#endif

///
/// @file  api_c.cpp
///        primecount's C API.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>
#include <primecount.hpp>

#include <exception>
#include <iostream>

int64_t primecount_pi(int64_t x)
{
  try
  {
    return primecount::pi(x);
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount: " << e.what() << std::endl;
    return -1;
  }
}

int64_t primecount_nth_prime(int64_t n)
{
  try
  {
    return primecount::nth_prime(n);
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount: " << e.what() << std::endl;
    return -1;
  }
}

int64_t primecount_phi(int64_t x, int64_t a)
{
  try
  {
    return primecount::phi(x, a);
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount: " << e.what() << std::endl;
    return -1;
  }
}

int primecount_get_num_threads()
{
  try
  {
    return primecount::get_num_threads();
  }
  catch(const std::exception& e)
  {
    return -1;
  }
}

void primecount_set_num_threads(int threads)
{
  try
  {
    primecount::set_num_threads(threads);
  }
  catch(const std::exception&)
  { }
}

const char* primecount_get_max_x()
{
  // 10^31
  return "10000000000000000000000000000000";
}

const char* primecount_version()
{
  return PRIMECOUNT_VERSION;
}

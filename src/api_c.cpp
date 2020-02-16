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
#include <int128_t.hpp>

#include <stdint.h>
#include <stddef.h>
#include <string>
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

int primecount_pi128(char* x, char* res, size_t len)
{
  try
  {
    // NULL pointer, nothing to do
    if (!x)
      return -1;

    std::string str(x);
    std::string pix = primecount::pi(str);

    // +1 required to add null at the end of the string
    if (len < pix.length() + 1)
    {
      std::cerr << "primecount: res buffer too small, res.len = " << len << " < required = " << pix.length() + 1 << std::endl;
      return -1;
    }

    pix.copy(res, pix.length());
    // std::string::copy does not append a null character
    // at the end of the copied content.
    res[pix.length()] = (char) 0;

    return (int) pix.length();
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
#ifdef HAVE_INT128_T
  // 10^31
  return "10000000000000000000000000000000";
#else
  // 2^63-1
  return "9223372036854775807";
#endif
}

const char* primecount_version()
{
  return PRIMECOUNT_VERSION;
}

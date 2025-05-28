///
/// @file  api_c.cpp
///        primecount's C API.
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primecount.h>
#include <primecount.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <stddef.h>
#include <sstream>
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
    std::cerr << "primecount_pi: " << e.what() << std::endl;
    return -1;
  }
}

pc_int128_t primecount_pi128(pc_int128_t x)
{
  try
  {
    primecount::pc_int128_t n;
    n.lo = x.lo;
    n.hi = x.hi;

    auto pix = primecount::pi(n);

    pc_int128_t res;
    res.lo = pix.lo;
    res.hi = pix.hi;
    return res;
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount_pi128: " << e.what() << std::endl;

    pc_int128_t res;
    res.lo = ~0ull;
    res.hi = -1;
    return res;
  }
}

int primecount_pi_str(const char* x, char* res, size_t len)
{
  try
  {
    if (!x)
      throw primecount::primecount_error("x must not be a NULL pointer");

    if (!res)
      throw primecount::primecount_error("res must not be a NULL pointer");

    std::string str(x);
    std::string pix = primecount::pi(str);

    // +1 required to add null at the end of the string
    if (len < pix.length() + 1)
    {
      std::ostringstream oss;
      oss << "res buffer too small, res.len = " << len << " < required = " << pix.length() + 1;
      throw primecount::primecount_error(oss.str());
    }

    pix.copy(res, pix.length());
    // std::string::copy does not append a null character
    // at the end of the copied content.
    res[pix.length()] = '\0';

    return (int) pix.length();
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount_pi_str: " << e.what() << std::endl;

    if (res && len > 0)
      res[0] = '\0';

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
    std::cerr << "primecount_nth_prime: " << e.what() << std::endl;
    return -1;
  }
}

pc_int128_t primecount_nth_prime128(pc_int128_t n)
{
  try
  {
    primecount::pc_int128_t x;
    x.lo = n.lo;
    x.hi = n.hi;

    auto nth_prime = primecount::nth_prime(x);

    pc_int128_t res;
    res.lo = nth_prime.lo;
    res.hi = nth_prime.hi;
    return res;
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount_nth_prime128: " << e.what() << std::endl;

    pc_int128_t res;
    res.lo = ~0ull;
    res.hi = -1;
    return res;
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
    std::cerr << "primecount_phi: " << e.what() << std::endl;
    return -1;
  }
}

int primecount_get_num_threads(void)
{
  try
  {
    return primecount::get_num_threads();
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount_get_num_threads: " << e.what() << std::endl;
    return -1;
  }
}

void primecount_set_num_threads(int threads)
{
  try
  {
    primecount::set_num_threads(threads);
  }
  catch(const std::exception& e)
  {
    std::cerr << "primecount_set_num_threads: " << e.what() << std::endl;
  }
}

const char* primecount_get_max_x(void)
{
#ifdef HAVE_INT128_T
  // 10^31
  return "10000000000000000000000000000000";
#else
  // 2^63-1
  return "9223372036854775807";
#endif
}

const char* primecount_version(void)
{
  return PRIMECOUNT_VERSION;
}

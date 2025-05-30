///
/// @file   nth_prime_sieve.cpp
/// @brief  Test the <nth_prime_sieve.hpp>
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <nth_prime_sieve.hpp>

#include <stdint.h>
#include <iostream>

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

int main()
{
  NthPrimeSieve<int64_t> nthPrimeSieve;
  uint64_t low = 7;
  uint64_t high = 3089;
  nthPrimeSieve.sieve(low, high);
  std::cout << "nthPrimeSieve.sieve(" << low << ", " << high << ")" << std::endl;
  uint64_t prime = nthPrimeSieve.find_nth_prime(1);
  std::cout << "nthPrimeSieve.find_nth_prime(1) = " << prime;
  check(prime == 7);
  prime = nthPrimeSieve.find_nth_prime(100);
  std::cout << "nthPrimeSieve.find_nth_prime(100) = " << prime;
  check(prime == 563);
  prime = nthPrimeSieve.find_nth_prime(439);
  std::cout << "nthPrimeSieve.find_nth_prime(439) = " << prime;
  check(prime == 3089);
  prime = nthPrimeSieve.find_nth_prime(440);
  std::cout << "nthPrimeSieve.find_nth_prime(440) failed=" << (prime == 0 ? "true" : "false");
  check(prime == 0);

  std::cout << std::endl;

  low = 8;
  high = 3088;
  nthPrimeSieve.sieve(low, high);
  std::cout << "nthPrimeSieve.sieve(" << low << ", " << high << ")" << std::endl;
  prime = nthPrimeSieve.find_nth_prime(1);
  std::cout << "nthPrimeSieve.find_nth_prime(1) = " << prime;
  check(prime == 11);
  prime = nthPrimeSieve.find_nth_prime(99);
  std::cout << "nthPrimeSieve.find_nth_prime(99) = " << prime;
  check(prime == 563);
  prime = nthPrimeSieve.find_nth_prime(437);
  std::cout << "nthPrimeSieve.find_nth_prime(437) = " << prime;
  check(prime == 3083);
  prime = nthPrimeSieve.find_nth_prime(3090);
  std::cout << "nthPrimeSieve.find_nth_prime(3090) failed=" << (prime == 0 ? "true" : "false");
  check(prime == 0);

  std::cout << std::endl;

  low = uint64_t(1e12);
  high = uint64_t(1e12) + (10000 * 30);
  nthPrimeSieve.sieve(low, high);
  std::cout << "nthPrimeSieve.sieve(" << low << ", " << high << ")" << std::endl;
  prime = nthPrimeSieve.find_nth_prime(1);
  std::cout << "nthPrimeSieve.find_nth_prime(1) = " << prime;
  check(prime == 1000000000039);
  prime = nthPrimeSieve.find_nth_prime(100);
  std::cout << "nthPrimeSieve.find_nth_prime(100) = " << prime;
  check(prime == 1000000002889);
  prime = nthPrimeSieve.find_nth_prime(10876);
  std::cout << "nthPrimeSieve.find_nth_prime(10876) = " << prime;
  check(prime == 1000000299997);
  prime = nthPrimeSieve.find_nth_prime(10877);
  std::cout << "nthPrimeSieve.find_nth_prime(10877) failed=" << (prime == 0 ? "true" : "false");
  check(prime == 0);

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

///
/// @file   FactorTableD.cpp
/// @brief  FactorTableD is a compressed lookup table of mu (moebius),
///         lpf (least prime factor) and mpf (max prime factor).
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

// factorTable.mu(n) = 0 is disabled by default for performance
// reasons, we only enable it for testing.
#define ENABLE_MU_0_TESTING

#include <FactorTableD.hpp>
#include <generate.hpp>

#include <stdint.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>

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
  std::uniform_int_distribution<int> dist_y(50000, 60000);
  std::uniform_int_distribution<int> dist_z(1200000, 1500000);

  auto y = dist_y(gen);
  auto z = dist_z(gen);
  auto threads = get_num_threads();
  auto lpf = generate_lpf(z);
  auto mpf = generate_mpf(z);
  auto mu = generate_moebius(z);

  FactorTableD<uint16_t> factorTable(y, z, threads);
  int64_t uint16_max = std::numeric_limits<uint16_t>::max();
  int64_t limit = factorTable.first_coprime();
  std::vector<int> small_primes = { 2, 3, 5, 7, 11, 13, 17, 19 };

  for (int n = 1; n <= z; n++)
  {
    int64_t i = factorTable.to_index(n);
    bool is_prime = (lpf[n] == n);

    // Check if n is coprime to the primes < limit
    for (int p : small_primes)
    {
      if (p >= limit)
        break;
      if (n % p == 0)
        goto not_coprime;
    }

    // primes > y and square free numbers with a prime factor > y
    // have been removed from the FactorTableD.
    if (mpf[n] > y)
    {
      std::cout << "prime_factor_larger_y(" << n << ") = " << (factorTable.is_leaf(i) == 0);
      check(factorTable.is_leaf(i) == 0);
      continue;
    }

    std::cout << "mu(" << n << ") = " << factorTable.mu(i);
    check(mu[n] == factorTable.mu(i));

    std::cout << "lpf(" << n << ") = " << lpf[n];

    // is_leaf(n) is a combination of the mu(n) (Möbius function),
    // lpf(n) (least prime factor) and mpf(n) (max prime factor)
    // functions. is_leaf(n) returns (with n = to_number(index)):
    //
    // 1) INT_MAX - 1  if n = 1
    // 2) INT_MAX      if n is a prime
    // 3) 0            if n has a prime factor > y
    // 4) 0            if moebius(n) = 0
    // 5) lpf - 1      if moebius(n) = 1
    // 6) lpf          if moebius(n) = -1

    if (n == 1)
      check(factorTable.is_leaf(i) == uint16_max - 1);
    else if (is_prime)
      check(factorTable.is_leaf(i) == uint16_max);
    else if (mu[n] == 0)
      check(factorTable.is_leaf(i) == 0);
    else
      check(lpf[n] == factorTable.is_leaf(i) + (factorTable.mu(i) == 1));

    not_coprime:;
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

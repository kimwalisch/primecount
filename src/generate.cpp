///
/// @file  generate.cpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <generate.hpp>

#include <stdint.h>
#include <limits>
#include <vector>

using namespace std;

namespace primecount {

/// Generate a vector with Möbius function values.
/// This implementation is based on code by Rick Sladkey:
/// https://mathoverflow.net/q/99545
///
vector<int32_t> generate_moebius(int64_t max)
{
  vector<int32_t> mu(max + 1, 1);

  for (int32_t i = 2; i * i <= max; i++)
  {
    if (mu[i] == 1)
    {
      for (int32_t j = i; j <= max; j += i)
        mu[j] *= -i;
      for (int32_t j = i * i; j <= max; j += i * i)
        mu[j] = 0;
    }
  }

  for (int32_t i = 2; i <= max; i++)
  {
    if (mu[i] == i)
      mu[i] = 1;
    else if (mu[i] == -i)
      mu[i] = -1;
    else if (mu[i] < 0)
      mu[i] = 1;
    else if (mu[i] > 0)
      mu[i] = -1;
  }

  return mu;
}

/// Generate a vector with the least prime factors
/// of the integers <= max
///
vector<int32_t> generate_lpf(int64_t max)
{
  vector<int32_t> lpf(max + 1, 1);

  // phi(x / 1, c) contributes to the sum in the
  // Lagarias-Miller-Odlyzko prime counting algorithm,
  // thus set lpf[1] = MAX (normally lpf[1] = 1)
  if (lpf.size() > 1)
    lpf[1] = numeric_limits<int32_t>::max();

  for (int32_t i = 2; i * i <= max; i++)
    if (lpf[i] == 1)
      for (int32_t j = i * 2; j <= max; j += i)
        if (lpf[j] == 1)
          lpf[j] = i;

  for (int32_t i = 2; i <= max; i++)
    if (lpf[i] == 1)
      lpf[i] = i;

  return lpf;
}

/// Generate a vector with the prime counts <= max
/// using the sieve of Eratosthenes
///
vector<int32_t> generate_pi(int64_t max)
{
  int64_t size = max + 1;
  vector<char> sieve(size, 1);

  for (int64_t i = 2; i * i <= max; i++)
    if (sieve[i])
      for (int64_t j = i * i; j <= max; j += i)
        sieve[j] = 0;

  vector<int32_t> pi(size, 0);
  int32_t pix = 0;

  for (int64_t x = 2; x <= max; x++)
  {
    pix += sieve[x];
    pi[x] = pix;
  }

  return pi;
}

} // namespace

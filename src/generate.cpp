///
/// @file  generate.cpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <generate.hpp>
#include <isqrt.hpp>

#include <stdint.h>
#include <limits>
#include <vector>

using namespace std;

namespace primecount {

/// Generate a vector with the prime counts <= max
/// using the sieve of Eratosthenes
///
vector<int32_t> generate_pi(int64_t max)
{
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  vector<char> sieve(size, 1);

  for (int64_t i = 2; i <= sqrt; i++)
    if (sieve[i])
      for (int64_t j = i * i; j < size; j += i)
        sieve[j] = 0;

  vector<int32_t> pi(size, 0);
  int32_t pix = 0;

  for (int64_t i = 2; i < size; i++)
  {
    pix += sieve[i];
    pi[i] = pix;
  }

  return pi;
}

/// Generate a vector with Möbius function values.
/// This implementation is based on code by Rick Sladkey:
/// https://mathoverflow.net/q/99545
///
vector<int32_t> generate_moebius(int64_t max)
{
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  vector<int32_t> mu(size, 1);

  for (int64_t i = 2; i <= sqrt; i++)
  {
    if (mu[i] == 1)
    {
      for (int64_t j = i; j < size; j += i)
        mu[j] *= (int32_t) -i;
      for (int64_t j = i * i; j < size; j += i * i)
        mu[j] = 0;
    }
  }

  for (int64_t i = 2; i < size; i++)
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
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  vector<int32_t> lpf(size, 1);

  for (int64_t i = 2; i <= sqrt; i++)
    if (lpf[i] == 1)
      for (int64_t j = i * i; j < size; j += i)
        if (lpf[j] == 1)
          lpf[j] = (int32_t) i;

  for (int64_t i = 2; i < size; i++)
    if (lpf[i] == 1)
      lpf[i] = (int32_t) i;

  // phi(x / 1, c) contributes to the sum in the
  // Lagarias-Miller-Odlyzko prime counting algorithm,
  // thus set lpf[1] = MAX (normally lpf[1] = 1)
  if (lpf.size() > 1)
    lpf[1] = numeric_limits<int32_t>::max();

  return lpf;
}

} // namespace

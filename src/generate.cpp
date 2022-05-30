///
/// @file  generate.cpp
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <generate.hpp>
#include <isqrt.hpp>
#include <pod_vector.hpp>

#include <stdint.h>
#include <limits>

using std::numeric_limits;

namespace primecount {

/// Generate a vector with the prime counts <= max
/// using the sieve of Eratosthenes
///
pod_vector<int32_t> generate_pi(int64_t max)
{
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  pod_vector<bool> sieve(size);
  std::fill(sieve.begin(), sieve.end(), 1);

  for (int64_t i = 2; i <= sqrt; i++)
    if (sieve[i])
      for (int64_t j = i * i; j < size; j += i)
        sieve[j] = 0;

  pod_vector<int32_t> pi(size);
  std::fill(pi.begin(), pi.end(), 0);
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
pod_vector<int32_t> generate_moebius(int64_t max)
{
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  pod_vector<int32_t> mu(size);
  std::fill(mu.begin(), mu.end(), 1);

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
/// of the integers <= max.
/// @Examples: lfp(2) = 2, lpf(15) = 3
///
pod_vector<int32_t> generate_lpf(int64_t max)
{
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  pod_vector<int32_t> lpf(size);
  std::fill(lpf.begin(), lpf.end(), 1);

  // By convention lfp(1) = +Infinity. Note that lpf(n) is
  // named pmin(n) in Tomás Oliveira e Silva's paper:
  // "Computing π(x): the combinatorial method".
  // The reason why lfp(1) is defined to be +Infinity is
  // that phi(x / 1, c) contributes to the ordinary leaves
  // (S1) in the Lagarias-Miller-Odlyzko and
  // Deleglise-Rivat prime counting algorithms. And
  // lfp(1) = +Infinity allows to simplify that algorithm.
  if (lpf.size() > 1)
    lpf[1] = numeric_limits<int32_t>::max();

  for (int64_t i = 2; i <= sqrt; i++)
    if (lpf[i] == 1)
      for (int64_t j = i * i; j < size; j += i)
        if (lpf[j] == 1)
          lpf[j] = (int32_t) i;

  for (int64_t i = 2; i < size; i++)
    if (lpf[i] == 1)
      lpf[i] = (int32_t) i;

  return lpf;
}

/// Generate a vector with the largest prime factors
/// of the integers <= max.
/// @Examples: mfp(2) = 2, mpf(15) = 5
///
pod_vector<int32_t> generate_mpf(int64_t max)
{
  int64_t size = max + 1;
  pod_vector<int32_t> mpf(size);
  std::fill(mpf.begin(), mpf.end(), 1);

  for (int64_t i = 2; i <= max; i++)
    if (mpf[i] == 1)
      for (int64_t j = i; j < size; j += i)
        mpf[j] = (int32_t) i;

  return mpf;
}

} // namespace

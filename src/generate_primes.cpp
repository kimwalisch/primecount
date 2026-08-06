///
/// @file  generate_primes.cpp
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <generate_primes.hpp>
#include <primecount.hpp>
#include <primecount-internal.hpp>
#include <int128_t.hpp>
#include <isqrt.hpp>
#include <macros.hpp>
#include <primesieve.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <stdint.h>

namespace primecount {

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
Vector<uint32_t> generate_primes_u32(int64_t max)
{
  Vector<uint32_t> primes;
  primes.resize(1);
  primes[0] = 0;
  primesieve::generate_primes(max, &primes);
  return primes;
}

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
Vector<int64_t> generate_primes_i64(int64_t max)
{
  Vector<int64_t> primes;
  primes.resize(1);
  primes[0] = 0;
  primesieve::generate_primes(max, &primes);
  return primes;
}

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
Vector<uint64_t> generate_primes_u64(int64_t max)
{
  Vector<uint64_t> primes;
  primes.resize(1);
  primes[0] = 0;
  primesieve::generate_primes(max, &primes);
  return primes;
}

/// Returns a vector with the first n primes.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
//
Vector<uint32_t> generate_n_primes_u32(int64_t n)
{
  Vector<uint32_t> primes;
  primes.reserve(n + 1);
  primes.push_back(0);
  primesieve::generate_n_primes(n, &primes);
  return primes;
}

/// Returns a vector with Möbius function values.
/// For each prime p <= max we flip the sign of the multiples
/// of p and set the multiples of p^2 to 0. Hence at the end
/// mu[n] is (-1)^omega(n) if n is square free and 0 otherwise.
///
Vector<int8_t> generate_moebius(int64_t max,
                                const Vector<uint32_t>& primes)
{
  ASSERT(pi_noprint(max, 1) == int64_t(primes.size() - 1));

  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  Vector<int8_t> mu(size);
  std::fill(mu.begin(), mu.end(), 1);

  for (std::size_t i = 1; i < primes.size(); i++)
  {
    int64_t prime = primes[i];

    for (int64_t j = prime; j < size; j += prime)
      mu[j] = -mu[j];

    if (prime <= sqrt)
      for (int64_t j = prime * prime; j < size; j += prime * prime)
        mu[j] = 0;
  }

  return mu;
}

/// Returns a vector with the largest prime factors
/// of the integers <= max.
/// @Examples: mfp(2) = 2, mpf(15) = 5
///
Vector<uint32_t> generate_mpf(int64_t max,
                              const Vector<uint32_t>& primes)
{
  if_unlikely(max > pstd::numeric_limits<uint32_t>::max())
    throw primecount_error("generate_mpf(max): max must be < 2^32");

  ASSERT(pi_noprint(max, 1) == int64_t(primes.size() - 1));

  int64_t size = max + 1;
  Vector<uint32_t> mpf(size);
  std::fill(mpf.begin(), mpf.end(), 1);

  for (std::size_t i = 1; i < primes.size(); i++)
    for (int64_t j = primes[i]; j < size; j += primes[i])
      mpf[j] = primes[i];

  return mpf;
}

/// Returns a vector with the least prime factors
/// of the integers <= max.
/// @Examples: lfp(2) = 2, lpf(15) = 3
///
Vector<uint32_t> generate_lpf(int64_t max,
                              const Vector<uint32_t>& primes)
{
  if_unlikely(max > pstd::numeric_limits<uint32_t>::max())
    throw primecount_error("generate_lpf(max): max must be < 2^32");

  ASSERT((max <= 1 && primes.size() <= 1) ||
         isqrt(max) < primes.back());

  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  Vector<uint32_t> lpf(size);
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
    lpf[1] = pstd::numeric_limits<uint32_t>::max();

  for (std::size_t i = 1; i < primes.size(); i++)
  {
    if (primes[i] > sqrt)
      break;

    int64_t prime_squared = primes[i] * primes[i];
    for (int64_t j = prime_squared; j < size; j += primes[i])
      if (lpf[j] == 1)
        lpf[j] = primes[i];
  }

  for (int64_t i = 2; i < size; i++)
    if (lpf[i] == 1)
      lpf[i] = uint32_t(i);

  return lpf;
}

} // namespace

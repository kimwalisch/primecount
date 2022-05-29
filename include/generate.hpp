///
/// @file  generate.hpp
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef GENERATE_HPP
#define GENERATE_HPP

#include <primesieve.hpp>
#include <pod_vector.hpp>

#include <stdint.h>

namespace primecount {

/// Generate a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
pod_vector<T> generate_primes(int64_t max)
{
  pod_vector<T> primes;
  primes.resize(1);
  primes[0] = 0;
  primesieve::generate_primes(max, &primes);
  return primes;
}

/// Generate a vector with the first n primes.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
//
template <typename T>
pod_vector<T> generate_n_primes(int64_t n)
{
  pod_vector<T> primes;
  primes.reserve(n + 1);
  primes.push_back(0);
  primesieve::generate_n_primes(n, &primes);
  return primes;
}

/// Generate a vector with Möbius function values
pod_vector<int32_t> generate_moebius(int64_t max);

/// Generate a vector with the least prime
/// factors of the integers <= max.
///
pod_vector<int32_t> generate_lpf(int64_t max);

/// Generate a vector with the largest prime
/// factors of the integers <= max.
///
pod_vector<int32_t> generate_mpf(int64_t max);

/// Generate a vector with the prime counts <= max
/// using the sieve of Eratosthenes.
///
pod_vector<int32_t> generate_pi(int64_t max);

} // namespace

#endif

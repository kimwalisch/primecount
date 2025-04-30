///
/// @file  generate_primes.hpp
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef GENERATE_PRIMES_HPP
#define GENERATE_PRIMES_HPP

#include <Vector.hpp>

#include <type_traits>
#include <stdint.h>

namespace primecount {

/// defined in generate_primes.cpp
Vector<int32_t> generate_primes_i32(int64_t max);
Vector<uint32_t> generate_primes_u32(int64_t max);
Vector<int64_t> generate_primes_i64(int64_t max);
Vector<uint64_t> generate_primes_u64(int64_t max);
Vector<int32_t> generate_n_primes_i32(int64_t n);

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
typename std::enable_if<std::is_same<T, int32_t>::value, Vector<int32_t>>::type
generate_primes(int64_t max)
{
  return generate_primes_i32(max);
}

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
typename std::enable_if<std::is_same<T, uint32_t>::value, Vector<uint32_t>>::type
generate_primes(int64_t max)
{
  return generate_primes_u32(max);
}

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
typename std::enable_if<std::is_same<T, int64_t>::value, Vector<int64_t>>::type
generate_primes(int64_t max)
{
  return generate_primes_i64(max);
}

/// Returns a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
template <typename T>
typename std::enable_if<std::is_same<T, uint64_t>::value, Vector<uint64_t>>::type
generate_primes(int64_t max)
{
  return generate_primes_u64(max);
}

/// Returns a vector with the first n primes.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
//
template <typename T>
typename std::enable_if<std::is_same<T, int32_t>::value, Vector<int32_t>>::type
generate_n_primes(int64_t n)
{
  return generate_n_primes_i32(n);
}

/// Returns a vector with Möbius function values
Vector<int32_t> generate_moebius(int64_t max);

/// Returns a vector with the least prime
/// factors of the integers <= max.
///
Vector<int32_t> generate_lpf(int64_t max);

/// Returns a vector with the largest prime
/// factors of the integers <= max.
///
Vector<int32_t> generate_mpf(int64_t max);

/// Returns a vector with the prime counts <= max
/// using the sieve of Eratosthenes.
///
Vector<int32_t> generate_pi(int64_t max);

} // namespace

#endif

///
/// @file   get_primes.cpp
/// @brief  Test the PiTable::get_primes() methods
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <PiTable.hpp>
#include <Vector.hpp>
#include <generate_primes.hpp>
#include <primecount.hpp>

#include <array>
#include <cstdlib>
#include <iostream>
#include <stdint.h>

using namespace primecount;

/// First 100 primes, generated independently of primecount.
const std::array<uint32_t, 100> primes_tiny =
{
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
  283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
  353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
  419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
  467, 479, 487, 491, 499, 503, 509, 521, 523, 541
};

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

template <typename T>
bool equals_tiny(const Vector<T>& primes,
                 std::size_t prime_count)
{
  if (primes.size() != prime_count + 1 || primes[0] != 0)
    return false;

  for (std::size_t i = 0; i < prime_count; i++)
    if (primes[i + 1] != primes_tiny[i])
      return false;

  return true;
}

std::size_t tiny_prime_count(uint64_t x)
{
  std::size_t count = 0;

  while (count < primes_tiny.size() && primes_tiny[count] <= x)
    count++;

  return count;
}

template <typename T>
bool equals(const Vector<T>& primes,
            const Vector<uint32_t>& expected)
{
  if (primes.size() != expected.size())
    return false;

  for (std::size_t i = 0; i < expected.size(); i++)
    if (primes[i] != expected[i])
      return false;

  return true;
}

int main()
{
  uint64_t max = 4000000;
  PiTable pi(max, 4);

  // First test against an independent table. This ensures these tests
  // still fail if PiTable and generate_primes() contain the same bug.
  const std::array<uint64_t, 12> tiny_limits =
  {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 539, 540, 541
  };

  for (uint64_t limit : tiny_limits)
  {
    std::size_t prime_count = tiny_prime_count(limit);

    auto primes_u32 = pi.get_primes<uint32_t>(limit, 1);
    std::cout << "get_primes<uint32_t>(" << limit << ").size() = " << primes_u32.size();
    check(equals_tiny(primes_u32, prime_count));

    auto primes_i64 = pi.get_primes<int64_t>(limit, 1);
    std::cout << "get_primes<int64_t>(" << limit << ").size() = " << primes_i64.size();
    check(equals_tiny(primes_i64, prime_count));
  }

  for (std::size_t n = 0; n <= primes_tiny.size(); n++)
  {
    auto first_n = pi.get_n_primes<uint32_t>(n);
    std::cout << "get_n_primes<uint32_t>(" << n << ").size() = " << first_n.size();
    check(equals_tiny(first_n, n));
  }

  // Each limit > 3 * PiTable's thread threshold, hence all requested
  // threads below are used by get_primes(). The limits use different
  // remainders modulo 240 in order to vary the bitmask that unsets the
  // bits > limit in the final thread's last 64-bit word.
  uint64_t boundary = max - max % 240;
  const std::array<uint64_t, 5> limits =
  {
    boundary - 1, boundary, boundary + 1, boundary + 120, max
  };
  const std::array<int, 4> thread_counts = { 1, 2, 3, 4 };

  for (uint64_t limit : limits)
  {
    auto expected = generate_primes<uint32_t>(limit);

    for (int threads : thread_counts)
    {
      auto primes_u32 = pi.get_primes<uint32_t>(limit, threads);
      std::cout << "get_primes<uint32_t>(" << limit << ", threads = " << threads
                << ").size() = " << primes_u32.size();
      check(equals(primes_u32, expected));

      auto primes_i64 = pi.get_primes<int64_t>(limit, threads);
      std::cout << "get_primes<int64_t>(" << limit << ", threads = " << threads
                << ").size() = " << primes_i64.size();
      check(equals(primes_i64, expected));
    }
  }

  auto expected = generate_primes<uint32_t>(max);
  auto first_n = pi.get_n_primes<uint32_t>(expected.size() - 1);
  std::cout << "get_n_primes<uint32_t>(" << expected.size() - 1 << ").size() = " << first_n.size();
  check(equals(first_n, expected));

  bool exception_thrown = false;
  try
  {
    pi.get_primes<uint32_t>(max + 1, 1);
  }
  catch (const primecount::primecount_error&)
  {
    exception_thrown = true;
  }
  std::cout << "get_primes<uint32_t>(" << max + 1 << ") throws primecount_error";
  check(exception_thrown);

  exception_thrown = false;
  try
  {
    pi.get_primes<int64_t>(max + 1, 1);
  }
  catch (const primecount::primecount_error&)
  {
    exception_thrown = true;
  }
  std::cout << "get_primes<int64_t>(" << max + 1 << ") throws primecount_error";
  check(exception_thrown);

  exception_thrown = false;
  try
  {
    pi.get_n_primes<uint32_t>(pi[max] + 1);
  }
  catch (const primecount::primecount_error&)
  {
    exception_thrown = true;
  }
  std::cout << "get_n_primes<uint32_t>(" << pi[max] + 1 << ") throws primecount_error";
  check(exception_thrown);

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

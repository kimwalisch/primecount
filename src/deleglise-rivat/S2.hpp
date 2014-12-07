///
/// @file  S2.hpp.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S2_HPP
#define S2_HPP

#include <PiTable.hpp>
#include <FactorTable.hpp>
#include <int128.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

int64_t S2_easy(int64_t x, int64_t y, int64_t z, int64_t c, std::vector<int32_t>& pi, std::vector<int32_t>& primes, int threads);
int64_t S2_easy(int64_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<int32_t>& primes, int threads);

#ifdef HAVE_INT128_T

int128_t S2_easy(uint128_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<uint32_t>& primes, int threads);
int128_t S2_easy(uint128_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<int64_t>& primes, int threads);

#endif

int64_t S2_sieve(int64_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<int32_t>& primes, FactorTable<uint16_t>& factors);
int64_t S2_sieve(int64_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<int32_t>& primes, FactorTable<uint16_t>& factors, int threads);

#ifdef HAVE_INT128_T

int128_t S2_sieve(uint128_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<uint32_t>& primes, FactorTable<uint16_t>& factors);
int128_t S2_sieve(uint128_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<int64_t>& primes, FactorTable<uint32_t>& factors);

int128_t S2_sieve(uint128_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<uint32_t>& primes, FactorTable<uint16_t>& factors, int threads);
int128_t S2_sieve(uint128_t x, int64_t y, int64_t z, int64_t c, PiTable& pi, std::vector<int64_t>& primes, FactorTable<uint32_t>& factors, int threads);

#endif

} // namespace primecount

#endif

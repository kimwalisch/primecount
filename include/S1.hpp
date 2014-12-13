///
/// @file  S1.hpp
/// @brief S1 function declarations.
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <FactorTable.hpp>
#include <inttypes.hpp>

#include <stdint.h>
#include <vector>

namespace primecount {

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int64_t prime_c,
           std::vector<int32_t>& lpf,
           std::vector<int32_t>& mu,
           int threads);

int64_t S1(int64_t x,
           int64_t y,
           int64_t c,
           int64_t prime_c,
           FactorTable<uint16_t>& factors,
           int threads);

#ifdef HAVE_INT128_T

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int64_t prime_c,
            std::vector<int32_t>& lpf,
            std::vector<int32_t>& mu,
            int threads);

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            uint32_t prime_c,
            FactorTable<uint16_t>& factors,
            int threads);

int128_t S1(int128_t x,
            int64_t y,
            int64_t c,
            int64_t prime_c,
            FactorTable<uint32_t>& factors,
            int threads);

#endif

} // namespace primecount

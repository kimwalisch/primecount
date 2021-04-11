///
/// @file  S.hpp
/// @brief The S1 and S2 functions are part of the Deleglise-Rivat
///        prime counting algorithm.
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef S_HPP
#define S_HPP

#include <int128_t.hpp>
#include <print.hpp>

#include <stdint.h>

namespace primecount {

int64_t S1(int64_t x, int64_t y, int64_t c, int threads, bool print = is_print());
int64_t S2_trivial(int64_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int64_t S2_easy(int64_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int64_t S2_hard(int64_t x, int64_t y, int64_t z, int64_t c, int64_t s2_hard_approx, int threads, bool print = is_print());

#ifdef HAVE_INT128_T

int128_t S1(int128_t x, int64_t y, int64_t c, int threads, bool print = is_print());
int128_t S2_trivial(int128_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int128_t S2_easy(int128_t x, int64_t y, int64_t z, int64_t c, int threads, bool print = is_print());
int128_t S2_hard(int128_t x, int64_t y, int64_t z, int64_t c, int128_t s2_hard_approx, int threads, bool print = is_print());

#endif

} // namespace

#endif

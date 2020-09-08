///
/// @file  print.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRINT_HPP
#define PRINT_HPP

#include <int128_t.hpp>
#include <noinline.hpp>

#include <stdint.h>
#include <string>

namespace primecount {

void set_print(bool print);
void set_print_variables(bool print_variables);

bool is_print();
bool is_print_combined_result();

NOINLINE void print(const std::string& str);
NOINLINE void print(const std::string& str, maxint_t res);
NOINLINE void print(const std::string& str, maxint_t res, double time);
NOINLINE void print(maxint_t x, int64_t y, int64_t z, int64_t c, int threads);
NOINLINE void print_vars(maxint_t x, int64_t y, int threads);
NOINLINE void print_vars(maxint_t x, int64_t y, int64_t c, int threads);
NOINLINE void print_seconds(double seconds);

NOINLINE void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads);
NOINLINE void print_gourdon_vars(maxint_t x, int64_t y, int threads);
NOINLINE void print_gourdon_vars(maxint_t x, int64_t y, int64_t z, int64_t k,  int threads);

} // namespace

#endif

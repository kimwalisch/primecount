///
/// @file  print.hpp
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRINT_HPP
#define PRINT_HPP

#include <int128_t.hpp>
#include <stdint.h>

#if !defined(__has_include)
  #define __has_include(x) 0
#endif

#if __cplusplus >= 201703L && __has_include(<string_view>)
  #include <string_view>
  using string_view_t = std::string_view;
#else
  using string_view_t = const char*;
#endif

namespace primecount {

void set_print(bool print);
void set_print_variables(bool print_variables);

bool is_print();
bool is_print_combined_result();

void print(string_view_t str);
void print(string_view_t str, maxint_t res);
void print(string_view_t str, maxint_t res, double time);
void print(maxint_t x, int64_t y, int64_t z, int64_t c, int threads);
void print_vars(maxint_t x, int64_t y, int threads);
void print_vars(maxint_t x, int64_t y, int64_t c, int threads);
void print_seconds(double seconds);

void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads);
void print_gourdon_vars(maxint_t x, int64_t y, int threads);
void print_gourdon_vars(maxint_t x, int64_t y, int64_t z, int64_t k,  int threads);

} // namespace

#endif

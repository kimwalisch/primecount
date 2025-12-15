///
/// @file  print.hpp
///
/// Copyright (C) 2025 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRINT_HPP
#define PRINT_HPP

#include <int128_t.hpp>
#include <macros.hpp>

#include <stdint.h>
#include <string>

#if __cplusplus >= 201703L && \
    __has_include(<string_view>)
  #include <string_view>
  using string_view_t = std::string_view;
#else
  using string_view_t = const char*;
#endif

namespace primecount {

std::string to_string(double x, int precision);

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
void print_status(const std::string& status);

void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads);
void print_gourdon_vars(maxint_t x, int64_t y, int threads);
void print_gourdon_vars(maxint_t x, int64_t y, int64_t z, int64_t k,  int threads);
void print_nth_prime_sieve(uint64_t n, bool sieve_forward, maxint_t nth_prime_approx, uint64_t dist_approx, uint64_t segment_size, int threads);

} // namespace

#endif

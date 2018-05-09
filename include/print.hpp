///
/// @file  print.hpp
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRINT_HPP
#define PRINT_HPP

#include <int128_t.hpp>

#include <stdint.h>
#include <string>

namespace primecount {

void set_print(bool print);

void set_print_variables(bool print_variables);

bool is_print();

bool print_result();

bool print_variables();

void print(const std::string& str);

void print(const std::string& str, maxint_t res);

void print(maxint_t x, int64_t y, int threads);

void print(maxint_t x, int64_t y, int64_t c, int threads);

void print(maxint_t x, int64_t y, int64_t z, int64_t c, double alpha, int threads);

void print(const std::string& str, maxint_t res, double time);

void print_seconds(double seconds);

} // namespace

#endif

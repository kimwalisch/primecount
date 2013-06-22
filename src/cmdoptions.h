///
/// @file  cmdoptions.h
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMECOUNT_CMDOPTIONS_H
#define PRIMECOUNT_CMDOPTIONS_H

#include <primecount.h>
#include <stdint.h>
#include <vector>

enum PrimeCountOptions {
  OPTION_HELP,
  OPTION_LEGENDRE,
  OPTION_LEHMER,
  OPTION_LI,
  OPTION_LIINV,
  OPTION_MEISSEL,
  OPTION_NTHPRIME,
  OPTION_NUMBER,
  OPTION_PHI,
  OPTION_PRIMESIEVE,
  OPTION_TEST,
  OPTION_THREADS,
  OPTION_VERSION
};

struct PrimeCountSettings
{
  std::vector<int64_t> n;
  int64_t option;
  int threads;
  PrimeCountSettings() :
    option(OPTION_LEHMER),
    threads(primecount::MAX_THREADS)
  { }
};

PrimeCountSettings processOptions(int, char**);

#endif

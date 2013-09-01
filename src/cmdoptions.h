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

namespace primecount {

enum OptionValues
{
  OPTION_HELP,
  OPTION_LEGENDRE,
  OPTION_LEHMER,
  OPTION_LMO,
  OPTION_LMOSIMPLE,
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

struct PrimeCountOptions
{
  int64_t x;
  int64_t a;
  int64_t option;
  int threads;
  PrimeCountOptions() :
    x(-1),
    a(-1),
    option(OPTION_LEHMER),
    threads(primecount::MAX_THREADS)
  { }
};

PrimeCountOptions parseOptions(int, char**);

} // namespace primecount

#endif

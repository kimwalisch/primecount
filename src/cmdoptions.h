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

struct PrimeCountSettings
{
  int64_t x;
  int64_t method;
  int threads;
  PrimeCountSettings() :
    x(-1), method(3), threads(primecount::MAX_THREADS)
  { }
};

PrimeCountSettings processOptions(int, char**);

#endif

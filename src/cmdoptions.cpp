///
/// @file   cmdoptions.cpp
/// @brief  Parse command-line options for the primecount console
///         (terminal) application.
///
/// Copyright (C) 2013 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.h"
#include "ExpressionParser.h"

#include <stdint.h>
#include <string>
#include <map>
#include <exception>
#include <cstdlib>
#include <cstddef>


namespace primecount {

void help();
void test();
void version();

}

using namespace std;

namespace {

/// e.g. id = "--threads", value = "4"
struct Option {
  string id;
  string value;
  template <typename T>
  T getValue() const
  {
    ExpressionParser<T> parser;
    T result = parser.eval(value);
    return result;
  }
};

/// Command-line options
map<string, PrimeCountOptions> cmdOptions;

void initCmdOptions()
{
  cmdOptions["-h"]           = OPTION_HELP;
  cmdOptions["--help"]       = OPTION_HELP;
  cmdOptions["-g"]           = OPTION_LEGENDRE;
  cmdOptions["--legendre"]   = OPTION_LEGENDRE;
  cmdOptions["-l"]           = OPTION_LEHMER;
  cmdOptions["--lehmer"]     = OPTION_LEHMER;
  cmdOptions["--Li"]         = OPTION_LI;
  cmdOptions["--Li_inverse"] = OPTION_LIINV;
  cmdOptions["-m"]           = OPTION_MEISSEL;
  cmdOptions["--meissel"]    = OPTION_MEISSEL;
  cmdOptions["-n"]           = OPTION_NTHPRIME;
  cmdOptions["--nth_prime"]  = OPTION_NTHPRIME;
  cmdOptions["--number"]     = OPTION_NUMBER;
  cmdOptions["--phi"]        = OPTION_PHI;
  cmdOptions["-p"]           = OPTION_PRIMESIEVE;
  cmdOptions["--primesieve"] = OPTION_PRIMESIEVE;
  cmdOptions["--test"]       = OPTION_TEST;
  cmdOptions["-t"]           = OPTION_THREADS;
  cmdOptions["--threads"]    = OPTION_THREADS;
  cmdOptions["-v"]           = OPTION_VERSION;
  cmdOptions["--version"]    = OPTION_VERSION;
}

/// e.g. "--threads=8" -> { id = "--threads", value = "8" }
Option makeOption(const string& str)
{
  Option option;
  size_t delimiter = str.find_first_of("=0123456789");

  if (delimiter == string::npos) {
    option.id = str;
  } else {
    option.id = str.substr(0, delimiter);
    option.value = str.substr(delimiter + (str.at(delimiter) == '=' ? 1 : 0));
  }
  if (option.id.empty() && !option.value.empty())
    option.id = "--number";
  if (cmdOptions.count(option.id) == 0)
    option.id = "--help";

  return option;
}

} // end namespace

PrimeCountSettings processOptions(int argc, char** argv)
{
  // skip program name in argv[0]
  argc--; argv++;
  PrimeCountSettings pcs;
  initCmdOptions();
  try {
    for (int i = 0; i < argc; i++) {
      Option option = makeOption(argv[i]);

      switch (cmdOptions[option.id]) {
        case OPTION_PRIMESIEVE: pcs.option = OPTION_PRIMESIEVE; break;
        case OPTION_LEGENDRE:   pcs.option = OPTION_LEGENDRE; break;
        case OPTION_LI:         pcs.option = OPTION_LI; break;
        case OPTION_LIINV:      pcs.option = OPTION_LIINV; break;
        case OPTION_PHI:        pcs.option = OPTION_PHI; break;
        case OPTION_LEHMER:     pcs.option = OPTION_LEHMER; break;
        case OPTION_MEISSEL:    pcs.option = OPTION_MEISSEL; break;
        case OPTION_NTHPRIME:   pcs.option = OPTION_NTHPRIME; break;
        case OPTION_NUMBER:     pcs.n.push_back( option.getValue<int64_t>() ); break;
        case OPTION_THREADS:    pcs.threads = option.getValue<int>(); break;
        case OPTION_HELP:       primecount::help(); break;
        case OPTION_TEST:       primecount::test(); break;
        case OPTION_VERSION:    primecount::version(); break;
      }
    }
  } catch (exception&) {
    primecount::help();
  }

  if (pcs.n.empty())
    primecount::help();

  // phi(x, a) takes two arguments
  if (pcs.option == OPTION_PHI && pcs.n.size() != 2)
    primecount::help();

  return pcs;
}

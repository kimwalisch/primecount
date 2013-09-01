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


using std::string;

namespace primecount {

void help();
void test();
void version();

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
std::map<string, OptionValues> optionMap;

void initOptionMap()
{
  optionMap["-h"]           = OPTION_HELP;
  optionMap["--help"]       = OPTION_HELP;
  optionMap["-g"]           = OPTION_LEGENDRE;
  optionMap["--legendre"]   = OPTION_LEGENDRE;
  optionMap["-l"]           = OPTION_LEHMER;
  optionMap["--lehmer"]     = OPTION_LEHMER;
  optionMap["--lmo"]        = OPTION_LMO;
  optionMap["--lmo_simple"] = OPTION_LMOSIMPLE;
  optionMap["--Li"]         = OPTION_LI;
  optionMap["--Li_inverse"] = OPTION_LIINV;
  optionMap["-m"]           = OPTION_MEISSEL;
  optionMap["--meissel"]    = OPTION_MEISSEL;
  optionMap["-n"]           = OPTION_NTHPRIME;
  optionMap["--nthprime"]   = OPTION_NTHPRIME;
  optionMap["--number"]     = OPTION_NUMBER;
  optionMap["-p"]           = OPTION_PRIMESIEVE;
  optionMap["--primesieve"] = OPTION_PRIMESIEVE;
  optionMap["--test"]       = OPTION_TEST;
  optionMap["-t"]           = OPTION_THREADS;
  optionMap["--threads"]    = OPTION_THREADS;
  optionMap["-v"]           = OPTION_VERSION;
  optionMap["--version"]    = OPTION_VERSION;
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
  if (optionMap.count(option.id) == 0)
    option.id = "--help";

  return option;
}

PrimeCountOptions parseOptions(int argc, char** argv)
{
  // skip program name in argv[0]
  argc--; argv++;
  PrimeCountOptions pco;
  initOptionMap();
  try {
    for (int i = 0; i < argc; i++)
    {
      Option option = makeOption(argv[i]);
      switch (optionMap[option.id])
      {
        case OPTION_NUMBER:     pco.x       = option.getValue<int64_t>(); break;
        case OPTION_THREADS:    pco.threads = option.getValue<int>();     break;
        case OPTION_HELP:       help(); break;
        case OPTION_TEST:       test(); break;
        case OPTION_VERSION:    version(); break;
        default:                pco.option = optionMap[option.id];
      }
    }
  } catch (std::exception&) {
    help();
  }
  if (pco.x < 0)
    help();

  return pco;
}

} // namespace primecount

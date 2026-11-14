///
/// @file   parallel_limit.cpp
/// @brief  Test parallel regions if creating worker threads
///         fails, e.g. due to resource limits or on WebAssembly
///         without thread support. The thread pool must then
///         run the regions with the workers it already has.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <threadpool/parallel.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace primecount;

void check(bool OK, const std::string& msg)
{
  std::cout << msg << "   " << (OK ? "OK" : "ERROR") << std::endl;
  if (!OK)
    std::exit(1);
}

/// Every id must run exactly once, also in nested regions
bool check_ids(int threads)
{
  Vector<int> counts(threads * threads);
  std::fill(counts.begin(), counts.end(), 0);

  parallel_region(threads, [&](int i) {
    parallel_region(threads, [&](int j) {
      counts[i * threads + j]++;
    });
  });

  return std::all_of(counts.begin(), counts.end(), [](int n) { return n == 1; });
}

bool check_sum(int threads)
{
  int64_t sum = parallel_sum<int64_t>(threads, [](int id) -> int64_t { return id; });
  return sum == (int64_t) threads * (threads - 1) / 2;
}

bool check_exception(int threads)
{
  try
  {
    parallel_region(threads, [](int id) {
      if (id == 1)
        throw std::runtime_error("parallel_region test exception");
    });
  }
  catch (const std::runtime_error&)
  {
    return true;
  }

  return false;
}

int main()
{
  int max_workers = max_threads() - 1;

  for (int limit : {0, 1, 3})
  {
    set_parallel_pool_limit(limit);
    bool OK = true;

    for (int threads : {2, 3, 8, 64})
    {
      OK = OK && check_ids(threads);
      OK = OK && check_sum(threads);
      OK = OK && check_exception(threads);
    }

    int workers = std::min(limit, max_workers);
    OK = OK && parallel_pool_size() == workers;
    check(OK, "Thread creation fails beyond " + std::to_string(workers) + " workers");
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

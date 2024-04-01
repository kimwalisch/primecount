///
/// @file   CPUID.cpp
/// @brief  Test CPUID code on x86 and x64 CPUs.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <CPUID.hpp>
#include <iostream>

int main()
{
#if defined(__x86_64__) || \
    defined(__i386__) || \
    defined(_M_X64) || \
    defined(_M_IX86)

  #if defined(__POPCNT__) && !defined(HAS_POPCNT)
    std::cerr << "Error: HAS_POPCNT must be defined if __POPCNT__ is defined!" << std::endl;
    return 1;
  #endif

  #if defined(__AVX__) && !defined(HAS_POPCNT)
    std::cerr << "Error: HAS_POPCNT must be defined if __AVX__ is defined!" << std::endl;
    return 1;
  #endif

  #if defined(__AVX2__) && !defined(HAS_POPCNT)
    std::cerr << "Error: HAS_POPCNT must be defined if __AVX2__ is defined!" << std::endl;
    return 1;
  #endif

  #if defined(HAS_POPCNT) && defined(ENABLE_CPUID_POPCNT)
    std::cerr << "Error: ENABLE_CPUID_POPCNT must not be defined if HAS_POPCNT is defined!" << std::endl;
    return 1;
  #endif

  #if !defined(HAS_POPCNT) && !defined(ENABLE_CPUID_POPCNT)
    std::cerr << "Error: ENABLE_CPUID_POPCNT must be defined if HAS_POPCNT is not defined!" << std::endl;
    return 1;
  #endif

  #if defined(ENABLE_CPUID_POPCNT)
    std::cout << "CPU supports POPCNT: " << (CPUID_POPCNT ? "yes" : "no") << std::endl;
  #endif

#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

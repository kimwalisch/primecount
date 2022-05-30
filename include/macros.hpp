///
/// @file  macros.hpp
///
/// Copyright (C) 2022 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MACROS_HPP
#define MACROS_HPP

#ifndef __has_attribute
  #define __has_attribute(x) 0
#endif

#ifndef __has_builtin
  #define __has_builtin(x) 0
#endif

#ifndef __has_cpp_attribute
  #define __has_cpp_attribute(x) 0
#endif

/// Unfortunately compilers cannot be trusted (especially GCC)
/// to inline performance critical functions. We must ensure
/// that e.g. pi[x] and segmentedPi[x] are inlined.
///
#if __has_attribute(always_inline)
  #define ALWAYS_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
  #define ALWAYS_INLINE __forceinline
#else
  #define ALWAYS_INLINE
#endif

#if __cplusplus >= 202002L && \
    __has_cpp_attribute(unlikely)
  #define if_unlikely(x) if (x) [[unlikely]]
#elif defined(__GNUC__) || \
      __has_builtin(__builtin_expect)
  #define if_unlikely(x) if (__builtin_expect(!!(x), 0))
#else
  #define if_unlikely(x) if (x)
#endif

#if __cplusplus >= 201703L && \
    __has_cpp_attribute(fallthrough)
  #define FALLTHROUGH [[fallthrough]]
#elif __has_attribute(fallthrough)
  #define FALLTHROUGH __attribute__((fallthrough))
#else
  #define FALLTHROUGH
#endif

#if __cplusplus >= 201703L && \
    __has_cpp_attribute(maybe_unused)
  #define MAYBE_UNUSED [[maybe_unused]]
#elif __has_attribute(unused)
  #define MAYBE_UNUSED __attribute__((unused))
#else
  #define MAYBE_UNUSED
#endif

// Silence GCC < 12 warning:
// warning: 'unused' attribute ignored [-Wattributes]
#if defined(__GNUC__) && \
   !defined(__clang__)
  #if __GNUC__ < 12
    #undef MAYBE_UNUSED
    #define MAYBE_UNUSED
  #endif
#endif

#if defined(__GNUC__) || \
    __has_builtin(__builtin_unreachable)
  #define UNREACHABLE __builtin_unreachable()
#elif defined(_MSC_VER)
  #define UNREACHABLE __assume(0)
#else
  #define UNREACHABLE
#endif

#endif

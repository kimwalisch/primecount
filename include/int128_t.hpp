///
/// @file   int128_t.hpp
/// @brief  Defines int128_t and uint128_t integer types and adds
///         functions for 128-bit integers that are missing
///         in the C++ standard library.
///
/// Copyright (C) 2024 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_T_HPP
#define INT128_T_HPP

#include <stdint.h>
#include <string>
#include <limits>
#include <type_traits>

/// If INT128_MAX is defined we know that int128_t and
/// uint128_t are available in <stdint.h>.
///
#if defined(INT128_MAX) && \
   !defined(DISABLE_INT128)

#define HAVE_INT128_T

namespace primecount {

using maxint_t = int128_t ;
using maxuint_t = uint128_t;
using std::to_string;

} // namespace

#elif defined(__SIZEOF_INT128__) && \
     !defined(DISABLE_INT128)

#define HAVE_INT128_T
#define ENABLE_INT128_TO_STRING

#include <ostream>

namespace primecount {

using int128_t = __int128_t ;
using uint128_t = __uint128_t;
using maxint_t = __int128_t ;
using maxuint_t = __uint128_t;

/// These functions are defined in util.cpp
std::string to_string(int128_t x);
std::string to_string(uint128_t x);

std::ostream& operator<<(std::ostream& stream, int128_t n);
std::ostream& operator<<(std::ostream& stream, uint128_t n);

} // namespace

#elif __has_include(<__msvc_int128.hpp>) && \
     !defined(DISABLE_INT128)

#define HAVE_INT128_T
#define ENABLE_INT128_TO_STRING

// Experimental MSVC int128_t support:
// https://github.com/microsoft/STL/blob/main/stl/inc/__msvc_int128.hpp
// https://developercommunity.visualstudio.com/t/support-for-128-bit-integer-type/879048
// https://stackoverflow.com/a/76440171
#include <__msvc_int128.hpp>
#include <ostream>

namespace primecount {

using int128_t = std::_Signed128;
using uint128_t = std::_Unsigned128;
using maxint_t = std::_Signed128;
using maxuint_t = std::_Unsigned128;

/// These functions are defined in util.cpp
std::string to_string(int128_t x);
std::string to_string(uint128_t x);

std::ostream& operator<<(std::ostream& stream, int128_t n);
std::ostream& operator<<(std::ostream& stream, uint128_t n);

} // namespace

#else // int128_t not supported

namespace primecount {

using maxint_t = int64_t ;
using maxuint_t = uint64_t;
using std::to_string;

} // namespace

#endif

// Portable C++ type traits that support int128_t and uint128_t.
// This is required for GCC/Clang if the user compiles with -std=c++*
// instead of -std=gnu++* and also for LLVM/Clang on Windows.
namespace {
namespace port {

using namespace primecount;

// port::is_same
template<class T, class U>
struct is_same : std::false_type {};

template<class T>
struct is_same<T, T> : std::true_type {};

// port::conditional
template <bool Cond, class T, class F>
struct conditional {
  using type = T;
};

template <class T, class F>
struct conditional<false, T, F> {
  using type = F;
};

// port::is_integral
template<typename T> struct is_integral {
  static constexpr bool value = std::is_integral<T>::value;
};

#if defined(HAVE_INT128_T)
  template<> struct is_integral<int128_t> : std::true_type {};
  template<> struct is_integral<uint128_t> : std::true_type {};
#endif

// port::is_floating_point
template<typename T> struct is_floating_point {
  static constexpr bool value = std::is_floating_point<T>::value;
};

#if defined(HAVE_INT128_T)
  template<> struct is_floating_point<int128_t> : std::false_type {};
  template<> struct is_floating_point<uint128_t> : std::false_type {};
#endif

// port::is_signed
template<typename T> struct is_signed {
  static constexpr bool value = std::is_signed<T>::value;
};

#if defined(HAVE_INT128_T)
  template<> struct is_signed<int128_t> : std::true_type {};
  template<> struct is_signed<uint128_t> : std::false_type {};
#endif

// port::is_unsigned
template<typename T> struct is_unsigned {
  static constexpr bool value = std::is_unsigned<T>::value;
};

#if defined(HAVE_INT128_T)
  template<> struct is_unsigned<int128_t> : std::false_type {};
  template<> struct is_unsigned<uint128_t> : std::true_type {};
#endif

// port::make_unsigned
template<typename T> struct make_unsigned {
  using type = typename std::make_unsigned<T>::type;
};

#if defined(HAVE_INT128_T)
  template<> struct make_unsigned<int128_t> { using type = uint128_t; };
  template<> struct make_unsigned<uint128_t> { using type = uint128_t; };
#endif

// port::numeric_limits
template<typename T> struct numeric_limits {
    static constexpr T min() { return std::numeric_limits<T>::min(); }
    static constexpr T max() { return std::numeric_limits<T>::max(); }
    static constexpr T infinity() { return std::numeric_limits<T>::infinity(); }
    static constexpr T epsilon() { return std::numeric_limits<T>::epsilon(); }
    static constexpr int digits = std::numeric_limits<T>::digits;
};

#if defined(HAVE_INT128_T)
  template<> struct numeric_limits<int128_t>
  {
    static constexpr int128_t min() { return int128_t(uint128_t(1) << 127); }
    static constexpr int128_t max() { return int128_t((uint128_t(1) << 127) - 1); }
    static constexpr int digits = 127;
  };

  template<> struct numeric_limits<uint128_t>
  {
    static constexpr uint128_t min() { return 0; }
    static constexpr uint128_t max() { return ~uint128_t(0); }
    static constexpr int digits = 128;
  };
#endif

} // namespace port
} // namespace

#if defined(HAVE_INT128_T) && \
    defined(_OPENMP) && \
    defined(ENABLE_INT128_OPENMP_PATCH)
  #include <int128_OpenMP_patch.hpp>
#endif

#endif

///
/// @file  min_max.hpp
/// @brief In order to make the code more readable primecount allows
///        to use the min and max functions with different types
///        without casting if both types are integral.
///
///        int64_t a = 100;
///        int32_t b = 999;
///
///        primecount::max(a, b); // compiles
///        primecount::max(b, a); // does not compile
///
///        std::max(a, b); // does not compile
///        std::max(b, a); // does not compile
///
/// Copyright (C) 2014 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MIN_MAX_HPP
#define MIN_MAX_HPP

#if __cplusplus >= 201103L

#include <inttypes.hpp>

#include <algorithm>
#include <cassert>
#include <stdint.h>
#include <type_traits>

namespace primecount {

template <typename T>
struct make_signed_workaround
{
#ifndef HAVE_INT128_T
  typedef typename std::make_signed<T>::type type;
#else
  typedef typename std::conditional<std::is_same<T, uint8_t>::value, int8_t,
          typename std::conditional<std::is_same<T, uint16_t>::value, int16_t,
          typename std::conditional<std::is_same<T, uint32_t>::value, int32_t,
          typename std::conditional<std::is_same<T, uint64_t>::value, int64_t,
          typename std::conditional<std::is_same<T, uint128_t>::value, int128_t,
          T>::type>::type>::type>::type>::type type;
#endif
};

template <typename T>
struct is_integral_workaround
{
  enum {
#ifndef HAVE_INT128_T
    value = std::is_integral<T>::value
#else
    value = std::is_integral<T>::value ||
            std::is_same<T, int128_t>::value ||
            std::is_same<T, uint128_t>::value
#endif
  };
};

template <typename A, typename B>
struct integral_not_same
{
  enum {
    value = !std::is_same<A, B>::value &&
             is_integral_workaround<A>::value &&
             is_integral_workaround<B>::value
  };
};

template <typename A, typename B,
          typename std::enable_if<std::is_same<A, B>::value>::type* = nullptr>
B min(A a, B b)
{
  return std::min(a, b);
}

template <typename A, typename B,
          typename std::enable_if<std::is_same<A, B>::value>::type* = nullptr>
inline B min3(A a, B b, B c)
{
  return std::min(a, std::min(b, c));
}

template <typename A, typename B,
          typename std::enable_if<std::is_same<A, B>::value>::type* = nullptr>
inline A max(A a, B b)
{
  return std::max(a, b);
}

template <typename A, typename B,
          typename std::enable_if<std::is_same<A, B>::value>::type* = nullptr>
inline A max3(A a, B b, B c)
{
  return std::max(a, std::max(b, c));
}

/// Convenience min function for different integer types.
/// @pre sizeof(A) >= sizeof(B) &&
///      a >= 0 && b >= 0.
///
template <typename A, typename B,
          typename std::enable_if<integral_not_same<A, B>::value>::type* = nullptr>
B min(A a, B b)
{
  static_assert(sizeof(A) >= sizeof(B), "Cannot compare types A and B");
  assert(a >= 0);
  assert(b >= 0);

  return (B) std::min(a, (A) b);
}

/// Convenience min function for different integer types.
/// @pre sizeof(A) >= sizeof(B) &&
///      a >= 0 && b >= 0 && c >= 0.
///
template <typename A, typename B,
          typename std::enable_if<integral_not_same<A, B>::value>::type* = nullptr>
inline B min3(A a, B b, B c)
{
  static_assert(sizeof(A) >= sizeof(B), "Cannot compare types A and B");
  assert(a >= 0);
  assert(b >= 0);
  assert(c >= 0);

  return (B) std::min(a, (A) std::min(b, c));
}

/// Convenience max function for different integer types.
/// @pre sizeof(A) >= sizeof(B) &&
///      a >= 0 && b >= 0.
///
template <typename A, typename B,
          typename std::enable_if<integral_not_same<A, B>::value>::type* = nullptr>
inline A max(A a, B b)
{
  static_assert(sizeof(A) >= sizeof(B), "Cannot compare types A and B");
  assert(a >= 0);
  assert(b >= 0);

  return std::max(a, (A) b);
}

/// Convenience max function for different integer types.
/// @pre sizeof(A) >= sizeof(B) &&
///      a >= 0 && b >= 0 && c >= 0.
///
template <typename A, typename B,
          typename std::enable_if<integral_not_same<A, B>::value>::type* = nullptr>
inline A max3(A a, B b, B c)
{
  static_assert(sizeof(A) >= sizeof(B), "Cannot compare types A and B");
  assert(a >= 0);
  assert(b >= 0);
  assert(c >= 0);

  return std::max(a, (A) std::max(b, c));
}

} // namespace

#else /* C++98 */

#include <algorithm>
#include <cassert>

namespace primecount {

template <typename A, typename B>
inline B min(A a, B b)
{
  assert(sizeof(A) >= sizeof(B));
  return (B) std::min(a, (A) b);
}

template <typename A, typename B>
inline B min3(A a, B b, B c)
{
  assert(sizeof(A) >= sizeof(B));
  return (B) std::min(a, (A) std::min(b, c));
}

template <typename A, typename B>
inline A max(A a, B b)
{
  assert(sizeof(A) >= sizeof(B));
  return std::max(a, (A) b);
}

template <typename A, typename B>
inline A max3(A a, B b, B c)
{
  assert(sizeof(A) >= sizeof(B));
  return std::max(a, (A) std::max(b, c));
}

} // namespace

#endif

#endif // MIN_MAX_HPP

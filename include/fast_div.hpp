///
/// @file  fast_div.hpp
/// @brief Integer division of small types is much faster than
///        integer division of large types on most CPUs. The
///        fast_div(x, y) function tries to take advantage of
///        this by casting x and y to smaller types (if possible)
///        before doing the division.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <int128_t.hpp>

#include <cassert>
#include <limits>
#include <type_traits>

namespace primecount {

template <typename T>
struct fastdiv
{
  typedef typename std::conditional<sizeof(T) / 2 <= sizeof(uint32_t), uint32_t,
          typename std::conditional<sizeof(T) / 2 <= sizeof(uint64_t), uint64_t,
          T>::type>::type type;
};

template <typename X, typename Y>
typename std::enable_if<(sizeof(X) == sizeof(Y)), X>::type
fast_div(X x, Y y)
{
  static_assert(prt::is_integral<X>::value &&
                prt::is_integral<Y>::value,
                "fast_div(x, y): types must be integral");

  using fastdiv_t = typename fastdiv<X>::type;

  if (x <= std::numeric_limits<fastdiv_t>::max() &&
      y <= std::numeric_limits<fastdiv_t>::max())
  {
    return (fastdiv_t) x / (fastdiv_t) y;
  }

  return x / y;
}

template <typename X, typename Y>
typename std::enable_if<(sizeof(X) > sizeof(Y)), X>::type
fast_div(X x, Y y)
{
  static_assert(prt::is_integral<X>::value &&
                prt::is_integral<Y>::value,
                "fast_div(x, y): types must be integral");

  using fastdiv_t = typename fastdiv<X>::type;

  if (x <= std::numeric_limits<fastdiv_t>::max())
    return (fastdiv_t) x / (fastdiv_t) y;

  return x / y;
}

template <typename R, typename X, typename Y>
typename std::enable_if<(sizeof(R) == sizeof(uint64_t) &&
                         sizeof(X) == sizeof(uint128_t) &&
                         sizeof(Y) <= sizeof(uint64_t)), R>::type
fast_div(X x, Y y)
{
#if defined(__x86_64__) && \
   (defined(__GNUC__) || defined(__clang__))

  // primecount does not need signed division so 
  // we use the unsigned division instruction further
  // down as DIV is usually faster than IDIV.
  assert(x >= 0 && y >= 0);

  uint64_t result;
  uint64_t remainder;
  uint64_t x0 = (uint64_t) x;
  uint64_t x1 = (uint64_t) (x >> 64);
  uint64_t d = y;

  // We know the result is 64-bit (even though the
  // numerator is 128-bit) so we can use the divq
  // instruction instead of doing a full 128-bit division.
  __asm__("divq %[d]"
          : "=a"(result), "=d"(remainder)
          : [d] "r"(d), "a"(x0), "d"(x1)
          );

  return result;
#else
  return (uint64_t) fast_div(x, y);
#endif
}

template <typename R, typename X, typename Y>
typename std::enable_if<!(sizeof(R) == sizeof(uint64_t) &&
                          sizeof(X) == sizeof(uint128_t) &&
                          sizeof(Y) <= sizeof(uint64_t)), R>::type
fast_div(X x, Y y)
{
  return fast_div(x, y);
}

} // namespace

#endif

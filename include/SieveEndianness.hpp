///
/// @file  SieveEndianness.hpp
/// @brief For performance reasons we cast a byte array to an
///        uint64_t array in Sieve.cpp. This is not endian-safe,
///        hence we need different lookup tables for little
///        and big-endian CPUs.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SIEVEENDIANNESS_HPP
#define SIEVEENDIANNESS_HPP

#include <stdint.h>
#include <array>

namespace {

#if defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__)
  #undef IS_BIG_ENDIAN
  #if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    #define IS_BIG_ENDIAN
  #endif
#endif

#if defined(IS_BIG_ENDIAN)

/// The values of the big-endian lookup tables have been computed by
/// converting the values of the little-endian lookup tables to
/// big-endian using __builtin_bswap64(unset_smaller[i]).

/// Unset bits < start
const std::array<uint64_t, 240> unset_smaller =
{
  0xffffffffffffffffull, 0xffffffffffffffffull, 0xfeffffffffffffffull, 0xfeffffffffffffffull,
  0xfeffffffffffffffull, 0xfeffffffffffffffull, 0xfeffffffffffffffull, 0xfeffffffffffffffull,
  0xfcffffffffffffffull, 0xfcffffffffffffffull, 0xfcffffffffffffffull, 0xfcffffffffffffffull,
  0xf8ffffffffffffffull, 0xf8ffffffffffffffull, 0xf0ffffffffffffffull, 0xf0ffffffffffffffull,
  0xf0ffffffffffffffull, 0xf0ffffffffffffffull, 0xe0ffffffffffffffull, 0xe0ffffffffffffffull,
  0xc0ffffffffffffffull, 0xc0ffffffffffffffull, 0xc0ffffffffffffffull, 0xc0ffffffffffffffull,
  0x80ffffffffffffffull, 0x80ffffffffffffffull, 0x80ffffffffffffffull, 0x80ffffffffffffffull,
  0x80ffffffffffffffull, 0x80ffffffffffffffull, 0x00ffffffffffffffull, 0x00ffffffffffffffull,
  0x00feffffffffffffull, 0x00feffffffffffffull, 0x00feffffffffffffull, 0x00feffffffffffffull,
  0x00feffffffffffffull, 0x00feffffffffffffull, 0x00fcffffffffffffull, 0x00fcffffffffffffull,
  0x00fcffffffffffffull, 0x00fcffffffffffffull, 0x00f8ffffffffffffull, 0x00f8ffffffffffffull,
  0x00f0ffffffffffffull, 0x00f0ffffffffffffull, 0x00f0ffffffffffffull, 0x00f0ffffffffffffull,
  0x00e0ffffffffffffull, 0x00e0ffffffffffffull, 0x00c0ffffffffffffull, 0x00c0ffffffffffffull,
  0x00c0ffffffffffffull, 0x00c0ffffffffffffull, 0x0080ffffffffffffull, 0x0080ffffffffffffull,
  0x0080ffffffffffffull, 0x0080ffffffffffffull, 0x0080ffffffffffffull, 0x0080ffffffffffffull,
  0x0000ffffffffffffull, 0x0000ffffffffffffull, 0x0000feffffffffffull, 0x0000feffffffffffull,
  0x0000feffffffffffull, 0x0000feffffffffffull, 0x0000feffffffffffull, 0x0000feffffffffffull,
  0x0000fcffffffffffull, 0x0000fcffffffffffull, 0x0000fcffffffffffull, 0x0000fcffffffffffull,
  0x0000f8ffffffffffull, 0x0000f8ffffffffffull, 0x0000f0ffffffffffull, 0x0000f0ffffffffffull,
  0x0000f0ffffffffffull, 0x0000f0ffffffffffull, 0x0000e0ffffffffffull, 0x0000e0ffffffffffull,
  0x0000c0ffffffffffull, 0x0000c0ffffffffffull, 0x0000c0ffffffffffull, 0x0000c0ffffffffffull,
  0x000080ffffffffffull, 0x000080ffffffffffull, 0x000080ffffffffffull, 0x000080ffffffffffull,
  0x000080ffffffffffull, 0x000080ffffffffffull, 0x000000ffffffffffull, 0x000000ffffffffffull,
  0x000000feffffffffull, 0x000000feffffffffull, 0x000000feffffffffull, 0x000000feffffffffull,
  0x000000feffffffffull, 0x000000feffffffffull, 0x000000fcffffffffull, 0x000000fcffffffffull,
  0x000000fcffffffffull, 0x000000fcffffffffull, 0x000000f8ffffffffull, 0x000000f8ffffffffull,
  0x000000f0ffffffffull, 0x000000f0ffffffffull, 0x000000f0ffffffffull, 0x000000f0ffffffffull,
  0x000000e0ffffffffull, 0x000000e0ffffffffull, 0x000000c0ffffffffull, 0x000000c0ffffffffull,
  0x000000c0ffffffffull, 0x000000c0ffffffffull, 0x00000080ffffffffull, 0x00000080ffffffffull,
  0x00000080ffffffffull, 0x00000080ffffffffull, 0x00000080ffffffffull, 0x00000080ffffffffull,
  0x00000000ffffffffull, 0x00000000ffffffffull, 0x00000000feffffffull, 0x00000000feffffffull,
  0x00000000feffffffull, 0x00000000feffffffull, 0x00000000feffffffull, 0x00000000feffffffull,
  0x00000000fcffffffull, 0x00000000fcffffffull, 0x00000000fcffffffull, 0x00000000fcffffffull,
  0x00000000f8ffffffull, 0x00000000f8ffffffull, 0x00000000f0ffffffull, 0x00000000f0ffffffull,
  0x00000000f0ffffffull, 0x00000000f0ffffffull, 0x00000000e0ffffffull, 0x00000000e0ffffffull,
  0x00000000c0ffffffull, 0x00000000c0ffffffull, 0x00000000c0ffffffull, 0x00000000c0ffffffull,
  0x0000000080ffffffull, 0x0000000080ffffffull, 0x0000000080ffffffull, 0x0000000080ffffffull,
  0x0000000080ffffffull, 0x0000000080ffffffull, 0x0000000000ffffffull, 0x0000000000ffffffull,
  0x0000000000feffffull, 0x0000000000feffffull, 0x0000000000feffffull, 0x0000000000feffffull,
  0x0000000000feffffull, 0x0000000000feffffull, 0x0000000000fcffffull, 0x0000000000fcffffull,
  0x0000000000fcffffull, 0x0000000000fcffffull, 0x0000000000f8ffffull, 0x0000000000f8ffffull,
  0x0000000000f0ffffull, 0x0000000000f0ffffull, 0x0000000000f0ffffull, 0x0000000000f0ffffull,
  0x0000000000e0ffffull, 0x0000000000e0ffffull, 0x0000000000c0ffffull, 0x0000000000c0ffffull,
  0x0000000000c0ffffull, 0x0000000000c0ffffull, 0x000000000080ffffull, 0x000000000080ffffull,
  0x000000000080ffffull, 0x000000000080ffffull, 0x000000000080ffffull, 0x000000000080ffffull,
  0x000000000000ffffull, 0x000000000000ffffull, 0x000000000000feffull, 0x000000000000feffull,
  0x000000000000feffull, 0x000000000000feffull, 0x000000000000feffull, 0x000000000000feffull,
  0x000000000000fcffull, 0x000000000000fcffull, 0x000000000000fcffull, 0x000000000000fcffull,
  0x000000000000f8ffull, 0x000000000000f8ffull, 0x000000000000f0ffull, 0x000000000000f0ffull,
  0x000000000000f0ffull, 0x000000000000f0ffull, 0x000000000000e0ffull, 0x000000000000e0ffull,
  0x000000000000c0ffull, 0x000000000000c0ffull, 0x000000000000c0ffull, 0x000000000000c0ffull,
  0x00000000000080ffull, 0x00000000000080ffull, 0x00000000000080ffull, 0x00000000000080ffull,
  0x00000000000080ffull, 0x00000000000080ffull, 0x00000000000000ffull, 0x00000000000000ffull,
  0x00000000000000feull, 0x00000000000000feull, 0x00000000000000feull, 0x00000000000000feull,
  0x00000000000000feull, 0x00000000000000feull, 0x00000000000000fcull, 0x00000000000000fcull,
  0x00000000000000fcull, 0x00000000000000fcull, 0x00000000000000f8ull, 0x00000000000000f8ull,
  0x00000000000000f0ull, 0x00000000000000f0ull, 0x00000000000000f0ull, 0x00000000000000f0ull,
  0x00000000000000e0ull, 0x00000000000000e0ull, 0x00000000000000c0ull, 0x00000000000000c0ull,
  0x00000000000000c0ull, 0x00000000000000c0ull, 0x0000000000000080ull, 0x0000000000000080ull,
  0x0000000000000080ull, 0x0000000000000080ull, 0x0000000000000080ull, 0x0000000000000080ull
};

/// Unset bits > stop
const std::array<uint64_t, 240> unset_larger =
{
  0x0000000000000000ull, 0x0100000000000000ull, 0x0100000000000000ull, 0x0100000000000000ull,
  0x0100000000000000ull, 0x0100000000000000ull, 0x0100000000000000ull, 0x0300000000000000ull,
  0x0300000000000000ull, 0x0300000000000000ull, 0x0300000000000000ull, 0x0700000000000000ull,
  0x0700000000000000ull, 0x0f00000000000000ull, 0x0f00000000000000ull, 0x0f00000000000000ull,
  0x0f00000000000000ull, 0x1f00000000000000ull, 0x1f00000000000000ull, 0x3f00000000000000ull,
  0x3f00000000000000ull, 0x3f00000000000000ull, 0x3f00000000000000ull, 0x7f00000000000000ull,
  0x7f00000000000000ull, 0x7f00000000000000ull, 0x7f00000000000000ull, 0x7f00000000000000ull,
  0x7f00000000000000ull, 0xff00000000000000ull, 0xff00000000000000ull, 0xff01000000000000ull,
  0xff01000000000000ull, 0xff01000000000000ull, 0xff01000000000000ull, 0xff01000000000000ull,
  0xff01000000000000ull, 0xff03000000000000ull, 0xff03000000000000ull, 0xff03000000000000ull,
  0xff03000000000000ull, 0xff07000000000000ull, 0xff07000000000000ull, 0xff0f000000000000ull,
  0xff0f000000000000ull, 0xff0f000000000000ull, 0xff0f000000000000ull, 0xff1f000000000000ull,
  0xff1f000000000000ull, 0xff3f000000000000ull, 0xff3f000000000000ull, 0xff3f000000000000ull,
  0xff3f000000000000ull, 0xff7f000000000000ull, 0xff7f000000000000ull, 0xff7f000000000000ull,
  0xff7f000000000000ull, 0xff7f000000000000ull, 0xff7f000000000000ull, 0xffff000000000000ull,
  0xffff000000000000ull, 0xffff010000000000ull, 0xffff010000000000ull, 0xffff010000000000ull,
  0xffff010000000000ull, 0xffff010000000000ull, 0xffff010000000000ull, 0xffff030000000000ull,
  0xffff030000000000ull, 0xffff030000000000ull, 0xffff030000000000ull, 0xffff070000000000ull,
  0xffff070000000000ull, 0xffff0f0000000000ull, 0xffff0f0000000000ull, 0xffff0f0000000000ull,
  0xffff0f0000000000ull, 0xffff1f0000000000ull, 0xffff1f0000000000ull, 0xffff3f0000000000ull,
  0xffff3f0000000000ull, 0xffff3f0000000000ull, 0xffff3f0000000000ull, 0xffff7f0000000000ull,
  0xffff7f0000000000ull, 0xffff7f0000000000ull, 0xffff7f0000000000ull, 0xffff7f0000000000ull,
  0xffff7f0000000000ull, 0xffffff0000000000ull, 0xffffff0000000000ull, 0xffffff0100000000ull,
  0xffffff0100000000ull, 0xffffff0100000000ull, 0xffffff0100000000ull, 0xffffff0100000000ull,
  0xffffff0100000000ull, 0xffffff0300000000ull, 0xffffff0300000000ull, 0xffffff0300000000ull,
  0xffffff0300000000ull, 0xffffff0700000000ull, 0xffffff0700000000ull, 0xffffff0f00000000ull,
  0xffffff0f00000000ull, 0xffffff0f00000000ull, 0xffffff0f00000000ull, 0xffffff1f00000000ull,
  0xffffff1f00000000ull, 0xffffff3f00000000ull, 0xffffff3f00000000ull, 0xffffff3f00000000ull,
  0xffffff3f00000000ull, 0xffffff7f00000000ull, 0xffffff7f00000000ull, 0xffffff7f00000000ull,
  0xffffff7f00000000ull, 0xffffff7f00000000ull, 0xffffff7f00000000ull, 0xffffffff00000000ull,
  0xffffffff00000000ull, 0xffffffff01000000ull, 0xffffffff01000000ull, 0xffffffff01000000ull,
  0xffffffff01000000ull, 0xffffffff01000000ull, 0xffffffff01000000ull, 0xffffffff03000000ull,
  0xffffffff03000000ull, 0xffffffff03000000ull, 0xffffffff03000000ull, 0xffffffff07000000ull,
  0xffffffff07000000ull, 0xffffffff0f000000ull, 0xffffffff0f000000ull, 0xffffffff0f000000ull,
  0xffffffff0f000000ull, 0xffffffff1f000000ull, 0xffffffff1f000000ull, 0xffffffff3f000000ull,
  0xffffffff3f000000ull, 0xffffffff3f000000ull, 0xffffffff3f000000ull, 0xffffffff7f000000ull,
  0xffffffff7f000000ull, 0xffffffff7f000000ull, 0xffffffff7f000000ull, 0xffffffff7f000000ull,
  0xffffffff7f000000ull, 0xffffffffff000000ull, 0xffffffffff000000ull, 0xffffffffff010000ull,
  0xffffffffff010000ull, 0xffffffffff010000ull, 0xffffffffff010000ull, 0xffffffffff010000ull,
  0xffffffffff010000ull, 0xffffffffff030000ull, 0xffffffffff030000ull, 0xffffffffff030000ull,
  0xffffffffff030000ull, 0xffffffffff070000ull, 0xffffffffff070000ull, 0xffffffffff0f0000ull,
  0xffffffffff0f0000ull, 0xffffffffff0f0000ull, 0xffffffffff0f0000ull, 0xffffffffff1f0000ull,
  0xffffffffff1f0000ull, 0xffffffffff3f0000ull, 0xffffffffff3f0000ull, 0xffffffffff3f0000ull,
  0xffffffffff3f0000ull, 0xffffffffff7f0000ull, 0xffffffffff7f0000ull, 0xffffffffff7f0000ull,
  0xffffffffff7f0000ull, 0xffffffffff7f0000ull, 0xffffffffff7f0000ull, 0xffffffffffff0000ull,
  0xffffffffffff0000ull, 0xffffffffffff0100ull, 0xffffffffffff0100ull, 0xffffffffffff0100ull,
  0xffffffffffff0100ull, 0xffffffffffff0100ull, 0xffffffffffff0100ull, 0xffffffffffff0300ull,
  0xffffffffffff0300ull, 0xffffffffffff0300ull, 0xffffffffffff0300ull, 0xffffffffffff0700ull,
  0xffffffffffff0700ull, 0xffffffffffff0f00ull, 0xffffffffffff0f00ull, 0xffffffffffff0f00ull,
  0xffffffffffff0f00ull, 0xffffffffffff1f00ull, 0xffffffffffff1f00ull, 0xffffffffffff3f00ull,
  0xffffffffffff3f00ull, 0xffffffffffff3f00ull, 0xffffffffffff3f00ull, 0xffffffffffff7f00ull,
  0xffffffffffff7f00ull, 0xffffffffffff7f00ull, 0xffffffffffff7f00ull, 0xffffffffffff7f00ull,
  0xffffffffffff7f00ull, 0xffffffffffffff00ull, 0xffffffffffffff00ull, 0xffffffffffffff01ull,
  0xffffffffffffff01ull, 0xffffffffffffff01ull, 0xffffffffffffff01ull, 0xffffffffffffff01ull,
  0xffffffffffffff01ull, 0xffffffffffffff03ull, 0xffffffffffffff03ull, 0xffffffffffffff03ull,
  0xffffffffffffff03ull, 0xffffffffffffff07ull, 0xffffffffffffff07ull, 0xffffffffffffff0full,
  0xffffffffffffff0full, 0xffffffffffffff0full, 0xffffffffffffff0full, 0xffffffffffffff1full,
  0xffffffffffffff1full, 0xffffffffffffff3full, 0xffffffffffffff3full, 0xffffffffffffff3full,
  0xffffffffffffff3full, 0xffffffffffffff7full, 0xffffffffffffff7full, 0xffffffffffffff7full,
  0xffffffffffffff7full, 0xffffffffffffff7full, 0xffffffffffffff7full, 0xffffffffffffffffull
};

#else // little-endian CPU architectures

/// Unset bits < start
const std::array<uint64_t, 240> unset_smaller =
{
  ~0ull << 0, ~0ull << 0, ~0ull << 1, ~0ull << 1, ~0ull << 1,
  ~0ull << 1, ~0ull << 1, ~0ull << 1, ~0ull << 2, ~0ull << 2,
  ~0ull << 2, ~0ull << 2, ~0ull << 3, ~0ull << 3, ~0ull << 4,
  ~0ull << 4, ~0ull << 4, ~0ull << 4, ~0ull << 5, ~0ull << 5,
  ~0ull << 6, ~0ull << 6, ~0ull << 6, ~0ull << 6, ~0ull << 7,
  ~0ull << 7, ~0ull << 7, ~0ull << 7, ~0ull << 7, ~0ull << 7,
  ~0ull << 8, ~0ull << 8, ~0ull << 9, ~0ull << 9, ~0ull << 9,
  ~0ull << 9, ~0ull << 9, ~0ull << 9, ~0ull << 10, ~0ull << 10,
  ~0ull << 10, ~0ull << 10, ~0ull << 11, ~0ull << 11, ~0ull << 12,
  ~0ull << 12, ~0ull << 12, ~0ull << 12, ~0ull << 13, ~0ull << 13,
  ~0ull << 14, ~0ull << 14, ~0ull << 14, ~0ull << 14, ~0ull << 15,
  ~0ull << 15, ~0ull << 15, ~0ull << 15, ~0ull << 15, ~0ull << 15,
  ~0ull << 16, ~0ull << 16, ~0ull << 17, ~0ull << 17, ~0ull << 17,
  ~0ull << 17, ~0ull << 17, ~0ull << 17, ~0ull << 18, ~0ull << 18,
  ~0ull << 18, ~0ull << 18, ~0ull << 19, ~0ull << 19, ~0ull << 20,
  ~0ull << 20, ~0ull << 20, ~0ull << 20, ~0ull << 21, ~0ull << 21,
  ~0ull << 22, ~0ull << 22, ~0ull << 22, ~0ull << 22, ~0ull << 23,
  ~0ull << 23, ~0ull << 23, ~0ull << 23, ~0ull << 23, ~0ull << 23,
  ~0ull << 24, ~0ull << 24, ~0ull << 25, ~0ull << 25, ~0ull << 25,
  ~0ull << 25, ~0ull << 25, ~0ull << 25, ~0ull << 26, ~0ull << 26,
  ~0ull << 26, ~0ull << 26, ~0ull << 27, ~0ull << 27, ~0ull << 28,
  ~0ull << 28, ~0ull << 28, ~0ull << 28, ~0ull << 29, ~0ull << 29,
  ~0ull << 30, ~0ull << 30, ~0ull << 30, ~0ull << 30, ~0ull << 31,
  ~0ull << 31, ~0ull << 31, ~0ull << 31, ~0ull << 31, ~0ull << 31,
  ~0ull << 32, ~0ull << 32, ~0ull << 33, ~0ull << 33, ~0ull << 33,
  ~0ull << 33, ~0ull << 33, ~0ull << 33, ~0ull << 34, ~0ull << 34,
  ~0ull << 34, ~0ull << 34, ~0ull << 35, ~0ull << 35, ~0ull << 36,
  ~0ull << 36, ~0ull << 36, ~0ull << 36, ~0ull << 37, ~0ull << 37,
  ~0ull << 38, ~0ull << 38, ~0ull << 38, ~0ull << 38, ~0ull << 39,
  ~0ull << 39, ~0ull << 39, ~0ull << 39, ~0ull << 39, ~0ull << 39,
  ~0ull << 40, ~0ull << 40, ~0ull << 41, ~0ull << 41, ~0ull << 41,
  ~0ull << 41, ~0ull << 41, ~0ull << 41, ~0ull << 42, ~0ull << 42,
  ~0ull << 42, ~0ull << 42, ~0ull << 43, ~0ull << 43, ~0ull << 44,
  ~0ull << 44, ~0ull << 44, ~0ull << 44, ~0ull << 45, ~0ull << 45,
  ~0ull << 46, ~0ull << 46, ~0ull << 46, ~0ull << 46, ~0ull << 47,
  ~0ull << 47, ~0ull << 47, ~0ull << 47, ~0ull << 47, ~0ull << 47,
  ~0ull << 48, ~0ull << 48, ~0ull << 49, ~0ull << 49, ~0ull << 49,
  ~0ull << 49, ~0ull << 49, ~0ull << 49, ~0ull << 50, ~0ull << 50,
  ~0ull << 50, ~0ull << 50, ~0ull << 51, ~0ull << 51, ~0ull << 52,
  ~0ull << 52, ~0ull << 52, ~0ull << 52, ~0ull << 53, ~0ull << 53,
  ~0ull << 54, ~0ull << 54, ~0ull << 54, ~0ull << 54, ~0ull << 55,
  ~0ull << 55, ~0ull << 55, ~0ull << 55, ~0ull << 55, ~0ull << 55,
  ~0ull << 56, ~0ull << 56, ~0ull << 57, ~0ull << 57, ~0ull << 57,
  ~0ull << 57, ~0ull << 57, ~0ull << 57, ~0ull << 58, ~0ull << 58,
  ~0ull << 58, ~0ull << 58, ~0ull << 59, ~0ull << 59, ~0ull << 60,
  ~0ull << 60, ~0ull << 60, ~0ull << 60, ~0ull << 61, ~0ull << 61,
  ~0ull << 62, ~0ull << 62, ~0ull << 62, ~0ull << 62, ~0ull << 63,
  ~0ull << 63, ~0ull << 63, ~0ull << 63, ~0ull << 63, ~0ull << 63
};

/// Unset bits > stop
const std::array<uint64_t, 240> unset_larger =
{
         0ull, ~0ull >> 63, ~0ull >> 63, ~0ull >> 63, ~0ull >> 63,
  ~0ull >> 63, ~0ull >> 63, ~0ull >> 62, ~0ull >> 62, ~0ull >> 62,
  ~0ull >> 62, ~0ull >> 61, ~0ull >> 61, ~0ull >> 60, ~0ull >> 60,
  ~0ull >> 60, ~0ull >> 60, ~0ull >> 59, ~0ull >> 59, ~0ull >> 58,
  ~0ull >> 58, ~0ull >> 58, ~0ull >> 58, ~0ull >> 57, ~0ull >> 57,
  ~0ull >> 57, ~0ull >> 57, ~0ull >> 57, ~0ull >> 57, ~0ull >> 56,
  ~0ull >> 56, ~0ull >> 55, ~0ull >> 55, ~0ull >> 55, ~0ull >> 55,
  ~0ull >> 55, ~0ull >> 55, ~0ull >> 54, ~0ull >> 54, ~0ull >> 54,
  ~0ull >> 54, ~0ull >> 53, ~0ull >> 53, ~0ull >> 52, ~0ull >> 52,
  ~0ull >> 52, ~0ull >> 52, ~0ull >> 51, ~0ull >> 51, ~0ull >> 50,
  ~0ull >> 50, ~0ull >> 50, ~0ull >> 50, ~0ull >> 49, ~0ull >> 49,
  ~0ull >> 49, ~0ull >> 49, ~0ull >> 49, ~0ull >> 49, ~0ull >> 48,
  ~0ull >> 48, ~0ull >> 47, ~0ull >> 47, ~0ull >> 47, ~0ull >> 47,
  ~0ull >> 47, ~0ull >> 47, ~0ull >> 46, ~0ull >> 46, ~0ull >> 46,
  ~0ull >> 46, ~0ull >> 45, ~0ull >> 45, ~0ull >> 44, ~0ull >> 44,
  ~0ull >> 44, ~0ull >> 44, ~0ull >> 43, ~0ull >> 43, ~0ull >> 42,
  ~0ull >> 42, ~0ull >> 42, ~0ull >> 42, ~0ull >> 41, ~0ull >> 41,
  ~0ull >> 41, ~0ull >> 41, ~0ull >> 41, ~0ull >> 41, ~0ull >> 40,
  ~0ull >> 40, ~0ull >> 39, ~0ull >> 39, ~0ull >> 39, ~0ull >> 39,
  ~0ull >> 39, ~0ull >> 39, ~0ull >> 38, ~0ull >> 38, ~0ull >> 38,
  ~0ull >> 38, ~0ull >> 37, ~0ull >> 37, ~0ull >> 36, ~0ull >> 36,
  ~0ull >> 36, ~0ull >> 36, ~0ull >> 35, ~0ull >> 35, ~0ull >> 34,
  ~0ull >> 34, ~0ull >> 34, ~0ull >> 34, ~0ull >> 33, ~0ull >> 33,
  ~0ull >> 33, ~0ull >> 33, ~0ull >> 33, ~0ull >> 33, ~0ull >> 32,
  ~0ull >> 32, ~0ull >> 31, ~0ull >> 31, ~0ull >> 31, ~0ull >> 31,
  ~0ull >> 31, ~0ull >> 31, ~0ull >> 30, ~0ull >> 30, ~0ull >> 30,
  ~0ull >> 30, ~0ull >> 29, ~0ull >> 29, ~0ull >> 28, ~0ull >> 28,
  ~0ull >> 28, ~0ull >> 28, ~0ull >> 27, ~0ull >> 27, ~0ull >> 26,
  ~0ull >> 26, ~0ull >> 26, ~0ull >> 26, ~0ull >> 25, ~0ull >> 25,
  ~0ull >> 25, ~0ull >> 25, ~0ull >> 25, ~0ull >> 25, ~0ull >> 24,
  ~0ull >> 24, ~0ull >> 23, ~0ull >> 23, ~0ull >> 23, ~0ull >> 23,
  ~0ull >> 23, ~0ull >> 23, ~0ull >> 22, ~0ull >> 22, ~0ull >> 22,
  ~0ull >> 22, ~0ull >> 21, ~0ull >> 21, ~0ull >> 20, ~0ull >> 20,
  ~0ull >> 20, ~0ull >> 20, ~0ull >> 19, ~0ull >> 19, ~0ull >> 18,
  ~0ull >> 18, ~0ull >> 18, ~0ull >> 18, ~0ull >> 17, ~0ull >> 17,
  ~0ull >> 17, ~0ull >> 17, ~0ull >> 17, ~0ull >> 17, ~0ull >> 16,
  ~0ull >> 16, ~0ull >> 15, ~0ull >> 15, ~0ull >> 15, ~0ull >> 15,
  ~0ull >> 15, ~0ull >> 15, ~0ull >> 14, ~0ull >> 14, ~0ull >> 14,
  ~0ull >> 14, ~0ull >> 13, ~0ull >> 13, ~0ull >> 12, ~0ull >> 12,
  ~0ull >> 12, ~0ull >> 12, ~0ull >> 11, ~0ull >> 11, ~0ull >> 10,
  ~0ull >> 10, ~0ull >> 10, ~0ull >> 10, ~0ull >> 9, ~0ull >> 9,
  ~0ull >> 9, ~0ull >> 9, ~0ull >> 9, ~0ull >> 9, ~0ull >> 8,
  ~0ull >> 8, ~0ull >> 7, ~0ull >> 7, ~0ull >> 7, ~0ull >> 7,
  ~0ull >> 7, ~0ull >> 7, ~0ull >> 6, ~0ull >> 6, ~0ull >> 6,
  ~0ull >> 6, ~0ull >> 5, ~0ull >> 5, ~0ull >> 4, ~0ull >> 4,
  ~0ull >> 4, ~0ull >> 4, ~0ull >> 3, ~0ull >> 3, ~0ull >> 2,
  ~0ull >> 2, ~0ull >> 2, ~0ull >> 2, ~0ull >> 1, ~0ull >> 1,
  ~0ull >> 1, ~0ull >> 1, ~0ull >> 1, ~0ull >> 1, ~0ull >> 0
};

#endif

} // namespace

#endif

#ifndef ISQRT_UTILS_H
#define ISQRT_UTILS_H

#include <stdint.h>
#include <cmath>

inline uint32_t isqrt(int64_t x)
{
  int64_t root = static_cast<int64_t>(std::sqrt(static_cast<double>(x)));
  // correct rounding errors
  while (root * root > x)
    root--;
  while ((root + 1) * (root + 1) <= x)
    root++;
  return static_cast<uint32_t>(root);
}

inline uint32_t isqrt3(int64_t x)
{
  int64_t root = static_cast<int64_t>(std::pow(static_cast<double>(x), 1.0 / 3.0));
  // correct rounding errors
  while (root * root * root > x)
    root--;
  while ((root + 1) * (root + 1) * (root + 1) <= x)
    root++;
  return static_cast<uint32_t>(root);
}

#endif

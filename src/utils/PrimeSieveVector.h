#ifndef PRIMESIEVEVECTOR_H
#define PRIMESIEVEVECTOR_H

#include <primesieve/soe/PrimeSieveCallback.h>
#include <stdint.h>
#include <vector>

template <typename T>
struct PrimeSieveVector : std::vector<T>, PrimeSieveCallback<uint64_t> {
  void callback(uint64_t prime)
  {
    this->push_back( static_cast<T>(prime) );
  }
};

#endif

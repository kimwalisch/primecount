#ifndef NEXT_N_PRIMES_VECTOR_H
#define NEXT_N_PRIMES_VECTOR_H

#include <primesieve/soe/PrimeSieve.h>
#include <primesieve/soe/PrimeSieveCallback.h>
#include <stdint.h>
#include <vector>
#include <exception>

template <typename T>
class Next_N_Primes_Vector : public std::vector<T>, 
                             public PrimeSieveCallback<uint64_t> {
public:
  void generatePrimes(uint64_t start, uint64_t n)
  {
    i_ = 0;
    n_ = n;
    uint64_t stop = start + n_ * 40;
    PrimeSieve ps;

    try {
      while (i_ < n_)
      {
        ps.generatePrimes(start, stop, this);
        start = stop + 1;
        stop = start + (i_ - n_) * 40;
      }
    } catch (stop_primesieve&) { }
  }

  void callback(uint64_t prime)
  {
    this->push_back( static_cast<T>(prime) );
    if (++i_ == n_)
      throw stop_primesieve();
  }
private:
  uint64_t i_;
  uint64_t n_;
  class stop_primesieve : public std::exception { };
};

#endif

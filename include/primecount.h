#ifndef PRIMECOUNT_H
#define PRIMECOUNT_H

#include <stdint.h>

namespace primecount {

enum {
  /// Uses all CPU cores, Do not modify.
  MAX_THREADS = -1,
  MIN_THREADS = 1
};

int64_t pi_primesieve(int64_t x, int threads = MAX_THREADS);
int64_t pi_legendre(int64_t x, int threads = MAX_THREADS);
int64_t pi_meissel(int64_t x, int threads = MAX_THREADS);
int64_t pi_lehmer(int64_t x, int threads = MAX_THREADS);
int64_t pi(int64_t x, int threads = MAX_THREADS);

int64_t phi(int64_t x, int64_t a, int threads = MAX_THREADS);
int64_t P2(int64_t x, int64_t a, int64_t b, int64_t pb, int threads = MAX_THREADS);
int64_t P3(int64_t x, int64_t a, int64_t c, int64_t pb, int threads = MAX_THREADS);

bool test();

} // namespace primecount

#endif

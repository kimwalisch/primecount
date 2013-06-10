#ifndef PRIMECOUNT_H
#define PRIMECOUNT_H

#include <stdint.h>

enum {
  /// Uses all CPU cores, Do not modify.
  MAX_THREADS = -1,
  MIN_THREADS = 1
};

namespace legendre {
int64_t pi (int64_t x, int threads = MAX_THREADS);
int64_t phi(int64_t x, int64_t a, int threads = MAX_THREADS);
}

namespace meissel {
int64_t pi(int64_t x, int threads = MAX_THREADS);
}

namespace Pkxa {
int64_t P2(int64_t x, int64_t a, int threads = MAX_THREADS);
}

#endif

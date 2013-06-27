#ifndef PHI_PRIMECOUNT_H
#define PHI_PRIMECOUNT_H

#include <primecount.h>
#include <stdint.h>

namespace primecount {

int64_t phi (int64_t x, int64_t a, int threads = MAX_THREADS);
int64_t phi2(int64_t x, int64_t a, int64_t b, int64_t pb, int threads = MAX_THREADS);
int64_t phi3(int64_t x, int64_t a, int64_t c, int64_t pb, int threads = MAX_THREADS);

} // namespace primecount

#endif

#ifndef PHI_LEGENDRE_H
#define PHI_LEGENDRE_H

#include "MAX_THREADS.h"
#include <stdint.h>

namespace legendre {

int64_t phi(int64_t x, int64_t a, int threads = MAX_THREADS);

} // namespace legendre

#endif

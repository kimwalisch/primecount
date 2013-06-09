#ifndef PI_LEGENDRE_H
#define PI_LEGENDRE_H

#include "MAX_THREADS.h"
#include <stdint.h>

namespace legendre {

int64_t pi(int64_t, int threads = MAX_THREADS);

} // namespace legendre

#endif

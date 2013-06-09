#ifndef PI_MEISSEL_H
#define PI_MEISSEL_H

#include <stdint.h>

namespace meissel {

/// By default use all CPU cores
enum { MAX_THREADS = -1 };

int64_t pi(int64_t, int threads = MAX_THREADS);

} // namespace meissel

#endif

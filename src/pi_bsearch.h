#ifndef PI_BSEARCH_H
#define PI_BSEARCH_H

#include <algorithm>
#include <iterator>

namespace primecount {

template <class ForwardIterator, typename T>
inline T pi_bsearch(ForwardIterator first, ForwardIterator last, const T& value)
{
  return static_cast<T>(
      std::distance(first, std::upper_bound(first, last, value)));
}

} // namespace primecount

#endif

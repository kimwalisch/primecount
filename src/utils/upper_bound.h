#ifndef PRIMECOUNT_UPPER_BOUND_H
#define PRIMECOUNT_UPPER_BOUND_H

#include <iterator>

namespace primecount {

/// This is a drop in replacement for std::upper_bound().
/// This implementation uses interpolation search O(log(log(N))) rather
/// than STL's binary search O(log(N)).
///
template <class ForwardIterator, typename T>
inline ForwardIterator
upper_bound(ForwardIterator first, ForwardIterator last, const T& value)
{
    typedef typename std::iterator_traits<ForwardIterator>::difference_type difference_type;
    ForwardIterator mid;
    std::advance(last, -1);

    while (*first <= value && *last >= value)
    {
        mid = first;
        std::advance(mid, static_cast<difference_type>(
                          static_cast<double>(std::distance(first, last)) * (value - *first) / (*last - *first)));

        if (*mid <= value)
        {
          first = mid;
          std::advance(first, 1);
        }
        else if (*mid > value)
        {
          last = mid;
          std::advance(last, -1);
        }
    }

    if (*mid <= value)
      std::advance(mid, 1);

    return mid;
}

} // namespace primecount

#endif

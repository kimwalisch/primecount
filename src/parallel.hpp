///
/// @file  parallel.hpp
/// @brief Fork-join parallelism for primecount. parallel_region()
///        calls body(id) for every id in [0, threads), the calling
///        thread executes id 0 and the other ids are executed by
///        worker threads of primecount's thread pool (see
///        ThreadPool.cpp). parallel_sum() additionally returns the
///        sum of all body(id) results, it replaces OpenMP's
///        reduction(+: sum) clause.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include <Vector.hpp>

#include <type_traits>

namespace primecount {

/// Non-owning reference to a callable with signature void(int).
/// Unlike std::function it never allocates memory.
///
class FunctionRef
{
public:
  template <typename F,
            typename std::enable_if<!std::is_same<typename std::decay<F>::type, FunctionRef>::value, int>::type = 0>
  FunctionRef(F& f) noexcept
    : obj_((void*) &f),
      call_(&call<F>)
  { }

  void operator()(int id) const
  {
    call_(obj_, id);
  }

private:
  template <typename F>
  static void call(void* obj, int id)
  {
    (*(F*) obj)(id);
  }

  void* obj_;
  void (*call_)(void*, int);
};

void parallel_region_impl(int threads, FunctionRef body);

/// Calls body(id) exactly once for every id in [0, threads) and
/// returns after all calls have finished. The calling thread
/// executes id 0, threads <= 1 calls body(0) inline. If calls
/// throw, the first exception is rethrown after all calls have
/// finished.
///
template <typename Body>
void parallel_region(int threads, Body&& body)
{
  if (threads <= 1)
    body(0);
  else
    parallel_region_impl(threads, FunctionRef(body));
}

/// Same as parallel_region(), but returns the sum of all
/// body(id) results. Each id accumulates its partial sum in a
/// local variable of type T and returns it. After all ids have
/// finished, the calling thread adds up the partial sums.
///
template <typename T, typename Body>
T parallel_sum(int threads, Body&& body)
{
  static_assert(std::is_same<decltype(body(0)), T>::value,
                "parallel_sum<T>(): body(id) must return T");

  if (threads <= 1)
    return body(0);

  // Partial sums, stored on the stack for up to 256 ids
  Array<T, 256> array;
  Vector<T> vector;
  T* sums = array.data();

  if (threads > (int) array.size())
  {
    vector.resize(threads);
    sums = vector.data();
  }

  parallel_region(threads, [&](int id) { sums[id] = body(id); });

  T sum = 0;
  for (int i = 0; i < threads; i++)
    sum += sums[i];

  return sum;
}

/// Returns true if the calling thread executes the body
/// of a parallel region with threads > 1.
///
bool in_parallel_region();

/// Maximum number of threads used by default
int max_threads();

// Hooks for tests and benchmarks, not used by algorithms

/// Number of worker threads in the thread pool
int parallel_pool_size();

/// Idle threads busy-wait for this time before they sleep
void set_parallel_spin_time(double seconds);

/// Creating more than max_workers worker threads fails as if
/// std::thread threw an exception. If thread creation has
/// failed before, the thread pool may grow again.
///
void set_parallel_pool_limit(int max_workers);

} // namespace

#endif

///
/// @file  parallel.hpp
/// @brief Fork-join parallelism. parallel_for() calls func(thread_id)
///        for every thread_id in [0, threads[, the calling thread
///        executes thread_id 0 and the other thread ids are executed
///        by worker threads of the thread pool (see ThreadPool.cpp).
///        parallel_sum() additionally returns the sum of all
///        func(thread_id) results, it replaces OpenMP's
///        reduction(+: sum) clause. A function that does not use its
///        thread_id may omit the parameter: func().
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include <macros.hpp>
#include <Vector.hpp>

#include <algorithm>
#include <type_traits>

namespace primecount {

/// Maximum number of threads used by default
int max_threads();

/// Number of worker threads in the thread pool
int thread_pool_size();

/// Idle threads busy-wait for at most ms milliseconds, then they sleep
void set_max_thread_spin_time(double ms);

/// Calls func(thread_id), or func() if func has no
/// thread_id parameter. For a given func only one of the two
/// return types compiles, the other overload is discarded
/// (SFINAE).
///
template <typename Func>
auto invoke_func(const Func& func, int thread_id) -> decltype(func(thread_id))
{
  return func(thread_id);
}

template <typename Func>
auto invoke_func(const Func& func, int) -> decltype(func())
{
  return func();
}

/// Non-owning reference to a callable with signature void(int)
/// or void(), which lives on the caller's stack until the job has
/// finished. Unlike std::function it never allocates memory. Each
/// call runs a copy of the callable on the calling thread's stack,
/// no other code can access the copy, hence the compiler can keep
/// its by-value captures in registers.
///
class FunctionRef
{
public:
  template <typename F>
  FunctionRef(const F& f) noexcept
    : obj_(&f),
      call_(&call<F>)
  {
    static_assert(std::is_copy_constructible<F>::value,
                  "parallel_for(): func must be copyable");
  }

  void operator()(int thread_id) const
  {
    call_(obj_, thread_id);
  }

private:
  template <typename F>
  static void call(const void* obj, int thread_id)
  {
    const F func = *(const F*) obj;
    invoke_func(func, thread_id);
  }

  const void* obj_;
  void (*call_)(const void*, int);
};

void parallel_for_impl(int threads, FunctionRef func);

/// Calls func(thread_id) exactly once for every thread_id in
/// [0, threads[ and returns after all calls have finished.
/// The calling thread executes thread_id 0, threads <= 1 calls
/// func(0) inline. If calls throw, the first exception is
/// rethrown after all calls have finished.
/// @pre threads <= max_threads()
///
template <typename Func>
void parallel_for(int threads, Func&& func)
{
  // Unlike other parallel_for() implementations, every iteration
  // (thread_id) uses its own thread, therefore we limit the number
  // of iterations to the number of CPU cores. This is designed
  // for a thread per CPU core architecture: each thread processes
  // many work items, usually assigned by a load balancer, which
  // provides the best performance for CPU-bound computations.
  ASSERT(threads <= max_threads());

  // threads <= 1 also runs a copy of func, hence func has a
  // single call site, which allows the compiler to inline it.
  FunctionRef funcRef(func);

  if (threads <= 1)
    funcRef(0);
  else
    parallel_for_impl(threads, funcRef);
}

/// Same as parallel_for(), but returns the sum of all
/// func(thread_id) results. Each thread_id accumulates its
/// partial sum in a local variable of type T and returns it,
/// hence hot loops need no atomics or locks. Once all thread
/// ids have finished, the calling thread adds up the sums.
/// @pre threads <= max_threads()
///
template <typename T, typename Func>
T parallel_sum(int threads, Func&& func)
{
  static_assert(std::is_same<decltype(invoke_func(func, 0)), T>::value,
                "parallel_sum<T>(): func(thread_id) must return T");

  // Partial sums, stored on the stack for up to 64 thread ids
  // to avoid a memory allocation per job.
  INDETERMINATE Array<T, 64> array;
  Vector<T> vector;
  T* sums = array.data();
  threads = std::max(threads, 1);

  if (threads > (int) array.size())
  {
    vector.resize(threads);
    sums = vector.data();
  }

  // func is captured by value, hence the caller's func
  // never escapes and func has a single call site.
  parallel_for(threads, [func, sums](int thread_id)
  {
    sums[thread_id] = invoke_func(func, thread_id);
  });

  T sum = 0;
  for (int i = 0; i < threads; i++)
    sum += sums[i];

  return sum;
}

} // namespace

#endif

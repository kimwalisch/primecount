///
/// @file  ThreadPool.cpp
/// @brief primecount's thread pool, it implements the parallel
///        regions of parallel.hpp.
///
///        The thread pool creates worker threads on demand, up to
///        max_threads() - 1, and never frees them. The leader of a
///        parallel region (the calling thread) leases up to
///        threads - 1 idle workers from the pool. The leader and its
///        leased workers claim the region's ids from an atomic
///        counter until no ids are left, hence a region is correct
///        no matter how many workers it gets. Once all leased
///        workers have checked out, the leader returns them to the
///        pool. A region started inside the body of another region
///        leases its own workers (nested parallelism), regions at
///        the third nesting level run on the calling thread.
///
///        Idle workers busy-wait for their next region for up to
///        spin_time before they sleep, this way consecutive
///        parallel regions start without OS thread wake-ups.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <parallel.hpp>
#include <primecount-config.hpp>
#include <macros.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <system_error>
#include <thread>

#if !defined(_WIN32) && \
    !defined(__EMSCRIPTEN__) && \
    __has_include(<pthread.h>)
  #include <pthread.h>
  #define HAVE_PTHREAD_ATFORK
#endif

#if defined(__linux__)
  #include <sched.h>
#endif

#if defined(_MSC_VER) && \
   (defined(_M_X64) || defined(_M_IX86) || defined(_M_ARM64))
  #include <intrin.h>
#endif

namespace {

using namespace primecount;
using Clock = std::chrono::steady_clock;

/// Idle threads check the clock every SPIN_CHECK iterations
constexpr unsigned SPIN_CHECK = 1024;

struct Region;

/// The leader of a parallel region sleeps on its LeaderSync if its
/// workers take long to check out. LeaderSync objects are never
/// freed: a worker may still notify a LeaderSync after its leader
/// has returned, which is a harmless spurious wake-up.
///
struct LeaderSync
{
  std::atomic<bool> sleeping{false};
  std::mutex mutex;
  std::condition_variable cond;
  LeaderSync* next = nullptr;
};

struct Worker
{
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  // Region assigned by a leader, nullptr while idle
  std::atomic<Region*> mailbox{nullptr};
  std::atomic<bool> sleeping{false};
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  std::mutex mutex;
  std::condition_variable cond;
  // Next idle worker in the pool or next worker of the
  // leader's team, only used by the pool and the leader.
  Worker* next = nullptr;
};

/// Region state, it lives on the leader's stack. This is safe
/// because the leader waits until all its workers have checked
/// out and workers never access the region afterwards.
///
struct Region
{
  Region(FunctionRef body, int threads, int depth)
    : body(body),
      threads(threads),
      depth(depth)
  { }

  FunctionRef body;
  LeaderSync* sync = nullptr;
  int threads;
  int depth;
  std::exception_ptr exception;
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  // The leader executes id 0
  std::atomic<int> next_id{1};
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  std::atomic<int> checked_out{0};
  std::atomic<bool> has_exception{false};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
};

/// Workers leased by the leader of a parallel region
struct Team
{
  // Linked list of workers using Worker::next
  Worker* first = nullptr;
  Worker* last = nullptr;
  LeaderSync* sync = nullptr;
  int size = 0;
  // New workers to create, their indexes in
  // the pool are [slot, slot + reserved[.
  int reserved = 0;
  int slot = 0;
  // Test hook: creating worker index >= limit fails
  int limit = 0;
};

class ThreadPool
{
public:
  ThreadPool();
  Team lease(int workers);
  void create_workers(Team& team, Region& region);
  void give_back(Team& team);
  int size();
  void set_limit(int max_workers);
  void set_spin_time(double seconds);
  Clock::time_point spin_deadline() const;

private:
  std::atomic<int> spin_time_us_{10000};
  std::mutex mutex_;
  Worker* idle_ = nullptr;
  LeaderSync* syncs_ = nullptr;
  int size_ = 0;
  int max_size_ = 0;
  int limit_ = INT_MAX;
  bool can_grow_ = true;
};

thread_local int region_depth = 0;
std::atomic<ThreadPool*> thread_pool{nullptr};

ALWAYS_INLINE void cpu_pause()
{
#if defined(_MSC_VER) && \
   (defined(_M_X64) || defined(_M_IX86))
  _mm_pause();
#elif defined(_MSC_VER) && \
      defined(_M_ARM64)
  __yield();
#elif (defined(__GNUC__) || defined(__clang__)) && \
      (defined(__x86_64__) || defined(__i386__))
  __builtin_ia32_pause();
#elif (defined(__GNUC__) || defined(__clang__)) && \
       defined(__aarch64__)
  __asm__ __volatile__("isb" ::: "memory");
#endif
}

#if defined(HAVE_PTHREAD_ATFORK)

/// After fork() the child process only contains the thread
/// that called fork(). We leak the parent's thread pool, the
/// child creates a new thread pool when needed.
///
void fork_child()
{
  thread_pool.store(nullptr, std::memory_order_relaxed);
}

#endif

ThreadPool& get_pool()
{
  ThreadPool* pool = thread_pool.load(std::memory_order_acquire);

  if (!pool)
  {
    ThreadPool* new_pool = new ThreadPool;

    if (thread_pool.compare_exchange_strong(pool, new_pool, std::memory_order_acq_rel))
    {
      pool = new_pool;
#if defined(HAVE_PTHREAD_ATFORK)
      static const int registered = pthread_atfork(nullptr, nullptr, fork_child);
      (void) registered;
#endif
    }
    else
      delete new_pool;
  }

  return *pool;
}

void run_id(Region& region, int id)
{
  try
  {
    region.body(id);
  }
  catch (...)
  {
    if (!region.has_exception.exchange(true, std::memory_order_relaxed))
      region.exception = std::current_exception();
  }
}

void claim_ids(Region& region)
{
  int threads = region.threads;

  for (int id = region.next_id.fetch_add(1, std::memory_order_relaxed);
       id < threads;
       id = region.next_id.fetch_add(1, std::memory_order_relaxed))
    run_id(region, id);
}

/// Wakes up a worker if it sleeps. The sequentially consistent
/// store of the mailbox and load of the sleeping flag pair with
/// the worker's store of its sleeping flag and load of its
/// mailbox: either the worker finds its region or we find the
/// worker sleeping, hence no wake-up is lost.
///
void post(Worker* worker, Region* region)
{
  worker->mailbox.store(region, std::memory_order_seq_cst);

  if (worker->sleeping.load(std::memory_order_seq_cst))
  {
    std::lock_guard<std::mutex> lock(worker->mutex);
    worker->cond.notify_one();
  }
}

Region* wait_for_region(ThreadPool* pool, Worker* worker)
{
  Clock::time_point deadline = pool->spin_deadline();
  Region* region = worker->mailbox.load(std::memory_order_acquire);

  for (unsigned i = 1; !region; i++)
  {
    cpu_pause();
    region = worker->mailbox.load(std::memory_order_acquire);

    if (!region &&
        i % SPIN_CHECK == 0 &&
        Clock::now() >= deadline)
    {
      std::unique_lock<std::mutex> lock(worker->mutex);
      worker->sleeping.store(true, std::memory_order_seq_cst);

      while (!(region = worker->mailbox.load(std::memory_order_seq_cst)))
        worker->cond.wait(lock);

      worker->sleeping.store(false, std::memory_order_relaxed);
    }
  }

  worker->mailbox.store(nullptr, std::memory_order_relaxed);
  return region;
}

/// Executes ids of the region until no ids are left, then checks
/// out. The region may be destroyed as soon as the last worker
/// has checked out, hence we read the leader's sync object
/// before checking out.
///
void run_worker(Region* region)
{
  region_depth = region->depth;
  claim_ids(*region);
  region_depth = 0;

  LeaderSync* sync = region->sync;
  region->checked_out.fetch_add(1, std::memory_order_seq_cst);

  if (sync->sleeping.load(std::memory_order_seq_cst))
  {
    std::lock_guard<std::mutex> lock(sync->mutex);
    sync->cond.notify_one();
  }
}

void worker_main(ThreadPool* pool, Worker* worker, Region* region)
{
  while (true)
  {
    run_worker(region);
    region = wait_for_region(pool, worker);
  }
}

/// The leader waits until all its workers have checked out.
/// The acquire load of checked_out makes all writes of the
/// workers visible to the leader.
///
void wait_for_team(ThreadPool* pool, Region& region, int workers)
{
  Clock::time_point deadline = pool->spin_deadline();

  for (unsigned i = 1; region.checked_out.load(std::memory_order_acquire) < workers; i++)
  {
    cpu_pause();

    if (i % SPIN_CHECK == 0 &&
        Clock::now() >= deadline)
    {
      LeaderSync* sync = region.sync;
      std::unique_lock<std::mutex> lock(sync->mutex);
      sync->sleeping.store(true, std::memory_order_seq_cst);

      while (region.checked_out.load(std::memory_order_seq_cst) < workers)
        sync->cond.wait(lock);

      sync->sleeping.store(false, std::memory_order_relaxed);
      break;
    }
  }
}

/// Runs all ids on the calling thread, used for
/// regions at the third nesting level.
///
void run_inline(int threads, FunctionRef body)
{
  std::exception_ptr exception;

  for (int id = 0; id < threads; id++)
  {
    try
    {
      body(id);
    }
    catch (...)
    {
      if (!exception)
        exception = std::current_exception();
    }
  }

  if (exception)
    std::rethrow_exception(exception);
}

ThreadPool::ThreadPool()
  : max_size_(max_threads() - 1)
{ }

/// Leases up to the requested number of idle workers and
/// reserves new workers if the pool may grow.
///
Team ThreadPool::lease(int workers)
{
  Team team;
  std::lock_guard<std::mutex> lock(mutex_);

  if (!idle_ && !(can_grow_ && size_ < max_size_))
    return team;

  if (syncs_)
  {
    team.sync = syncs_;
    syncs_ = syncs_->next;
  }
  else
    team.sync = new LeaderSync;

  while (idle_ && team.size < workers)
  {
    Worker* worker = idle_;
    idle_ = worker->next;
    worker->next = team.first;
    team.first = worker;
    if (!team.last)
      team.last = worker;
    team.size++;
  }

  if (team.size < workers && can_grow_)
  {
    team.reserved = std::min(workers - team.size, max_size_ - size_);
    team.slot = size_;
    team.limit = limit_;
    size_ += team.reserved;
  }

  return team;
}

/// Creates the reserved workers, they start working on the region
/// immediately. If thread creation fails, the pool stops growing.
///
void ThreadPool::create_workers(Team& team, Region& region)
{
  int created = 0;

  try
  {
    for (; created < team.reserved; created++)
    {
      if (team.slot + created >= team.limit)
        throw std::system_error(std::make_error_code(std::errc::resource_unavailable_try_again));

      Worker* worker = new Worker;

      try
      {
        std::thread(worker_main, this, worker, &region).detach();
      }
      catch (...)
      {
        delete worker;
        throw;
      }

      worker->next = team.first;
      team.first = worker;
      if (!team.last)
        team.last = worker;
      team.size++;
    }
  }
  catch (...)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    size_ -= team.reserved - created;
    can_grow_ = false;
  }
}

/// Returns the leased workers to the pool, the most recently
/// used workers are leased first next time.
///
void ThreadPool::give_back(Team& team)
{
  std::lock_guard<std::mutex> lock(mutex_);

  if (team.last)
  {
    team.last->next = idle_;
    idle_ = team.first;
  }

  team.sync->next = syncs_;
  syncs_ = team.sync;
}

int ThreadPool::size()
{
  std::lock_guard<std::mutex> lock(mutex_);
  return size_;
}

void ThreadPool::set_limit(int max_workers)
{
  std::lock_guard<std::mutex> lock(mutex_);
  limit_ = max_workers;
  can_grow_ = true;
}

void ThreadPool::set_spin_time(double seconds)
{
  double us = std::min(seconds * 1e6, (double) INT_MAX);
  spin_time_us_.store((int) std::max(us, 0.0), std::memory_order_relaxed);
}

/// Idle threads busy-wait until the deadline, then they sleep
Clock::time_point ThreadPool::spin_deadline() const
{
  int us = spin_time_us_.load(std::memory_order_relaxed);
  return Clock::now() + std::chrono::microseconds(us);
}

} // namespace

namespace primecount {

void parallel_region_impl(int threads, FunctionRef body)
{
  if (threads <= 1)
  {
    body(0);
    return;
  }

  int depth = region_depth;

  if (depth >= 2)
  {
    run_inline(threads, body);
    return;
  }

  ThreadPool& pool = get_pool();
  INDETERMINATE Region region(body, threads, depth + 1);
  Team team = pool.lease(threads - 1);
  region.sync = team.sync;

  // Idle workers start immediately, new
  // workers start once they are created.
  for (Worker* worker = team.first; worker; worker = worker->next)
    post(worker, &region);

  if (team.reserved > 0)
    pool.create_workers(team, region);

  region_depth = depth + 1;
  run_id(region, 0);
  claim_ids(region);
  region_depth = depth;

  if (team.size > 0)
    wait_for_team(&pool, region, team.size);
  if (team.sync)
    pool.give_back(team);

  if (region.exception)
    std::rethrow_exception(region.exception);
}

bool in_parallel_region()
{
  return region_depth > 0;
}

int max_threads()
{
  static const int threads = []() -> int
  {
    int n = (int) std::thread::hardware_concurrency();

#if defined(__linux__) && \
    defined(CPU_COUNT)
    // Respect the CPU affinity mask, e.g. taskset
    cpu_set_t cpus;
    if (sched_getaffinity(0, sizeof(cpus), &cpus) == 0)
      n = CPU_COUNT(&cpus);
#endif

    return std::max(n, 1);
  }();

  return threads;
}

int parallel_pool_size()
{
  return get_pool().size();
}

void set_parallel_spin_time(double seconds)
{
  get_pool().set_spin_time(seconds);
}

void set_parallel_pool_limit(int max_workers)
{
  get_pool().set_limit(max_workers);
}

} // namespace

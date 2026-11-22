///
/// @file  ThreadPool.cpp
/// @brief The thread pool creates worker threads on demand, up to
///        2 * max_threads() - 1, and never frees them. The main
///        thread of a parallel job (the calling thread) borrows idle
///        workers from the pool, the threads that no borrowed worker
///        runs are claimed from an atomic counter. Jobs may be
///        nested.
///
///        The borrowed workers start each other in a binomial tree
///        and report their completion back up the same tree, hence
///        starting and joining a job take O(log threads) steps. Idle
///        workers busy-wait before they sleep, this way consecutive
///        jobs start without OS thread wake-ups.
///
///        Each waiting thread measures its average wake-up time,
///        starting from 20 us. It busy-waits for 4x the time needed
///        to wake up all levels of its last job's binomial tree, but
///        at most for 3 ms (max_spin_time). While a job creates
///        workers, waiting threads keep spinning, otherwise Linux
///        gives their CPUs to the new threads and after waking up
///        two threads may share a CPU. When more threads are busy
///        than there are CPU cores, spinning is disabled and
///        waiting threads sleep instead.
///
///        Like OpenMP, a job runs with all its threads, even if other
///        jobs already use all CPUs. If the pool is full, the job
///        creates temporary workers that exit after the job.
///
///        fork() is not supported: if the parent process has used the
///        thread pool, parallel jobs in a child process created by
///        fork() may hang (similar to GCC OpenMP).
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "parallel.hpp"

#include <primecount-config.hpp>
#include <imath.hpp>
#include <macros.hpp>
#include <popcnt.hpp>

#include <stdint.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>

#if defined(__linux__) && \
    __has_include(<sched.h>)
  #include <sched.h>
  #if defined(CPU_COUNT)
    #define HAVE_SCHED_GETAFFINITY
  #endif
#endif

#if defined(_WIN32) && \
    __has_include(<windows.h>)
  #if !defined(NOMINMAX)
    #define NOMINMAX
  #endif
  #include <windows.h>
  // get_cpus_win32() requires Windows 7
  #if defined(_WIN32_WINNT) && \
      _WIN32_WINNT >= 0x0601
    #define HAVE_GET_CPUS_WIN32
  #endif
  // WaitOnAddress() requires Windows 8
  #if defined(_WIN32_WINNT) && \
      _WIN32_WINNT >= 0x0602
    #define HAVE_WAIT_ON_ADDRESS
  #endif
#endif

#if defined(_MSC_VER) && \
   !defined(__clang__) && \
   (defined(_M_X64) || defined(_M_IX86) || defined(_M_ARM64))
  #include <intrin.h>
#endif

// See cmake/pin_library.cmake
#if defined(HAVE_PIN_LIBRARY)
  #include <dlfcn.h>
#endif

namespace {

using namespace primecount;
using Clock = std::chrono::steady_clock;

// Native atomic waits support 4 bytes on Linux (futex) and
// 8 bytes elsewhere, e.g. with libc++ on macOS and FreeBSD.
#if defined(__linux__)
  using epoch_t = int32_t;
#else
  using epoch_t = int64_t;
#endif

class Job;
struct Worker;

/// A waiting thread uses its Sleeper to sleep if spinning takes
/// too long. The thread sets its sleeping flag before its last
/// check of ready(), the thread that wakes it up checks the flag
/// after making ready() true, hence no wake-up is lost. The epoch
/// is read before ready() is checked, hence a wake-up between
/// the check and wait() makes wait() return immediately.
///
class Sleeper
{
public:
  template <typename Ready>
  void sleep_until(Ready ready)
  {
    bool slept = false;

    while (true)
    {
      epoch_t epoch = epoch_.load(std::memory_order_seq_cst);
      sleeping_.store(true, std::memory_order_seq_cst);
      if (ready())
        break;
      wait(epoch);
      slept = true;
    }

    sleeping_.store(false, std::memory_order_relaxed);

    if (slept)
    {
      // Read before Clock::now(), hence the sample is >= 0
      int64_t wake_time = wake_time_.load(std::memory_order_relaxed);
      std::chrono::nanoseconds now = Clock::now().time_since_epoch();
      add_sample(now.count() - wake_time);
    }
  }

  /// Wakes up the thread after its condition became true
  void wake_up()
  {
    if (sleeping_.load(std::memory_order_seq_cst))
    {
      std::chrono::nanoseconds now = Clock::now().time_since_epoch();
      wake_time_.store(now.count(), std::memory_order_relaxed);
      epoch_.fetch_add(1, std::memory_order_seq_cst);
      notify();
    }
  }

  /// Average time from wake_up() until the thread runs again.
  /// Only the thread that sleeps on the Sleeper accesses it.
  ///
  int64_t wake_ns() const
  {
    return wake_ns_;
  }

private:
  /// Waits until epoch_ != epoch, may return spuriously:
  /// WaitOnAddress() may return while epoch_ is unchanged.
  /// std::atomic::wait() and the condition variable loop do
  /// not, they return only once epoch_ has changed.
  ///
  void wait(epoch_t epoch)
  {
#if defined(HAVE_WAIT_ON_ADDRESS)
    WaitOnAddress(&epoch_, &epoch, sizeof(epoch), INFINITE);
#elif __cplusplus >= 202002L
    epoch_.wait(epoch, std::memory_order_seq_cst);
#else
    std::unique_lock<std::mutex> lock(mutex_);
    while (epoch_.load(std::memory_order_seq_cst) == epoch)
      cond_.wait(lock);
#endif
  }

  void notify()
  {
#if defined(HAVE_WAIT_ON_ADDRESS)
    WakeByAddressSingle(&epoch_);
#elif __cplusplus >= 202002L
    epoch_.notify_one();
#else
    std::lock_guard<std::mutex> lock(mutex_);
    cond_.notify_one();
#endif
  }

  /// Clipped running average, spurious wake-ups (e.g. a late
  /// wake-up of a reused MainThreadSync) may cause outliers.
  ///
  void add_sample(int64_t ns)
  {
    // An outlier counts at most 4x the average, the 1 us
    // keeps the average from getting stuck at 0.
    int64_t max_ns = 4 * wake_ns_ + 1000;
    ns = in_between(0, ns, max_ns);
    wake_ns_ += (ns - wake_ns_) / 8;
  }

  std::atomic<bool> sleeping_{false};
  std::atomic<epoch_t> epoch_{0};
  std::atomic<int64_t> wake_time_{0};
  // Until the thread has slept once, we assume
  // a typical wake-up latency of 20 us.
  int64_t wake_ns_ = 20000;

#if !defined(HAVE_WAIT_ON_ADDRESS) && \
    __cplusplus < 202002L
  std::mutex mutex_;
  std::condition_variable cond_;
#endif
};

/// The main thread of a parallel job sleeps on its MainThreadSync
/// if its workers take long to finish. The last worker of a job
/// may still wake up the MainThreadSync after the main thread has
/// returned from parallel_for(), hence MainThreadSync objects are
/// never freed. If a later job already uses the MainThreadSync,
/// this is a harmless spurious wake-up.
///
struct MainThreadSync
{
  /// Writes the entries of the job's workers, hence their memory
  /// pages are mapped before borrow() takes the pool lock.
  MainThreadSync(int max_workers, int job_workers)
    : workers(max_workers)
  {
    initialized = in_between(0, job_workers, max_workers);
    std::fill_n(workers.begin(), initialized, nullptr);
  }

  // Every job's last report_finished() reads the sleeper, the
  // padding separates it from next and from neighboring heap
  // memory, which the main thread writes in give_back().
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  Sleeper sleeper;
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  MainThreadSync* next = nullptr;
  // Workers of the main thread's group, see ThreadPool::borrow().
  // Only workers[0, initialized[ have been written, zeroing the
  // whole array would take microseconds on large servers.
  Vector<Worker*> workers;
  int initialized = 0;
  // The main thread writes next in give_back(), the padding
  // separates it from the worker array that the workers read.
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
};

struct Worker
{
  // Each idle worker spins on its own job pointer, hence
  // starting a worker only touches its cache line. The padding
  // prevents false sharing between workers.
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  // Job to work on, nullptr while idle
  std::atomic<Job*> job{nullptr};
  // Position in the main thread's group, set together with job
  int position = 0;
  Sleeper sleeper;
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  // Unfinished parts of the worker's subtree: 1 for its own
  // threads plus 1 per child subtree, see Job::report_finished().
  std::atomic<int> pending{0};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
};

/// Workers borrowed by the main thread of a parallel job
struct WorkerGroup
{
  // The workers are sync->workers[0, size[
  MainThreadSync* sync = nullptr;
  int size = 0;
  // New workers to create
  int new_workers = 0;
  // Temporary workers to create, they exit after the job
  int temp_workers = 0;
};

/// Job state, it lives on the main thread's stack. This is
/// safe because the main thread waits until all its workers have
/// finished and workers never access the job afterwards.
///
class Job
{
public:
  Job(FunctionRef func, int threads, const WorkerGroup& group);
  void start();
  void start_children(Worker& worker);
  void run_thread(int thread_id) noexcept;
  void run_remaining_threads();
  void report_finished(int position);
  bool finished() const;
  int threads() const;
  void rethrow_exception() const;

private:
  Worker* worker_at(int position) const;
  void start_worker(int position);

  FunctionRef func_;
  MainThreadSync* sync_;
  // The borrowed workers are workers_[0, borrowed_[,
  // workers_[i] has the position i + 1.
  Worker** workers_;
  int borrowed_;
  int threads_;
  std::exception_ptr exception_;
  // Many threads update next_thread_id_ and pending_, the
  // padding separates them from each other and from the
  // read-only fields above to avoid false sharing.
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<int> next_thread_id_;
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  // Unfinished parts of the job: 1 per subtree of the
  // main thread plus 1 per new or temporary worker.
  std::atomic<int> pending_;
  std::atomic<bool> has_exception_{false};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
};

/// The thread pool is never destroyed: its detached persistent
/// workers run until the process exits. Hence it never frees
/// its Worker and MainThreadSync objects either.
///
class ThreadPool
{
public:
  ThreadPool();
  void run(int threads, FunctionRef func);
  int size();
  void set_max_spin_time(double ms);

private:
  WorkerGroup borrow(int workers);
  void create_workers(WorkerGroup& group, Job& job) noexcept;
  void wait_for_workers(const WorkerGroup& group, const Job& job);
  void give_back(WorkerGroup& group);
  void worker_loop(Worker* worker, Job* job) const;
  Job* wait_for_job(Worker& worker, int spin_factor) const;
  template <typename Ready>
  void wait_until(Sleeper& sleeper, int spin_factor, Ready ready) const;
  static int get_spin_factor(int threads);
  Clock::time_point spin_deadline(const Sleeper& sleeper, int spin_factor) const;

  // Every waiting thread reads max_spin_time_us_, the padding
  // separates it from the fields that main threads write in
  // every job.
  MAYBE_UNUSED char pad1[MAX_CACHE_LINE_SIZE];
  std::atomic<int> max_spin_time_us_{3000};
  // Waiting threads sleep instead of spinning while
  // busy_workers_ >= cpus_, it changes only when
  // busy_workers_ crosses cpus_.
  std::atomic<bool> oversubscribed_{false};
  // Number of jobs creating workers, see wait_until(). It
  // changes rarely, only when a job creates workers.
  std::atomic<int> creating_workers_{0};
  MAYBE_UNUSED char pad2[MAX_CACHE_LINE_SIZE];
  // Spinning threads read waiting_main_threads_, see wait_until()
  std::atomic<int> waiting_main_threads_{0};
  MAYBE_UNUSED char pad3[MAX_CACHE_LINE_SIZE];
  std::mutex mutex_;
  // Unused MainThreadSync objects, linked by next. borrow()
  // also reads it right before locking mutex_, hence both
  // should be in the same cache line.
  std::atomic<MainThreadSync*> unused_syncs_{nullptr};
  // Stack of idle workers, the most recently used on top
  Vector<Worker*> idle_workers_;
  int size_ = 0;
  int max_size_ = 0;
  // Number of CPUs, see oversubscribed_
  int cpus_ = 0;
  // Borrowed, new and temporary workers of all jobs
  int busy_workers_ = 0;
  // Main threads write the fields above in every job, the
  // padding separates them from neighboring heap memory.
  MAYBE_UNUSED char pad4[MAX_CACHE_LINE_SIZE];
};

/// The pause instruction makes busy-waiting cheaper: it saves
/// power and leaves more resources to the other hyper-thread of
/// the same CPU core. MSVC is checked last because clang-cl
/// also defines _MSC_VER but should use the same code as Clang.
///
ALWAYS_INLINE void cpu_pause()
{
#if (defined(__GNUC__) || defined(__clang__)) && \
    (defined(__x86_64__) || defined(__i386__))
  __builtin_ia32_pause();
#elif (defined(__GNUC__) || defined(__clang__)) && \
       defined(__aarch64__)
  __asm__ __volatile__("isb" ::: "memory");
#elif defined(_MSC_VER) && \
     (defined(_M_X64) || defined(_M_IX86))
  _mm_pause();
#elif defined(_MSC_VER) && \
      defined(_M_ARM64)
  __isb(_ARM64_BARRIER_SY);
#endif
}

/// The main thread runs thread 0 and the borrowed workers run
/// the threads of their positions, the remaining threads are
/// claimed, e.g. by new workers. Each subtree of the main thread
/// reports its completion once, new and temporary workers report
/// individually.
///
Job::Job(FunctionRef func,
         int threads,
         const WorkerGroup& group)
  : func_(func),
    sync_(group.sync),
    workers_(group.sync->workers.data()),
    borrowed_(group.size),
    threads_(threads),
    next_thread_id_(group.size + 1),
    pending_(group.new_workers + group.temp_workers +
        (group.size > 0 ? ilog2(group.size) + 1 : 0))
{ }

/// The main thread starts the positions 1, 2, 4, 8, ...
void Job::start()
{
  for (int bit = 1; bit <= borrowed_; bit *= 2)
    start_worker(bit);
}

/// Binomial tree: the main thread has position 0 and the borrowed
/// workers have the positions 1 to borrowed_. The worker at
/// position p starts the workers at the positions p + 2^i for
/// all 2^i > p, e.g. the worker at position 1 starts the workers
/// at the positions 3, 5, 9, 17, ... It starts the worker with
/// the largest subtree first. Each started worker immediately
/// starts its own workers, hence the number of started workers
/// doubles in each step.
///
void Job::start_children(Worker& worker)
{
  int p = worker.position;
  int first = 2 << ilog2(p);
  int children = 0;

  for (int bit = first; p + bit <= borrowed_; bit *= 2)
    children++;

  // Must be set before the children can finish
  worker.pending.store(1 + children, std::memory_order_relaxed);

  for (int bit = first; p + bit <= borrowed_; bit *= 2)
    start_worker(p + bit);
}

/// The main thread rethrows the first exception after all
/// threads have finished. Only std::exception subclasses are
/// caught, other exceptions call std::terminate() due to
/// noexcept.
///
void Job::run_thread(int thread_id) noexcept
{
  try
  {
    func_(thread_id);
  }
  catch (const std::exception&)
  {
    if (!has_exception_.exchange(true, std::memory_order_relaxed))
      exception_ = std::current_exception();
  }
}

void Job::run_remaining_threads()
{
  int threads = threads_;

  // All threads may be assigned to worker positions,
  // then a load avoids contention on next_thread_id_.
  if (next_thread_id_.load(std::memory_order_relaxed) >= threads)
    return;

  for (int thread_id = next_thread_id_.fetch_add(1, std::memory_order_relaxed);
       thread_id < threads;
       thread_id = next_thread_id_.fetch_add(1, std::memory_order_relaxed))
    run_thread(thread_id);
}

/// In the tree of borrowed workers, the last thread to finish
/// in a subtree reports its completion to the parent, the
/// position without its highest bit. Hence each subtree of
/// the main thread reports only once.
///
void Job::report_finished(int position)
{
  for (int p = position; p > 0; p -= 1 << ilog2(p))
    if (worker_at(p)->pending.fetch_sub(1, std::memory_order_acq_rel) != 1)
      return;

  // The job may be destroyed as soon as pending_ is 0,
  // hence we read the main thread's sync object before.
  MainThreadSync* sync = sync_;

  if (pending_.fetch_sub(1, std::memory_order_seq_cst) == 1)
    sync->sleeper.wake_up();
}

/// The read-modify-writes of the pending counts pass
/// the workers' writes on to the main thread.
///
bool Job::finished() const
{
  return pending_.load(std::memory_order_seq_cst) == 0;
}

int Job::threads() const
{
  return threads_;
}

void Job::rethrow_exception() const
{
  if (exception_)
    std::rethrow_exception(exception_);
}

Worker* Job::worker_at(int position) const
{
  return workers_[position - 1];
}

/// Starts the worker at the given position and wakes it up if
/// it sleeps. The sequentially consistent store waits until
/// the worker's cache line arrives, see Job::start_children().
///
void Job::start_worker(int position)
{
  Worker* worker = worker_at(position);
  worker->position = position;
  worker->job.store(this, std::memory_order_seq_cst);
  worker->sleeper.wake_up();
}

ThreadPool& get_pool()
{
  static ThreadPool* pool = new ThreadPool;
  return *pool;
}

/// The thread pool's workers run the library's code until the
/// process exits. If the shared library were unloaded, e.g. by
/// dlclose() when a language binding unloads its C module at the
/// end of a script, workers that are still busy-waiting after
/// the last parallel job would crash the process. Idle workers
/// would leak, and loading the library again would create another
/// thread pool.
///
/// Hence we prevent the library (or executable) containing the
/// thread pool from being unloaded: dlopen() with RTLD_NOLOAD
/// increments the reference count of the already loaded library
/// and we never close the handle. RTLD_NODELETE additionally
/// marks the library as never to be unloaded. Like Rust's Rayon
/// library, the thread pool is never shut down.
///
void pin_library()
{
#if defined(HAVE_PIN_LIBRARY)
  static char address;
  Dl_info info;

  if (dladdr(&address, &info) && info.dli_fname)
    dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD | RTLD_NODELETE);
#endif
}

ThreadPool::ThreadPool()
  : max_size_(2 * max_threads() - 1),
    cpus_(max_threads())
{
  idle_workers_.reserve(max_size_);
  // Before the first worker is created
  pin_library();
}

/// Busy-waits until ready() returns true, the thread sleeps after
/// spin_deadline(). Linux may queue a woken worker on a CPU that
/// runs one of our spinning threads and then not run it for a
/// whole time slice (~3 ms). Hence spinning threads yield while a
/// main thread waits, yielding is cheap if no other thread waits
/// for their CPU.
///
template <typename Ready>
void ThreadPool::wait_until(Sleeper& sleeper, int spin_factor, Ready ready) const
{
  // Reading the clock is slower than spinning
  constexpr unsigned SPIN_CHECK = 1024;
  // Yield while a main thread waits for its workers
  constexpr unsigned YIELD_CHECK = 128;

  Clock::time_point deadline = spin_deadline(sleeper, spin_factor);

  for (unsigned i = 1; !ready(); i++)
  {
    cpu_pause();

    if (i % YIELD_CHECK == 0)
    {
      // Spinning would take CPU time from the busy threads
      if (oversubscribed_.load(std::memory_order_relaxed))
      {
        sleeper.sleep_until(ready);
        return;
      }

      if (waiting_main_threads_.load(std::memory_order_relaxed) > 0)
        std::this_thread::yield();
    }

    if (i % SPIN_CHECK == 0 &&
        Clock::now() >= deadline)
    {
      // Don't sleep while a job creates workers (50 - 100 us
      // per thread), otherwise a new thread may take our CPU.
      // max_spin_time = 0 disables all spinning.
      if (creating_workers_.load(std::memory_order_relaxed) > 0 &&
          max_spin_time_us_.load(std::memory_order_relaxed) > 0)
      {
        deadline = spin_deadline(sleeper, spin_factor);
        continue;
      }

      sleeper.sleep_until(ready);
      return;
    }
  }
}

/// Each worker thread runs this loop until the process exits,
/// temporary workers (worker = nullptr) exit after their job.
///
void ThreadPool::worker_loop(Worker* worker, Job* job) const
{
  // Used when waiting for the next job, which
  // likely has the same number of threads.
  int spin_factor = get_spin_factor(job->threads());

  // New workers have no position
  job->run_remaining_threads();
  job->report_finished(0);

  if (!worker)
    return;

  while (true)
  {
    job = wait_for_job(*worker, spin_factor);
    int position = worker->position;
    spin_factor = get_spin_factor(job->threads());
    job->start_children(*worker);
    job->run_thread(position);
    job->run_remaining_threads();
    job->report_finished(position);
  }
}

/// Waits until the worker is started, see Job::start_worker()
Job* ThreadPool::wait_for_job(Worker& worker, int spin_factor) const
{
  Job* job = nullptr;

  wait_until(worker.sleeper, spin_factor, [&]()
  {
    job = worker.job.load(std::memory_order_seq_cst);
    return job != nullptr;
  });

  worker.job.store(nullptr, std::memory_order_relaxed);
  return job;
}

/// The calling thread is the job's main thread and runs thread 0
void ThreadPool::run(int threads, FunctionRef func)
{
  WorkerGroup group = borrow(threads - 1);
  INDETERMINATE Job job(func, threads, group);
  job.start();

  if (group.new_workers + group.temp_workers > 0)
    create_workers(group, job);

  job.run_thread(0);
  job.run_remaining_threads();

  wait_for_workers(group, job);
  give_back(group);
  job.rethrow_exception();
}

/// Borrows up to the requested number of idle workers, reserves
/// new workers if the pool may grow and temporary workers for
/// the rest.
///
WorkerGroup ThreadPool::borrow(int workers)
{
  // Creating a MainThreadSync may enter the kernel (malloc,
  // page faults), hence we create it before taking the lock
  // if no unused MainThreadSync seems to be available.
  std::unique_ptr<MainThreadSync> new_sync;
  if (!unused_syncs_.load(std::memory_order_relaxed))
    new_sync.reset(new MainThreadSync(max_size_, workers));

  WorkerGroup group;
  std::lock_guard<std::mutex> lock(mutex_);
  MainThreadSync* sync = unused_syncs_.load(std::memory_order_relaxed);

  if (!sync)
  {
    // Rare race: after our check before taking the lock, another
    // main thread took the last unused MainThreadSync.
    if_unlikely(!new_sync)
      new_sync.reset(new MainThreadSync(max_size_, workers));
    sync = new_sync.release();
  }
  else if (new_sync)
  {
    // Another main thread returned its MainThreadSync meanwhile,
    // we use it and keep the new one for a later job.
    new_sync->next = sync->next;
    sync->next = new_sync.release();
  }

  group.sync = sync;
  unused_syncs_.store(sync->next, std::memory_order_relaxed);
  group.size = std::min(workers, (int) idle_workers_.size());

  // Pop workers from the idle stack. Workers read the group's
  // worker array, hence we only write entries that changed in
  // order not to move cache lines back to the main thread.
  Worker** idle_top = idle_workers_.end();
  for (int i = 0; i < group.size; i++)
  {
    Worker* worker = *--idle_top;
    if (i >= group.sync->initialized ||
        group.sync->workers[i] != worker)
      group.sync->workers[i] = worker;
  }

  idle_workers_.resize(idle_workers_.size() - group.size);

  int missing = workers - group.size;
  int room = max_size_ - size_;
  group.new_workers = std::min(missing, room);
  group.temp_workers = missing - group.new_workers;
  size_ += group.new_workers;

  // The busy workers and a main thread exceed the CPUs
  int busy = busy_workers_ + workers;
  if (busy >= cpus_ && busy_workers_ < cpus_)
    oversubscribed_.store(true, std::memory_order_relaxed);
  busy_workers_ = busy;

  return group;
}

/// Creates the new workers, they start working on the job
/// immediately. Creating threads is slow, hence we create them
/// without holding the pool lock, the job's work proceeds and
/// waiting threads keep spinning meanwhile. The workers already
/// use the job, hence a failure cannot be thrown and we abort,
/// like OpenMP runtimes do.
///
void ThreadPool::create_workers(WorkerGroup& group, Job& job) noexcept
{
  creating_workers_.fetch_add(1, std::memory_order_relaxed);

  try
  {
    for (int i = 0; i < group.new_workers; i++)
    {
      std::unique_ptr<Worker> worker(new Worker);
      std::thread(&ThreadPool::worker_loop, this, worker.get(), &job).detach();
      group.sync->workers[group.size++] = worker.release();
    }

    for (int i = 0; i < group.temp_workers; i++)
      std::thread(&ThreadPool::worker_loop, this, nullptr, &job).detach();
  }
  catch (const std::exception& e)
  {
    std::cerr << "primecount: Failed to create thread: " << e.what() << std::endl;
    std::abort();
  }

  creating_workers_.fetch_sub(1, std::memory_order_relaxed);
}

/// Waits until all workers of the job have finished
void ThreadPool::wait_for_workers(const WorkerGroup& group, const Job& job)
{
  waiting_main_threads_.fetch_add(1, std::memory_order_relaxed);

  wait_until(group.sync->sleeper, get_spin_factor(job.threads()), [&]()
  {
    return job.finished();
  });

  waiting_main_threads_.fetch_sub(1, std::memory_order_relaxed);
}

/// Returns the borrowed workers to the pool. The most recently
/// used workers are borrowed first next time because they are
/// the most likely to still spin. Pushing them in reverse order
/// gives them the same positions in the next job.
///
void ThreadPool::give_back(WorkerGroup& group)
{
  std::lock_guard<std::mutex> lock(mutex_);

  std::size_t old_size = idle_workers_.size();
  ASSERT(old_size + group.size <= idle_workers_.capacity());
  idle_workers_.resize(old_size + group.size);

  std::reverse_copy(group.sync->workers.begin(),
                    group.sync->workers.begin() + group.size,
                    idle_workers_.begin() + old_size);

  group.sync->initialized = std::max(group.sync->initialized, group.size);
  group.sync->next = unused_syncs_.load(std::memory_order_relaxed);
  unused_syncs_.store(group.sync, std::memory_order_relaxed);

  int busy = busy_workers_ - (group.size + group.temp_workers);
  if (busy < cpus_ && busy_workers_ >= cpus_)
    oversubscribed_.store(false, std::memory_order_relaxed);
  busy_workers_ = busy;
}

int ThreadPool::size()
{
  std::lock_guard<std::mutex> lock(mutex_);
  return size_;
}

void ThreadPool::set_max_spin_time(double ms)
{
  double us = std::min(ms * 1e3, (double) INT_MAX);
  max_spin_time_us_.store((int) std::max(us, 0.0), std::memory_order_relaxed);
}

/// Starting a job wakes up the levels of its binomial tree one
/// after another, which takes about tree_depth times the wake-up
/// latency. Wake-up latencies vary a lot, hence waiting threads
/// busy-wait 4x longer, see spin_deadline().
///
int ThreadPool::get_spin_factor(int threads)
{
  int tree_depth = ilog2(std::max(threads - 1, 1)) + 1;
  return 4 * tree_depth;
}

/// Waiting threads busy-wait for spin_factor times
/// their measured wake-up latency, but at most for
/// max_spin_time. Then they sleep.
///
Clock::time_point ThreadPool::spin_deadline(const Sleeper& sleeper, int spin_factor) const
{
  int64_t max_ns = max_spin_time_us_.load(std::memory_order_relaxed) * 1000ll;
  int64_t ns = std::min(spin_factor * sleeper.wake_ns(), max_ns);
  return Clock::now() + std::chrono::nanoseconds(ns);
}

#if defined(HAVE_GET_CPUS_WIN32)

/// Returns the number of CPUs that the process may use. Before
/// Windows 11, or if restricted e.g. by start /affinity, the
/// process is assigned to a single processor group and may use
/// the CPUs of its affinity mask. On Windows 11 a process spans
/// all processor groups by default, then its affinity mask is
/// unreliable (LLVM D138747) and we count the CPUs of all groups.
/// .NET and LLVM use the same method.
///
int get_cpus_win32()
{
  USHORT groups = 0;

  // Fails with ERROR_INSUFFICIENT_BUFFER and sets groups
  if (GetProcessGroupAffinity(GetCurrentProcess(), &groups, nullptr) ||
      GetLastError() != ERROR_INSUFFICIENT_BUFFER)
    groups = 1;

  DWORD_PTR process_mask = 0;
  DWORD_PTR system_mask = 0;

  if (groups == 1 &&
      GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask) &&
      process_mask != 0)
    return (int) popcnt64(process_mask);

  // As of 11/2026, std::thread::hardware_concurrency() returns
  // at most 64 (one processor group) with MSYS2/MinGW-w64,
  // hence we count the CPUs of all processor groups.
  return (int) GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
}

#endif

} // namespace

namespace primecount {

void parallel_for_impl(int threads, FunctionRef func)
{
  ASSERT(threads > 1);
  get_pool().run(threads, func);
}

int max_threads()
{
  static const int threads = []() -> int
  {
#if defined(HAVE_GET_CPUS_WIN32)
    int n = get_cpus_win32();
#else
    int n = (int) std::thread::hardware_concurrency();
#endif

#if defined(HAVE_SCHED_GETAFFINITY)
    // Respect the CPU affinity mask, e.g. taskset
    cpu_set_t cpus;
    if (sched_getaffinity(0, sizeof(cpus), &cpus) == 0)
      n = CPU_COUNT(&cpus);
#endif

    return std::max(n, 1);
  }();

  return threads;
}

int thread_pool_size()
{
  return get_pool().size();
}

void set_max_thread_spin_time(double ms)
{
  get_pool().set_max_spin_time(ms);
}

} // namespace

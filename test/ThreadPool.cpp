///
/// @file   ThreadPool.cpp
/// @brief  Test parallel_for(), parallel_sum() and
///         primecount's thread pool.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <threadpool/parallel.hpp>
#include <int128_t.hpp>
#include <Vector.hpp>

#include <stdint.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

using namespace primecount;

std::atomic<long> progress{0};

void check(bool OK, const std::string& msg)
{
  progress++;
  std::cout << msg << "   " << (OK ? "OK" : "ERROR") << std::endl;
  if (!OK)
    std::exit(1);
}

/// Unlike sleep_for(), busy-waiting is precise on all operating
/// systems, e.g. sleep_for(1ms) may take 15.6 ms on Windows.
///
void busy_wait(std::chrono::microseconds duration)
{
  auto stop = std::chrono::steady_clock::now() + duration;
  while (std::chrono::steady_clock::now() < stop)
    std::this_thread::yield();
}

/// Fails the test if a parallel job hangs
void start_watchdog()
{
  std::thread([]()
  {
    long last = -1;

    while (true)
    {
      std::this_thread::sleep_for(std::chrono::seconds(60));
      long current = progress.load();
      if (current == last)
      {
        std::cout << "No progress for 60 seconds, parallel job hangs   ERROR" << std::endl;
        std::_Exit(1);
      }
      last = current;
    }
  }).detach();
}

/// Checks that every id in [0, threads[ runs exactly once
/// and that the calling thread executes id 0.
///
bool check_ids(int threads)
{
  progress++;
  Vector<int> counts(threads);
  std::fill(counts.begin(), counts.end(), 0);
  std::atomic<int> total{0};
  std::thread::id main_thread = std::this_thread::get_id();
  bool main_runs_id0 = false;

  parallel_for(threads, [&](int id)
  {
    counts[id]++;
    total.fetch_add(1, std::memory_order_relaxed);
    if (id == 0)
      main_runs_id0 = (std::this_thread::get_id() == main_thread);
  });

  bool OK = (total == threads) && main_runs_id0;
  for (int i = 0; i < threads; i++)
    OK = OK && (counts[i] == 1);

  return OK;
}

int64_t sum_ids(int threads)
{
  return parallel_sum<int64_t>(threads, [](int id) -> int64_t { return id; });
}

int64_t expected_sum_ids(int threads)
{
  return (int64_t) threads * (threads - 1) / 2;
}

/// The exception must be rethrown after all ids have finished
void check_exception(int threads, int throw_id)
{
  std::atomic<int> finished{0};
  bool caught = false;
  int finished_at_catch = -1;

  try
  {
    parallel_for(threads, [&](int id)
    {
      if (id == throw_id)
        throw std::runtime_error("parallel_for test exception");
      busy_wait(std::chrono::microseconds(1000));
      finished++;
    });
  }
  catch (const std::runtime_error&)
  {
    caught = true;
    finished_at_catch = finished;
  }

  check(caught && finished_at_catch == threads - 1,
        "Exception from id " + std::to_string(throw_id) + " of " + std::to_string(threads) + " threads");
}

/// Jobs in which every id waits until all ids have started.
/// These hang if a worker misses a wake-up or if a job runs
/// with fewer threads than ids.
///
bool barrier_jobs(int threads, int jobs)
{
  std::atomic<int> arrived{0};

  for (int i = 0; i < jobs; i++)
  {
    progress++;

    parallel_for(threads, [&]()
    {
      arrived.fetch_add(1);
      // Relaxed loads, with ThreadSanitizer many threads that
      // acquire-load the same atomic can starve its writers.
      while (arrived.load(std::memory_order_relaxed) < threads)
        std::this_thread::yield();
    });

    if (arrived.load() != threads)
      return false;

    arrived = 0;
  }

  return true;
}

/// Nested jobs whose threads all wait for each other. The outer
/// job has threads ids and each inner job has inner ids, hence
/// threads * inner threads run at once and the pool is full if
/// threads * inner >= 2 * max_threads(). Every job has at most
/// max_threads() ids.
///
bool nested_barrier_job(int threads, int inner)
{
  std::atomic<int> arrived{0};
  int total = threads * inner;

  parallel_for(threads, [&, inner, total]()
  {
    parallel_for(inner, [&, total]()
    {
      arrived.fetch_add(1);
      while (arrived.load(std::memory_order_relaxed) < total)
        std::this_thread::yield();
    });
  });

  return arrived.load() == total;
}

int main()
{
  start_watchdog();
  int cpus = max_threads();
  // Jobs with 1 thread never create workers
  int max_workers = (cpus >= 2) ? 2 * cpus - 1 : 0;
  // Nested jobs run up to 4 * nested_cpus threads at once. Each
  // thread reserves 8 MiB of stack on Linux, hence a 32-bit
  // process fits only a few hundred threads.
  int nested_cpus = std::min(cpus, 48);
  // nested_barrier_job() fills the pool if it needs at least max_workers
  bool fills_pool = nested_cpus * std::min(nested_cpus, 3) - 1 >= max_workers;

  {
    int calls = 0;
    parallel_for(1, [&](int id) { calls += (id == 0); });
    parallel_for(0, [&](int id) { calls += (id == 0); });
    int64_t sum = parallel_sum<int64_t>(1, [](int id) -> int64_t { return id + 7; });
    bool OK = calls == 2 && sum == 7;
    OK = OK && thread_pool_size() == 0;
    check(OK, "threads <= 1 runs inline without creating threads");
  }

  {
    bool OK = true;
    for (int threads = 1; threads <= cpus; threads++)
      OK = OK && check_ids(threads);
    check(OK, "Every id runs exactly once for threads = 1 to " + std::to_string(cpus));

    // Jobs whose threads finish early may skip creating workers,
    // a job whose threads wait for each other creates all of them.
    OK = thread_pool_size() <= cpus - 1;
    OK = OK && barrier_jobs(cpus, 1);
    check(OK && thread_pool_size() == cpus - 1, "Thread pool size = " + std::to_string(thread_pool_size()));
  }

  {
    int a = std::min(cpus, 4);
    int b = std::min(cpus, 5);
    int c = std::min(cpus, 3);
    Vector<int> counts2(a * b);
    Vector<int> counts3(a * b * c);
    std::fill(counts2.begin(), counts2.end(), 0);
    std::fill(counts3.begin(), counts3.end(), 0);

    parallel_for(a, [&](int i) {
      parallel_for(b, [&](int j) {
        counts2[i * b + j]++;
      });
    });

    parallel_for(a, [&](int i) {
      parallel_for(b, [&](int j) {
        parallel_for(c, [&](int k) {
          counts3[(i * b + j) * c + k]++;
        });
      });
    });

    bool OK = std::all_of(counts2.begin(), counts2.end(), [](int n) { return n == 1; });
    check(OK, "Nested parallel jobs with 2 levels");
    OK = std::all_of(counts3.begin(), counts3.end(), [](int n) { return n == 1; });
    check(OK, "Nested parallel jobs with 3 levels");

    // The outer job uses all CPUs (at most 48),
    // hence the inner jobs oversubscribe the CPUs.
    int inner = std::min(nested_cpus, 4);
    int64_t sum = parallel_sum<int64_t>(nested_cpus, [inner]() -> int64_t { return sum_ids(inner); });
    check(sum == nested_cpus * expected_sum_ids(inner), "Nested parallel jobs with more threads than CPUs");
  }

  {
    // Like OpenMP, a job runs with all its threads, even if
    // the pool is full or other jobs already use all CPUs.
    bool OK = nested_barrier_job(nested_cpus, std::min(nested_cpus, 3));
    std::atomic<bool> nested_OK{true};
    int threads = std::min(cpus, 8);

    parallel_for(std::min(cpus, 3), [&, threads]()
    {
      if (!barrier_jobs(threads, 10))
        nested_OK = false;
    });

    OK = OK && nested_OK;
    check(OK, "Oversubscribed jobs in which all threads wait for each other");
    int size = thread_pool_size();
    check(size == max_workers || (!fills_pool && size < max_workers),
          "Thread pool size = " + std::to_string(size));
  }

  {
    int threads = std::min(cpus, 8);
    check_exception(std::min(cpus, 2), 0);
    check_exception(std::min(cpus, 2), std::min(cpus, 2) - 1);
    check_exception(threads, 0);
    check_exception(threads, threads - 1);

    int n = std::min(cpus, 4);
    bool caught = false;
    try
    {
      parallel_for(n, [n]() {
        parallel_for(n, [n](int id) {
          if (id == n - 1)
            throw std::runtime_error("nested exception");
        });
      });
    }
    catch (const std::runtime_error&)
    {
      caught = true;
    }
    check(caught, "Exception from a nested parallel job");
  }

  {
    std::atomic<bool> OK{true};
    std::thread threads[8];

    for (int t = 0; t < 8; t++)
    {
      threads[t] = std::thread([&OK, t, cpus]()
      {
        int m = std::min(cpus, 3);

        for (int i = 0; i < 10; i++)
        {
          int n = std::min(cpus, 2 + (t + i) % 3);
          if (sum_ids(n) != expected_sum_ids(n))
            OK = false;

          int64_t sum = parallel_sum<int64_t>(m, [m]() -> int64_t { return sum_ids(m); });
          if (sum != m * expected_sum_ids(m))
            OK = false;
        }
      });
    }

    for (auto& thread : threads)
      thread.join();

    // The pool never grows beyond max_workers
    int size = thread_pool_size();
    OK = OK && (size == max_workers || (!fills_pool && size < max_workers));
    check(OK, "Parallel jobs from 8 std::threads at once");
  }

  {
    bool OK = true;
    int threads = std::min(cpus, 16);
    set_max_thread_spin_time(0.2);

    for (int gap_us : {0, 100, 300, 1000})
    {
      for (int i = 0; i < 10; i++)
      {
        busy_wait(std::chrono::microseconds(gap_us));
        OK = OK && sum_ids(threads) == expected_sum_ids(threads);
        OK = OK && check_ids(cpus);
      }
    }

    // Workers sleep after each job
    set_max_thread_spin_time(0);
    for (int i = 0; i < 500; i++)
      OK = OK && sum_ids(threads) == expected_sum_ids(threads);

    for (int i = 0; i < 3; i++)
    {
      busy_wait(std::chrono::microseconds(1000));
      OK = OK && check_ids(cpus);
    }

    set_max_thread_spin_time(3);
    check(OK, "Parallel jobs with idle times shorter and longer than spin_time");
  }

  // Jobs with more threads than CPUs oversubscribe the CPUs,
  // then waiting threads sleep, which is slow for tiny jobs.
  if (cpus >= 2)
  {
    int max_job_threads = std::min(cpus, 8);
    int64_t sum = 0;
    int64_t expected = 0;
    for (int i = 0; i < 10000; i++)
    {
      if (i % 1000 == 0)
        progress++;
      int threads = 2 + i % (max_job_threads - 1);
      sum += sum_ids(threads);
      expected += expected_sum_ids(threads);
    }
    check(sum == expected, "Stress test with 10,000 tiny parallel jobs");
  }

  if (cpus >= 2)
  {
    int threads = std::min(cpus, 4);
    bool OK = barrier_jobs(threads, 4000);
    set_max_thread_spin_time(0);
    OK = OK && barrier_jobs(threads, 1000);
    set_max_thread_spin_time(3);
    check(OK, "4,000 + 1,000 jobs in which " + std::to_string(threads) + " threads wait for each other");
  }

  if (cpus >= 64)
    check(barrier_jobs(64, 2000), "2,000 jobs in which 64 threads wait for each other");

  {
    int threads = std::min(cpus, 64);
    int64_t s1 = parallel_sum<int64_t>(threads, [](int id) -> int64_t { return -1000000000000 * (id + 1); });
    uint64_t s2 = parallel_sum<uint64_t>(threads, [](int id) -> uint64_t { return (1ull << 57) + id; });
    maxint_t s3 = parallel_sum<maxint_t>(threads, [](int id) -> maxint_t { return (maxint_t) id * 3; });
    int64_t e1 = 0;
    uint64_t e2 = 0;
    maxint_t e3 = 0;

    for (int id = 0; id < threads; id++)
    {
      e1 += -1000000000000 * (id + 1);
      e2 += (1ull << 57) + id;
      e3 += (maxint_t) id * 3;
    }

    bool OK = (s1 == e1) && (s2 == e2) && (s3 == e3);

#if defined(HAVE_INT128_T)
    int128_t s4 = parallel_sum<int128_t>(threads, [](int id) -> int128_t { return -(int128_t(1) << 100) + id; });
    uint128_t s5 = parallel_sum<uint128_t>(threads, [](int id) -> uint128_t { return (uint128_t(1) << 90) + id; });
    int128_t e4 = 0;
    uint128_t e5 = 0;

    for (int id = 0; id < threads; id++)
    {
      e4 += -(int128_t(1) << 100) + id;
      e5 += (uint128_t(1) << 90) + id;
    }

    OK = OK && (s4 == e4) && (s5 == e5);
#endif

    // More ids than partial sums on the stack
    if (cpus > 128)
      OK = OK && sum_ids(cpus) == expected_sum_ids(cpus);

    // Dynamically claimed chunks of 16 iterations
    int64_t iters = 100000;
    std::atomic<int64_t> next_chunk{0};
    int64_t s6 = parallel_sum<int64_t>(threads, [&]() -> int64_t
    {
      int64_t sum = 0;
      for (int64_t i = next_chunk.fetch_add(16, std::memory_order_relaxed);
           i < iters;
           i = next_chunk.fetch_add(16, std::memory_order_relaxed))
        for (int64_t j = i; j < std::min(i + 16, iters); j++)
          sum += j;
      return sum;
    });

    OK = OK && s6 == iters * (iters - 1) / 2;

    bool caught = false;
    try
    {
      parallel_sum<int64_t>(threads, [](int id) -> int64_t
      {
        if (id == 1)
          throw std::runtime_error("parallel_sum test exception");
        return id;
      });
    }
    catch (const std::runtime_error&)
    {
      caught = true;
    }

    OK = OK && (caught || threads < 2);
    check(OK, "parallel_sum() with 64-bit and 128-bit integers");
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

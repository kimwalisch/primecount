///
/// @file   parallel.cpp
/// @brief  Test parallel_region(), parallel_sum() and
///         primecount's thread pool.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <parallel.hpp>
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

#if defined(__has_feature)
  #if __has_feature(thread_sanitizer)
    #define THREAD_SANITIZER
  #endif
#endif

#if defined(__SANITIZE_THREAD__)
  #define THREAD_SANITIZER
#endif

#if (defined(__unix__) || defined(__APPLE__)) && \
    !defined(__EMSCRIPTEN__) && \
    !defined(THREAD_SANITIZER)
  #include <sys/wait.h>
  #include <unistd.h>
  #define TEST_FORK
#endif

using namespace primecount;

std::atomic<long> progress{0};

void check(bool OK, const std::string& msg)
{
  progress++;
  std::cout << msg << "   " << (OK ? "OK" : "ERROR") << std::endl;
  if (!OK)
    std::exit(1);
}

/// Fails the test if a parallel region hangs
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
        std::cout << "No progress for 60 seconds, parallel region hangs   ERROR" << std::endl;
        std::_Exit(1);
      }
      last = current;
    }
  }).detach();
}

/// Checks that every id in [0, threads) runs exactly once
/// and that the calling thread executes id 0.
///
bool check_ids(int threads)
{
  progress++;
  Vector<int> counts(threads);
  std::fill(counts.begin(), counts.end(), 0);
  std::atomic<int> total{0};
  std::thread::id leader = std::this_thread::get_id();
  bool leader_id0 = false;

  parallel_region(threads, [&](int id)
  {
    counts[id]++;
    total.fetch_add(1, std::memory_order_relaxed);
    if (id == 0)
      leader_id0 = (std::this_thread::get_id() == leader);
  });

  bool OK = (total == threads) && leader_id0;
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
    parallel_region(threads, [&](int id)
    {
      if (id == throw_id)
        throw std::runtime_error("parallel_region test exception");
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
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

/// Regions in which every id waits until all ids have started.
/// These hang if a worker misses a wake-up.
///
bool barrier_regions(int threads, int regions)
{
  std::atomic<int> arrived{0};

  for (int i = 0; i < regions; i++)
  {
    progress++;

    parallel_region(threads, [&](int)
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

int main()
{
  start_watchdog();
  int cpus = max_threads();
  int max_workers = cpus - 1;

  {
    int calls = 0;
    parallel_region(1, [&](int id) { calls += (id == 0); });
    parallel_region(0, [&](int id) { calls += (id == 0); });
    int64_t sum = parallel_sum<int64_t>(1, [](int id) -> int64_t { return id + 7; });
    bool OK = calls == 2 && sum == 7;
    OK = OK && parallel_pool_size() == 0;
    check(OK, "threads <= 1 runs inline without creating threads");
  }

  {
    bool OK = true;
    for (int threads = 1; threads <= 2 * cpus; threads++)
      OK = OK && check_ids(threads);
    OK = OK && check_ids(256);
    check(OK, "Every id runs exactly once for threads = 1 to " + std::to_string(2 * cpus) + " and 256");
    check(parallel_pool_size() == max_workers, "Thread pool size = " + std::to_string(parallel_pool_size()));
  }

  {
    int a = 4, b = 5, c = 3;
    Vector<int> counts2(a * b);
    Vector<int> counts3(a * b * c);
    std::fill(counts2.begin(), counts2.end(), 0);
    std::fill(counts3.begin(), counts3.end(), 0);

    parallel_region(a, [&](int i) {
      parallel_region(b, [&](int j) {
        counts2[i * b + j]++;
      });
    });

    parallel_region(a, [&](int i) {
      parallel_region(b, [&](int j) {
        parallel_region(c, [&](int k) {
          counts3[(i * b + j) * c + k]++;
        });
      });
    });

    bool OK = std::all_of(counts2.begin(), counts2.end(), [](int n) { return n == 1; });
    check(OK, "Nested parallel regions with 2 levels");
    OK = std::all_of(counts3.begin(), counts3.end(), [](int n) { return n == 1; });
    check(OK, "Nested parallel regions with 3 levels");

    // All workers are busy in the outer region, hence
    // the inner regions run without workers.
    int64_t sum = parallel_sum<int64_t>(cpus, [](int) -> int64_t { return sum_ids(4); });
    check(sum == cpus * expected_sum_ids(4), "Nested parallel regions without idle workers");
  }

  {
    check_exception(2, 0);
    check_exception(2, 1);
    check_exception(8, 0);
    check_exception(8, 7);

    bool caught = false;
    try
    {
      parallel_region(4, [](int) {
        parallel_region(4, [](int id) {
          if (id == 3)
            throw std::runtime_error("nested exception");
        });
      });
    }
    catch (const std::runtime_error&)
    {
      caught = true;
    }
    check(caught, "Exception from a nested parallel region");
  }

  {
    std::atomic<bool> OK{true};
    std::thread threads[8];

    for (int t = 0; t < 8; t++)
    {
      threads[t] = std::thread([&OK, t]()
      {
        for (int i = 0; i < 300; i++)
        {
          int n = 2 + (t + i) % 15;
          if (sum_ids(n) != expected_sum_ids(n))
            OK = false;

          int64_t sum = parallel_sum<int64_t>(3, [](int) -> int64_t { return sum_ids(3); });
          if (sum != 3 * expected_sum_ids(3))
            OK = false;
        }
      });
    }

    for (auto& thread : threads)
      thread.join();

    check(OK, "Parallel regions from 8 std::threads at once");
  }

  {
    bool OK = true;
    int threads = std::min(cpus, 16);
    set_parallel_spin_time(0.001);

    for (int gap_us : {0, 200, 5000})
    {
      for (int i = 0; i < 50; i++)
      {
        std::this_thread::sleep_for(std::chrono::microseconds(gap_us));
        OK = OK && sum_ids(threads) == expected_sum_ids(threads);
        OK = OK && check_ids(cpus);
      }
    }

    // Workers sleep after each region
    set_parallel_spin_time(0);
    for (int i = 0; i < 2000; i++)
      OK = OK && sum_ids(threads) == expected_sum_ids(threads);

    for (int i = 0; i < 3; i++)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      OK = OK && check_ids(cpus);
    }

    set_parallel_spin_time(0.003);
    check(OK, "Parallel regions with idle times shorter and longer than spin_time");
  }

  {
    int64_t sum = 0;
    int64_t expected = 0;
    for (int i = 0; i < 100000; i++)
    {
      if (i % 10000 == 0)
        progress++;
      int threads = 2 + i % 7;
      sum += sum_ids(threads);
      expected += expected_sum_ids(threads);
    }
    check(sum == expected, "Stress test with 100,000 tiny parallel regions");
  }

  if (cpus >= 2)
  {
    int threads = std::min(cpus, 4);
    bool OK = barrier_regions(threads, 20000);
    set_parallel_spin_time(0);
    OK = OK && barrier_regions(threads, 5000);
    set_parallel_spin_time(0.003);
    check(OK, "20,000 + 5,000 regions in which " + std::to_string(threads) + " threads wait for each other");
  }

  if (cpus >= 64)
    check(barrier_regions(64, 200000), "200,000 regions in which 64 threads wait for each other");

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
    OK = OK && sum_ids(300) == expected_sum_ids(300);

    // Dynamically claimed chunks of 16 iterations
    int64_t iters = 100000;
    std::atomic<int64_t> next_chunk{0};
    int64_t s6 = parallel_sum<int64_t>(threads, [&](int) -> int64_t
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

#if defined(TEST_FORK)
  {
    pid_t pid = fork();

    // The child process does not contain the parent's
    // worker threads, it must create new ones.
    if (pid == 0)
    {
      int threads = std::min(cpus, 16);
      bool OK = sum_ids(threads) == expected_sum_ids(threads);
      OK = OK && check_ids(threads);
      std::_Exit(OK ? 0 : 1);
    }

    int status = 0;
    bool OK = pid > 0 &&
              waitpid(pid, &status, 0) == pid &&
              WIFEXITED(status) &&
              WEXITSTATUS(status) == 0;

    check(OK, "Parallel regions in a child process after fork()");
  }
#endif

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}

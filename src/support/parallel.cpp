// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/support/parallel.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace quicksynth {

namespace {

class ThreadPool {
  public:
    static ThreadPool& instance();

    // Runs fn() on `threads` workers, caller included; returns when all finish.
    void run(int threads, const std::function<void()>& fn);

  private:
    ThreadPool() = default;
    void grow(size_t want);
    void worker();

    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable start_, done_;
    const std::function<void()>* job_ = nullptr;
    unsigned long generation_ = 0;
    size_t pending_ = 0, wanted_ = 0, taken_ = 0;
    bool busy_ = false;
};

// Leaked on purpose: workers park for the process lifetime, and joining them
// in a static destructor at interpreter shutdown can hang.
ThreadPool& ThreadPool::instance() {
    static ThreadPool* p = new ThreadPool();
    return *p;
}

void ThreadPool::worker() {
    unsigned long seen = 0;
    for (;;) {
        const std::function<void()>* job = nullptr;
        {
            std::unique_lock<std::mutex> lock(m_);
            start_.wait(lock, [&] { return generation_ != seen && taken_ < wanted_; });
            seen = generation_;
            ++taken_;
            job = job_;
        }
        (*job)();
        {
            std::lock_guard<std::mutex> lock(m_);
            --pending_;
        }
        done_.notify_one();
    }
}

void ThreadPool::grow(size_t want) {
    while (workers_.size() < want) workers_.emplace_back([this] { worker(); });
}

void ThreadPool::run(int threads, const std::function<void()>& fn) {
    if (threads <= 1) {
        fn();
        return;
    }
    const size_t helpers = static_cast<size_t>(threads) - 1;
    {
        std::unique_lock<std::mutex> lock(m_);
        if (busy_) {  // nested or concurrent: run on the caller
            lock.unlock();
            fn();
            return;
        }
        busy_ = true;
        grow(helpers);
        job_ = &fn;
        wanted_ = helpers;
        taken_ = 0;
        pending_ = helpers;
        ++generation_;
    }
    start_.notify_all();
    fn();  // the caller is a worker too
    {
        std::unique_lock<std::mutex> lock(m_);
        done_.wait(lock, [&] { return pending_ == 0; });
        job_ = nullptr;
        wanted_ = 0;
        busy_ = false;
    }
}

}  // namespace

void parallel_for_on_pool(size_t items, int threads, const std::function<void(size_t)>& fn) {
    std::atomic<size_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_lock;
    const std::function<void()> work = [&]() {
        while (!failed.load(std::memory_order_relaxed)) {
            const size_t k = next.fetch_add(1);
            if (k >= items) return;
            try {
                fn(k);
            } catch (...) {
                const std::lock_guard<std::mutex> hold(error_lock);
                if (!error) error = std::current_exception();
                failed.store(true, std::memory_order_relaxed);
            }
        }
    };
    ThreadPool::instance().run(static_cast<int>(std::min(static_cast<size_t>(threads), items)),
                               work);
    if (error) std::rethrow_exception(error);
}

}  // namespace quicksynth

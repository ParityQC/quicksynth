// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Parallel loops on a process-wide thread pool. Workers are created once and
// parked between loops, since a scan opens one per scale and determinant class.
#pragma once

#include <cstddef>
#include <functional>

namespace quicksynth {

// parallel_for's pool path; needs at least two items and two threads.
void parallel_for_on_pool(size_t items, int threads, const std::function<void(size_t)>& fn);

// Runs fn(0..items-1) on min(threads, items) workers, caller included; items go
// out in order via a shared counter. Nested or concurrent loops run serially, so
// fn must not depend on the worker count. On a throw, unstarted items are
// skipped and the first exception is rethrown after all workers return.
template <class Fn>
void parallel_for(size_t items, int threads, const Fn& fn) {
    if (threads <= 1 || items < 2) {
        for (size_t k = 0; k < items; ++k) fn(k);
        return;
    }
    parallel_for_on_pool(items, threads, std::cref(fn));
}

}  // namespace quicksynth

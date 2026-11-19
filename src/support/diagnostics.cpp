// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/support/diagnostics.hpp"

#include <algorithm>
#include <cstdio>

namespace quicksynth {

namespace {

std::string format(const char* fmt, double x) {
    char buf[64];
    std::snprintf(buf, sizeof buf, fmt, x);
    return buf;
}

}  // namespace

std::string walk_line(int K, int det_class, int64_t points, std::vector<int64_t> costs,
                      double estimated, double seconds) {
    std::sort(costs.begin(), costs.end());
    std::string hits;
    for (size_t i = 0; i < std::min<size_t>(costs.size(), 8); ++i)
        hits += (i ? ", " : "") + std::to_string(costs[i]);
    return "    K=" + std::to_string(K) + " parity " + std::to_string(det_class) + ": " +
           std::to_string(points) + " points, hits [" + hits + "]" +
           (costs.size() > 8 ? "..." : "") + "  est " + format("%.0f", estimated) + "  [" +
           format("%.1f", seconds) + "s]";
}

}  // namespace quicksynth

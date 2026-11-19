// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Replay helpers for golden_scan, golden_split and golden_combined. Replays
// are single-threaded and require the recorded element.
#pragma once

#include "fixtures.hpp"

#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/price_candidates/price.hpp"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace quicksynth::test {

// Uncalibrated throughput the fixtures were recorded at.
double recorded_throughput(const std::string& ring);
// Shared radius-6 prefix ball at the catalysis price.
std::shared_ptr<const CostBall> recorded_prefix_ball(const std::string& ring);

// Recorded run, expected result and target.
struct Recorded {
    const json& run;
    const json& want;
    Unitary V;
};
// Runs of golden_<group> on one ring and price with `key` == `value` (null if absent).
std::vector<Recorded> recorded_runs(const std::string& group, const std::string& ring,
                                    const Price& price, const std::string& key, const json& value);

Price run_price(const json& run);
GateSet run_gateset(const json& run);
// Cost under the run's gate set and price.
Cost run_cost(const json& run, const CyclotomicMatrix& e);
// Cost of a search hit; nullopt if none.
std::optional<int64_t> best_cost(const std::optional<Enumerator::Hit>& best);

// Checks `got` against `want` ({"element", "n_t", "n_q", "best_cost"}).
void check_same_circuit(const json& run, const json& want,
                        const std::optional<CyclotomicMatrix>& got);

}  // namespace quicksynth::test

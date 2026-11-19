// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "search_scales/replay.hpp"

#include "oracles.hpp"

#include <doctest/doctest.h>

#include <map>
#include <mutex>

namespace quicksynth::test {

double recorded_throughput(const std::string& ring) { return ring == "sqrtT" ? 1.4e7 : 2.5e7; }

std::shared_ptr<const CostBall> recorded_prefix_ball(const std::string& ring) {
    static std::mutex mutex;
    static std::map<std::string, std::shared_ptr<const CostBall>> balls;
    std::lock_guard<std::mutex> lock(mutex);
    std::shared_ptr<const CostBall>& slot = balls[ring];
    if (!slot) {
        auto ball = std::make_shared<CostBall>(Price{1, 3}, ring == "sqrtT");
        ball->build(6);
        slot = ball;
    }
    return slot;
}

std::vector<Recorded> recorded_runs(const std::string& group, const std::string& ring,
                                    const Price& price, const std::string& key, const json& value) {
    const json& fixture = golden(group);
    const std::vector<Unitary> targets = decode_targets(fixture["meta"]["targets"]);
    std::vector<Recorded> out;
    for (const json& record : fixture["runs"]) {
        const json& run = record["run"];
        if (run["ring"] == ring && run["price"] == json::array({price.t, price.q}) &&
            run.value(key, json()) == value)
            out.push_back({run, record["result"], targets[run["target"].get<size_t>()]});
    }
    return out;
}

Price run_price(const json& run) { return {run["price"][0], run["price"][1]}; }

GateSet run_gateset(const json& run) { return ring_gateset(run["ring"]); }

Cost run_cost(const json& run, const CyclotomicMatrix& e) {
    return Pricer::for_gateset(run_gateset(run), run_price(run)).price(e);
}

std::optional<int64_t> best_cost(const std::optional<Enumerator::Hit>& best) {
    if (!best) return std::nullopt;
    return best->cost.value;
}

void check_same_circuit(const json& run, const json& want,
                        const std::optional<CyclotomicMatrix>& got) {
    REQUIRE(!want["element"].is_null());
    REQUIRE(got.has_value());
    CHECK(row_of(*got) == want["element"]);
    const Cost c = run_cost(run, *got);
    CHECK(c.value == want["best_cost"]);
    CHECK(c.n_t == want["n_t"]);
    CHECK(c.n_q == want["n_q"]);
}

}  // namespace quicksynth::test

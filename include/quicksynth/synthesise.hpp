// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// One-call synthesis: a checked gate word within epsilon of a 2x2 unitary.
// The target is already replaced
// by its nearest unitary, and every search, distance and final check starts
// from its bits.
//
// Runs scan_then_split.hpp through stages() until one answers. On Clifford+T
// an answer at the scan's scale bound is a
// certified optimum. Every answer is spelt in the normal form, multiplied out
// exactly, and its distance, in double precision from the exact element
// (distance.hpp), checked against epsilon - kDistanceError; the searches run
// at epsilon - kEpsilonMargin. Long steps poll the host's
// interrupt check (host.hpp), and progress goes to the host's printer: under
// verbose a line for each walk of QuickSynth's scan and for each stage, and
// none while a cost ball is built.
//
// On Clifford+sqrt(T) a Clifford+T search then looks for a pure-T word
// cheaper than the answer: stage 0 and the stages, with the answer's cost as
// the incumbent, so each stops where its scale bound reaches that cost.
// A pure-T word lies far deeper in the sqrt(T) scales than a sqrt(T) word of
// its cost, so the sqrt(T) search stops before it; at a high sqrt(T) price
// it can be the cheaper. Its walks are not logged; under verbose one line
// reports a cheaper word.
//
// Stage 0 (Settings::rotation_path) comes first: a target Clifford-equivalent
// to a z-rotation is searched by the z-rotation path (synthesise_rotation/),
// on Clifford+T certified T-optimal when exact, on Clifford+sqrt(T) up to its
// count rule; it builds no prefix ball.
#pragma once

#include "quicksynth/policy.hpp"
#include "quicksynth/search_scales/scan_then_split.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/support/diagnostics.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth {

struct Synthesised {
    std::string word;
    Cost cost;
    double distance = 0.0;
    // Clifford+T: scan at its scale bound, or stage 0 at its bound on an exact
    // Clifford-equivalent z-rotation
    bool certified = false;
    int stages = 1;  // stages run; 0: stage 0 answered (of the Clifford+T search's, if it answered)
    // (stage, walk) of every walk, for time-model calibration
    std::vector<std::pair<int, SearchLog>> log;
    std::optional<SearchSummary> search;  // answering stage
};

// Thrown when the time limit passes before the next stage.
struct SynthesisTimeout : std::runtime_error {
    using std::runtime_error::runtime_error;
};

class QuickSynth {
  public:
    struct Settings {
        GateSet gateset = GateSet::sqrt_t();
        Price price;
        std::optional<int> extra_scales;     // scan scales past the first hit
        std::optional<double> point_budget;  // scan point budget, stage 1
        int workers = 1;
        std::optional<double> time_limit;  // seconds; no stage starts after it
        std::string time_limit_text;       // time_limit as given, for the message
        bool verbose = false;              // print progress
        // stage 1 (stopping_rule.hpp)
        StoppingRule stopping_rule = StoppingRule::fixed;
        // rate rule at kFastRate in the scan and the split, every stage, and
        // LevelSearch at kFineFastRate instead of kFineDefaultRate
        bool fast = false;
        // rate rule at this rate in the scan, the split and LevelSearch instead (over `fast`)
        std::optional<double> rate;
        bool partial_walk = false;  // PrefixSplit::Settings::partial
        // Clifford+sqrt(T), fixed stopping rule: finer levels (ScaleScan::Settings::fine_levels);
        // none: on
        std::optional<bool> fine_levels;
        // finer levels on Clifford+sqrt(T): LevelSearch; false: the fixed depth below
        bool level_search = true;
        bool skip_setup_levels = true;                  // LevelSearch::Settings
        std::optional<int> extra_levels, split_levels;  // none: the same depth as scales
        // Stage 0, the z-rotation path (synthesise_rotation/): a target within
        // rotation_tolerance * epsilon of a Clifford-equivalent z-rotation is
        // searched there first, at epsilon less its off-axis part; no answer
        // falls through to the stages. A Clifford+sqrt(T) word found before a
        // cap stopped the search is returned if the stages find nothing
        // cheaper or time_limit stops them. Counts against time_limit; stage 1
        // starts unless stage 0 left a word.
        bool rotation_path = true;
        double rotation_tolerance = 1e-2;
        // Clifford+sqrt(T) count rule: solvable diagonal entries before stage 0
        // stops (rotation_solvable()).
        std::optional<int> rotation_solvable;
    };

    explicit QuickSynth(Settings settings);

    // Stage 1 settings for the requested epsilon; `fallback` forces the fixed stopping rule.
    ScanThenSplit::Settings first_stage(double epsilon, bool fallback) const;

    // Cheapest word found within epsilon of U. Throws invalid_argument (price
    // never uses T), SynthesisTimeout, or runtime_error (no stage answers).
    Synthesised synthesise(const Unitary& U, double epsilon) const;

    // Builds all prefix balls and prefix sets ahead of timed calls.
    void prepare(double epsilon) const;

    // Stage 0's Clifford+sqrt(T) count rule: Settings::rotation_solvable, else
    // kRotationSolvableFast under fast, else kRotationSolvable.
    int rotation_solvable() const;

  private:
    // The stages and stage 0, without the Clifford+T comparison.
    Synthesised search(const Unitary& U, double epsilon) const;
    // Clifford+T: a word below `cost` from stage 0 and the stages, each with
    // `cost` as its incumbent; a stage that reaches it ends the search.
    std::optional<Synthesised> cheaper_than(const Unitary& U, double epsilon, int64_t cost,
                                            std::optional<double> deadline) const;
    std::optional<Synthesised> rotation_stage(const Unitary& U, double epsilon,
                                              std::optional<double> deadline,
                                              std::vector<std::pair<int, SearchLog>>& log,
                                              bool& final) const;
    std::shared_ptr<CostBall> stopping_ball(const ScanThenSplit::Settings& stage) const;
    std::shared_ptr<CostBall> stage_ball(double epsilon, int prefix_scale) const;

    Settings settings_;
    Price price_;  // price used by the searches and balls
};

// Process-wide LRU caches of prefix balls (eight, by price and radius) and
// prefix radii (32, by price and scale); prefix sets stay with their ball.
// Thread-safe: one lock, held while a ball is built.

// Prefix ball built to `radius` levels; the prefix sets depend on it.
std::shared_ptr<CostBall> prefix_ball(const Price& price, bool use_q, int radius);

// Least radius >= kPrefixRadius whose ball holds kMaxPrefixes elements of
// `scale`; at most 16, or 2 * scale on Clifford+T (no scale-4 element beyond
// T-count 8, none of scale 5 beyond 10; radius 8 holds all 384 of scale 4).
// Price-dependent: radius 6 holds 4800 scale-3 elements at the gate count, 288
// at catalysis, 0 at power.
int prefix_radius(const Price& price, bool use_q, int scale);

// kPrefixRadius, or prefix_radius() below kDeepEpsilon (Clifford+sqrt(T),
// scale 3: 6 gate count, 7 catalysis, 10 power; Clifford+T: 8 at scale 4, 9 at 5).
// Clifford+T below kDeeperEpsilon: every coset of the scale, T-count <= 2 scale
// (std::invalid_argument if 2 scale t_cost passes int).
int split_radius(const Price& price, const GateSet& gateset, double epsilon, int scale);

}  // namespace quicksynth

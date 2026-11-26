// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Scale scan [paper: alg:synth]. Walks every determinant class of each scale
// from the first that can hold an element within epsilon. The stopping rule sets how
// far past the first hit: extra_scales plus a predicted saving of kScaleSaving
// x best cost ("fixed"), or the counting model of stopping_rule.hpp ("model").
#pragma once

#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"
#include "quicksynth/support/diagnostics.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

namespace quicksynth {

constexpr int kExtraScales = 2;
constexpr double kPointBudget = 3e6;

class ScaleScan {
  public:
    struct Settings {
        int extra_scales = kExtraScales;     // scales walked after the first hit (fixed)
        double point_budget = kPointBudget;  // max estimated points per class walk
        // rate rule, magic states per second; none: off
        std::optional<double> rate = kDefaultRate;
        // fixed stopping rule: estimated points per second; none: reference_throughput.
        // The model stopping rule uses time_model instead.
        std::optional<double> throughput;
        StoppingRule stopping_rule = StoppingRule::fixed;
        // fixed stopping rule: walk finer levels (SearchRegion::first_level), L per scale
        bool fine_levels = false;
        std::optional<int> extra_levels;  // levels past the first hit; none: L extra_scales
        // cost to beat: stop where the scale bound reaches it, uncertified
        std::optional<int64_t> incumbent;
    };
    // Called after each walk with the scale's estimated points.
    using Progress = std::function<void(const Enumerator::Result&, double)>;
    struct Result {
        std::optional<Enumerator::Hit> best;  // first of the cheapest hits
        std::optional<int> first_hit_scale;
        bool certified = false;
        double lower_bound = 0.0;
        std::optional<int> last_full_scale;  // deepest scale walked in full
        // in levels (L per scale; L = 1 without finer levels)
        std::optional<int> first_hit_level, last_level;
        // stopped where the scale bound reached Settings::incumbent
        bool incumbent_reached = false;
        double estimated_points = 0.0;
        std::vector<SearchLog> log;  // one entry per walk
    };

    ScaleScan(double epsilon, GateSet gateset, Price price, int threads, Settings settings);

    const Settings& settings() const { return settings_; }
    Enumerator& enumerator() { return enumerator_; }
    // null under the fixed stopping rule
    const StoppingModel* model() const { return model_.get(); }
    // levels per scale: L with finer levels, else 1
    int levels() const;

    // Cost lower bound at sqrt2-scale K: integer_scale_bound (model stopping rule) or
    // scan_scale_bound (fixed).
    double scale_lower_bound(int K) const {
        return bound(K, settings_.stopping_rule == StoppingRule::model);
    }

    Result run(const Unitary& V, const Progress& progress = {});

    // Walks both determinant classes of level l (`est` points per class) into
    // `result`: hits, log, estimated points, last and first-hit levels.
    void walk_level(const Unitary& V, int l, double est, Result& result,
                    const Progress& progress = {});
    // Sets the lower bound and certification from the last full scale walked.
    void conclude(Result& result) const;

  private:
    double bound(int K, bool model) const;
    bool model_worth_it(const StoppingModel& model, int K, const Result& result) const;
    bool provably_optimal(int K, const Result& result, bool model) const;

    Settings settings_;
    Price price_;
    Enumerator enumerator_;
    std::shared_ptr<const StoppingModel> model_;
};

// The summary of a scan's walks, bounds and answer; the split fields keep their
// defaults.
SearchSummary scan_summary(const ScaleScan::Result& scan);

// Lowers a log entry's best cost to the other search's answer where cheaper.
void count_other(SearchLog& entry, const std::optional<Enumerator::Hit>& other);

}  // namespace quicksynth

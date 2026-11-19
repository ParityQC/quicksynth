// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/scan_then_split.hpp"

#include "quicksynth/data_structures/prefix_sets.hpp"

#include <utility>

namespace quicksynth {

namespace {

// Returns the split's first composite level (none: the split's default): one
// past the scan's last level, or with `from_hit` the scan's first hit level if
// it has one. A walked level without a hit has no element within epsilon, and
// one with a hit holds every element the split would find there.
std::optional<int> split_start(const ScaleScan::Result& s, bool from_hit) {
    if (!s.last_level) return std::nullopt;
    if (from_hit && s.first_hit_level) return s.first_hit_level;
    return *s.last_level + 1;
}

}  // namespace

int prefix_scale_for(double epsilon, const GateSet& gateset) {
    if (epsilon >= kDeepEpsilon) return 2;
    if (gateset.step == 1) return 3;
    return epsilon >= kDeeperEpsilon ? 4 : 6;
}

std::vector<int> split_scales(const ScanThenSplit::Settings& settings, int scale) {
    return settings.fallback_scales ? prefix_scales(scale) : std::vector<int>{scale};
}

double split_budget_for(double epsilon) {
    return epsilon >= kDeepEpsilon ? kSplitPointBudget : 3e8;
}

ScanThenSplit::ScanThenSplit(double epsilon, GateSet gateset, Price price, int threads,
                             std::shared_ptr<const CostBall> ball, Settings settings,
                             std::shared_ptr<const CostBall> stopping_ball)
    : price_(price),
      split_below_(settings.split_below),
      split_levels_(settings.split_levels),
      level_search_rate_(
          settings.level_search
              ? (settings.level_search_rate ? settings.level_search_rate : settings.scan.rate)
              : std::nullopt),
      skip_setup_levels_(settings.skip_setup_levels),
      scan_(epsilon, gateset, price, threads, settings.scan) {
    if (!ball) return;
    const int prefix_scale =
        settings.prefix_scale ? *settings.prefix_scale : prefix_scale_for(epsilon, gateset);
    const double split_budget =
        settings.split_budget ? *settings.split_budget : split_budget_for(epsilon);
    const PrefixSets& fixed =
        ball->prefix_sets(split_scales(settings, prefix_scale), settings.max_prefixes);
    // model stopping rule only: extra prefix scales after the fixed sets
    PrefixSets extra;
    const std::vector<int> more = stopping_extra_scales(prefix_scale);
    if (settings.scan.stopping_rule == StoppingRule::model && !more.empty())
        extra = (stopping_ball ? *stopping_ball : *ball).prefix_sets(more, settings.max_prefixes);
    split_ = std::make_unique<PrefixSplit>(
        epsilon, std::move(gateset), price, threads, fixed, extra,
        PrefixSplit::Settings{split_budget, settings.split_rate, settings.split_partial,
                              scan_.levels()});
}

ScanThenSplit::Result ScanThenSplit::run(const Unitary& V, const ScaleScan::Progress& progress) {
    if (level_search_rate_ && scan_.levels() > 1 && scan_.enumerator().gateset().step == 1) {
        LevelSearch search(scan_, split_.get(), price_,
                           {*level_search_rate_, split_below_, skip_setup_levels_});
        LevelSearch::Result r = search.run(V, progress);
        return {std::move(r.best), std::move(r.summary), std::move(r.log)};
    }
    const ScaleScan::Result s = scan_.run(V, progress);
    Result out{s.best, scan_summary(s), s.log};
    SearchSummary& summary = out.summary;
    // A certified scan leaves nothing for the split.
    const double epsilon = scan_.enumerator().epsilon();
    const StoppingModel* model = scan_.model();
    if (!split_ || s.certified || !(model || !s.best || epsilon < split_below_)) return out;
    summary.split_ran = true;
    PrefixSplit::Result r;
    if (model) {
        r = split_->run_model(V, *model, s, scan_.settings().rate.value_or(0.0),
                              epsilon < split_below_);
    } else {
        const int L = scan_.levels();
        // Finer levels start past the scan; Clifford+T starts at the first hit
        // level (its pinned replays depend on it).
        const bool from_hit = L == 1 || scan_.enumerator().gateset().step != 1;
        r = split_->run(V, split_levels_ ? *split_levels_ : L * (kFixedSplitScales - 1) + 1,
                        split_start(s, from_hit),
                        s.best ? std::optional<int64_t>(s.best->cost.value) : std::nullopt);
    }
    // log best cost includes the scan's answer (the fixed split's log omits it)
    for (SearchLog e : r.log) {
        count_other(e, s.best);
        out.log.push_back(e);
    }
    if (r.best && (!out.best || r.best->cost.value < out.best->cost.value)) {
        out.best = r.best;
        summary.source = Source::split;
    }
    return out;
}

}  // namespace quicksynth

// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/level_search.hpp"

#include "quicksynth/search_scales/scale_bound.hpp"
#include "quicksynth/search_scales/stopping_rule.hpp"
#include "quicksynth/enumerate_lattice/search_region.hpp"

#include <algorithm>
#include <memory>
#include <utility>

namespace quicksynth {

namespace {

// Predicted LLL set-up seconds per walk for one search's directions. The
// enumerator reduces sub 0 afresh, derives each finer sub from the next coarser
// and keeps every sub on the way, so subs 0..finest are built.
class Setups {
  public:
    double seconds(int sub) const {
        if (sub <= finest_) return 0.0;
        const double fresh = finest_ < 0 ? kFreshSetupSeconds : 0.0;
        return fresh + kDerivedSetupSeconds * (sub - std::max(finest_, 0));
    }
    void built(int sub) { finest_ = std::max(finest_, sub); }

  private:
    int finest_ = -1;
};

}  // namespace

LevelSearch::LevelSearch(ScaleScan& scan, PrefixSplit* split, Price price, Settings settings)
    : scan_(scan), split_(split), price_(price), settings_(settings) {}

LevelSearch::Result LevelSearch::run(const Unitary& V, const ScaleScan::Progress& progress) {
    const int det_classes = GateSet::kDeterminantClasses;
    Enumerator& enumerator = scan_.enumerator();
    const GateSet& gateset = enumerator.gateset();
    const int L = scan_.levels();
    const double scan_budget = scan_.settings().point_budget;
    const double split_budget = split_ ? split_->settings().point_budget : 0.0;
    Result out;
    ScaleScan::Result scan;
    std::optional<Enumerator::Hit> split_best;
    std::optional<PrefixSplit::Prefixed> prefixed;  // the split's prefix set
    Setups scan_setups, split_setups;
    bool split_ran = false;
    const auto walk_scan = [&](int l, double est) {
        const size_t logged = scan.log.size();
        scan_.walk_level(V, l, est, scan, progress);
        scan_setups.built(level_at(l, L).sub);
        for (size_t i = logged; i < scan.log.size(); ++i) {
            out.log.push_back(scan.log[i]);
            count_other(out.log.back(), split_best);
        }
    };
    const auto walk_split = [&](int l) {
        out.log.push_back(split_->walk_level(*prefixed, l, split_best));
        split_setups.built(out.log.back().sub);
        count_other(out.log.back(), scan.best);
        split_ran = true;
    };
    // The cheaper answer (the scan's on a tie) and the summary of what was walked.
    const auto finish = [&] {
        scan_.conclude(scan);
        out.summary = scan_summary(scan);
        out.summary.split_ran = split_ran;
        out.best = scan.best;
        if (split_best && (!out.best || split_best->cost.value < out.best->cost.value)) {
            out.best = split_best;
            out.summary.source = Source::split;
        }
        return std::move(out);
    };

    // Whether the split leaves out level y before the first answer: its set-ups
    // take longer than its walks and the next level, which holds its elements,
    // stays a factor kSkipHeadroom under the split budget.
    const auto setup_bound = [&](int y, const PrefixSplit::Level& at) {
        if (!(split_setups.seconds(at.remainder.sub) > walk_seconds(at.est))) return false;
        const std::optional<PrefixSplit::Level> next = split_->level(*prefixed, y + 1);
        return next && !(next->walks * next->est * kSkipHeadroom > split_budget);
    };

    // Before the first answer: the scan's levels in order.
    const int first = enumerator.region().first_level(L);
    for (int l = first;; ++l) {
        const auto [K, sub] = level_at(l, L);
        const double est = enumerator.estimate_points(K, sub);
        if (est > scan_budget) break;
        walk_scan(l, est);
        if (scan.best) break;
    }
    int cov = scan.last_level.value_or(first - 1);  // the scan's last level
    int scov = cov;                                 // the split's; none while <= cov
    std::optional<int> answer_level = scan.first_hit_level;
    // Without a scan answer: the split's levels past the scan's, each fixed set in turn.
    for (size_t set = 0; !scan.best && split_ && !split_best && set < split_->fixed_sets(); ++set) {
        prefixed = split_->prefixed(V, set);
        split_setups = Setups{};
        scov = cov;
        for (int y = cov + 1; !prefixed->problems.empty() && y <= cov + L * kSplitScaleGuard; ++y) {
            const std::optional<PrefixSplit::Level> at = split_->level(*prefixed, y);
            if (!at) continue;
            if (at->walks * at->est > split_budget) break;
            if (settings_.skip_setup_levels && setup_bound(y, *at)) continue;
            walk_split(y);
            scov = y;
            if (split_best) {
                answer_level = y;
                break;
            }
        }
    }
    if (!answer_level) return finish();

    // After it: the best-value candidate level while one is worth its time.
    const bool split_after = split_ && (!scan.best || enumerator.epsilon() < settings_.split_below);
    for (size_t set = 0; split_after && !prefixed && set < split_->fixed_sets(); ++set) {
        PrefixSplit::Prefixed p = split_->prefixed(V, set);
        if (!p.problems.empty()) prefixed = std::move(p);
    }
    const std::unique_ptr<LevelGain> gain = level_gain(gateset, price_, enumerator.epsilon(), V,
                                                       enumerator.first_scale(), *answer_level);
    const LevelRates rates = gain->rates(settings_.rate);
    const auto best_cost = [&] {
        if (!split_best) return scan.best->cost.value;
        if (!scan.best) return split_best->cost.value;
        return std::min(scan.best->cost.value, split_best->cost.value);
    };
    for (;;) {
        const int64_t best = best_cost();
        struct Pick {
            bool split;
            int level;
            double est;
        };
        std::optional<Pick> pick;
        double value = 0.0;
        for (int a = cov + 1;; ++a) {
            const auto [K, sub] = level_at(a, L);
            if (scan_.scale_lower_bound(K) >= static_cast<double>(best)) break;
            const double est = enumerator.estimate_points(K, sub);
            if (est > scan_budget) break;
            const double seconds = det_classes * (walk_seconds(est) + scan_setups.seconds(sub));
            const double v = gain->scan(cov, a, scov, best) - rates.scan * seconds;
            if (v > value) {
                value = v;
                pick = Pick{false, a, est};
            }
        }
        const int lo = std::max(cov, scov);
        for (int y = lo + 1; split_after && prefixed && y <= lo + L * kSplitScaleGuard; ++y) {
            const std::optional<PrefixSplit::Level> at = split_->level(*prefixed, y);
            if (!at) continue;
            if (at->walks * at->est > split_budget) break;
            const std::optional<double> bound =
                split_scale_bound(gateset, price_, level_at(y, L).K);
            if (bound && *bound >= static_cast<double>(best)) break;
            const double seconds =
                at->walks * (walk_seconds(at->est) + split_setups.seconds(at->remainder.sub));
            const double v = gain->split(lo, y, best) - rates.split * seconds;
            if (v > value) {
                value = v;
                pick = Pick{true, y, at->est};
            }
        }
        if (!pick) break;
        if (pick->split) {
            walk_split(pick->level);
            scov = pick->level;
        } else {
            walk_scan(pick->level, pick->est);
            cov = pick->level;
        }
    }
    return finish();
}

}  // namespace quicksynth

// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/scan.hpp"

#include "quicksynth/search_scales/scale_bound.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <utility>

namespace quicksynth {

namespace {

SearchLog walk_log(int scale, int sub, int det_class, double seconds, double estimated,
                   int64_t points, const std::optional<Enumerator::Hit>& best) {
    SearchLog e;
    e.scale = scale;
    e.sub = sub;
    e.det_class = det_class;
    e.seconds = seconds;
    e.estimated = estimated;
    e.points = points;
    if (best) e.best_cost = best->cost.value;
    return e;
}

}  // namespace

ScaleScan::ScaleScan(double epsilon, GateSet gateset, Price price, int threads, Settings settings)
    : settings_(settings),
      price_(price),
      enumerator_(epsilon, std::move(gateset), price, threads, Enumerator::Settings{}),
      model_(settings.stopping_rule == StoppingRule::model
                 ? stopping_model(enumerator_.gateset(), price, enumerator_.epsilon(),
                                  enumerator_.first_scale())
                 : nullptr) {
    if (!settings_.throughput) settings_.throughput = reference_throughput(enumerator_.gateset());
}

double ScaleScan::bound(int K, bool model) const {
    const GateSet& gateset = enumerator_.gateset();
    if (model) return static_cast<double>(integer_scale_bound(gateset, price_, K));
    return scan_scale_bound(gateset, price_, K);
}

bool ScaleScan::model_worth_it(const StoppingModel& model, int K, const Result& result) const {
    // Model test from state (K - 1, best). Scale time o + det_classes est /
    // theta; infinite over the point budget.
    const TimeModel& time = time_model(enumerator_.gateset());
    const double theta = time.scan_throughput;
    const double det_classes = GateSet::kDeterminantClasses;
    const auto seconds = [&](int g) {
        const double est = enumerator_.estimate_points(g);
        if (est > settings_.point_budget) return std::numeric_limits<double>::infinity();
        return time.scan_overhead + det_classes * est / theta;
    };
    // std::cref avoids a std::function allocation
    return model.walk(K - 1, result.best->cost.value, *settings_.rate, std::cref(seconds));
}

bool ScaleScan::provably_optimal(int K, const Result& result, bool model) const {
    if (!result.best) return false;
    return bound(K, model) >= static_cast<double>(result.best->cost.value);
}

int ScaleScan::levels() const {
    return settings_.fine_levels && !model_ ? levels_per_scale(enumerator_.gateset()) : 1;
}

ScaleScan::Result ScaleScan::run(const Unitary& V, const Progress& progress) {
    const int L = levels();
    const int extra = settings_.extra_levels ? *settings_.extra_levels : L * settings_.extra_scales;
    Result result;
    const StoppingModel* stopping = model_.get();
    const bool model = stopping != nullptr;
    int levels_after_hit = 0;
    for (int l = L == 1 ? enumerator_.first_scale() : enumerator_.region().first_level(L);; ++l) {
        // K: the scale of every element not yet walked
        const auto [K, sub] = level_at(l, L);
        if (provably_optimal(K, result, model)) {
            result.certified = true;
            break;
        }
        const double est = enumerator_.estimate_points(K, sub);
        if (settings_.rate && result.best &&
            !(model ? model_worth_it(*stopping, K, result)
                    : fixed_worth_it(est, result.best->cost.value, *settings_.rate,
                                     *settings_.throughput)))
            break;
        if (est > settings_.point_budget) break;
        walk_level(V, l, est, result, progress);
        if (!model && result.first_hit_scale && ++levels_after_hit > extra) break;
    }
    conclude(result);
    return result;
}

void ScaleScan::walk_level(const Unitary& V, int l, double est, Result& result,
                           const Progress& progress) {
    const int det_classes = GateSet::kDeterminantClasses;
    const auto [K, sub] = level_at(l, levels());
    bool hit_here = false;
    for (int det_class = 0; det_class < det_classes; ++det_class) {
        const Enumerator::Result r = enumerator_.enumerate(V, K, det_class, sub);
        hit_here = hit_here || !r.hits.empty();
        for (const Enumerator::Hit& h : r.hits)
            if (!result.best || h.cost.value < result.best->cost.value) result.best = h;
        result.log.push_back(walk_log(K, sub, det_class, r.seconds, est, r.points, result.best));
        if (progress) progress(r, est);
    }
    result.estimated_points += det_classes * est;
    // a level contains every element of the scales below K
    result.last_full_scale = sub == 0 ? K : K - 1;
    result.last_level = l;
    if (hit_here && !result.first_hit_scale) {
        result.first_hit_scale = K;
        result.first_hit_level = l;
    }
}

void ScaleScan::conclude(Result& result) const {
    if (!result.last_full_scale) return;
    const double bound_next = bound(*result.last_full_scale + 1, model_ != nullptr);
    if (!result.best) {
        result.lower_bound = bound_next;
    } else {
        const double best = static_cast<double>(result.best->cost.value);
        result.lower_bound = std::min(best, bound_next);
        result.certified = result.certified || bound_next >= best;
    }
}

SearchSummary scan_summary(const ScaleScan::Result& scan) {
    SearchSummary summary;
    for (const SearchLog& e : scan.log) summary.scales.push_back(e.scale);
    std::sort(summary.scales.begin(), summary.scales.end());
    summary.scales.erase(std::unique(summary.scales.begin(), summary.scales.end()),
                         summary.scales.end());
    summary.first_hit_scale = scan.first_hit_scale;
    summary.last_full_scale = scan.last_full_scale;
    summary.lower_bound = scan.lower_bound;
    summary.estimated_points = scan.estimated_points;
    if (scan.best) summary.scan_cost = scan.best->cost.value;
    summary.scan_certified = scan.certified;
    return summary;
}

void count_other(SearchLog& entry, const std::optional<Enumerator::Hit>& other) {
    if (other && (!entry.best_cost || other->cost.value < *entry.best_cost))
        entry.best_cost = other->cost.value;
}

}  // namespace quicksynth

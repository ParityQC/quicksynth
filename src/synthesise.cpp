// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/synthesise.hpp"

#include "quicksynth/synthesise_rotation/diagonal_form.hpp"
#include "quicksynth/synthesise_rotation/rotation_search.hpp"
#include "quicksynth/exact_synthesis/descent.hpp"
#include "quicksynth/exact_synthesis/normal_form.hpp"
#include "quicksynth/data_structures/distance.hpp"
#include "quicksynth/data_structures/prefix_sets.hpp"
#include "quicksynth/support/clock.hpp"
#include "quicksynth/support/host.hpp"
#include "quicksynth/support/lru.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <tuple>

namespace quicksynth {

namespace {

// Searched price; Clifford+T uses t only.
Price searched_price(const GateSet& gateset, const Price& price) {
    return gateset.step == 1 ? price : Price::clifford_t(price.t);
}

std::string format(const char* fmt, double x) {
    char buf[64];
    std::snprintf(buf, sizeof buf, fmt, x);
    return buf;
}

// Respells a hit in normal form at its cost (magic gates kept, no Q on
// Clifford+T: test_synthesise.cpp). Accepted if its distance is <= epsilon -
// kDistanceError, independent of the search margin.
std::optional<Synthesised> accept(const Unitary& U, double epsilon, const std::string& raw,
                                  const Cost& cost) {
    std::string word = normal_form(raw);
    const double distance = precise_distance(word, U);
    if (!(distance <= epsilon - kDistanceError)) return std::nullopt;
    Synthesised out;
    out.word = std::move(word);
    out.cost = cost;
    out.distance = distance;
    return out;
}

static_assert(kEpsilonMargin > kDistanceError);

// (t, q, use_q, radius or scale)
using BallKey = std::tuple<int64_t, int64_t, bool, int>;

struct Caches {
    std::mutex mutex;
    Lru<BallKey, std::shared_ptr<CostBall>> balls{8};
    Lru<BallKey, int> radius{32};  // keyed by scale, not radius
};

Caches& caches() {
    static Caches c;
    return c;
}

// prefix_ball with the lock held.
std::shared_ptr<CostBall> locked_prefix_ball(Caches& c, const Price& price, bool use_q,
                                             int radius) {
    const BallKey key{price.t, price.q, use_q, radius};
    if (auto* hit = c.balls.find(key)) return *hit;
    auto ball = std::make_shared<CostBall>(price, use_q);
    ball->build(radius);
    c.balls.put(key, ball);
    return ball;
}

}  // namespace

// --- Prefix-ball cache -----------------------------------------------------------------

std::shared_ptr<CostBall> prefix_ball(const Price& price, bool use_q, int radius) {
    Caches& c = caches();
    std::lock_guard<std::mutex> lock(c.mutex);
    return locked_prefix_ball(c, price, use_q, radius);
}

int prefix_radius(const Price& price, bool use_q, int scale) {
    Caches& c = caches();
    std::lock_guard<std::mutex> lock(c.mutex);
    const BallKey key{price.t, price.q, use_q, scale};
    if (auto* hit = c.radius.find(key)) return *hit;
    const int limit = use_q ? 16 : 2 * scale;
    int radius = kPrefixRadius;
    for (;; ++radius) {
        poll_interrupt();
        auto ball = locked_prefix_ball(c, price, use_q, radius);
        if (ball->scale_count(scale) >= kMaxPrefixes || radius >= limit) break;
    }
    c.radius.put(key, radius);
    return radius;
}

int split_radius(const Price& price, const GateSet& gateset, double epsilon, int scale) {
    if (epsilon >= kDeepEpsilon) return kPrefixRadius;
    // every Clifford+T coset of the scale: none of scale k beyond T-count 2k
    if (gateset.step != 1 && epsilon < kDeeperEpsilon) {
        if (price.t > std::numeric_limits<int>::max() / (2 * scale))
            throw std::invalid_argument("split radius: the T price is too large");
        return 2 * scale * static_cast<int>(price.t);
    }
    return prefix_radius(price, gateset.step == 1, scale);
}

// --- QuickSynth -----------------------------------------------------------------------

QuickSynth::QuickSynth(Settings settings)
    : settings_(std::move(settings)), price_(searched_price(settings_.gateset, settings_.price)) {}

ScanThenSplit::Settings QuickSynth::first_stage(double epsilon, bool fallback) const {
    const Settings& s = settings_;
    ScanThenSplit::Settings out;
    if (s.extra_scales) out.scan.extra_scales = *s.extra_scales;
    if (s.point_budget) out.scan.point_budget = *s.point_budget;
    // model stopping rule only if the price has a model at the search target
    out.scan.stopping_rule = s.stopping_rule == StoppingRule::model && !fallback &&
                                     has_stopping_model(s.gateset, price_, epsilon - kEpsilonMargin)
                                 ? StoppingRule::model
                                 : StoppingRule::fixed;
    // Split settings follow the requested epsilon, not the search target: at
    // exactly 1e-7 or kSplitBelow the target would cross their thresholds.
    out.prefix_scale = prefix_scale_for(epsilon, s.gateset);
    out.split_budget = split_budget_for(epsilon);
    out.split_below = kSplitBelow - kEpsilonMargin;
    // every coset of the prefix scale (split_radius), no fallback sets
    if (s.gateset.step != 1 && epsilon < kDeeperEpsilon) {
        out.max_prefixes = std::nullopt;
        out.fallback_scales = false;
    }
    if (s.fast) {
        out.scan.rate = kFastRate;
        out.split_rate = kFastRate;
    }
    out.level_search_rate = s.fast ? kFineFastRate : kFineDefaultRate;
    if (s.rate) {
        out.scan.rate = *s.rate;
        out.split_rate = *s.rate;
        out.level_search_rate = *s.rate;
    }
    out.split_partial = s.partial_walk;
    // Clifford+T gains nothing from finer levels at its two levels per scale
    out.scan.fine_levels = s.fine_levels.value_or(true) && s.gateset.step == 1;
    out.scan.extra_levels = s.extra_levels;
    out.split_levels = s.split_levels;
    out.level_search = s.level_search;
    out.skip_setup_levels = s.skip_setup_levels;
    return out;
}

// Prefix ball for stopping_extra_scales, sized for kMaxPrefixes of the deepest;
// null for the fixed stopping rule or no extra scales.
std::shared_ptr<CostBall> QuickSynth::stopping_ball(const ScanThenSplit::Settings& stage) const {
    if (stage.scan.stopping_rule != StoppingRule::model) return nullptr;
    const std::vector<int> more = stopping_extra_scales(stage.prefix_scale.value());
    if (more.empty()) return nullptr;
    const bool use_q = settings_.gateset.step == 1;
    return prefix_ball(price_, use_q, prefix_radius(price_, use_q, more.back()));
}

// A stage's split prefix ball at split_radius(epsilon).
std::shared_ptr<CostBall> QuickSynth::stage_ball(double epsilon, int prefix_scale) const {
    const GateSet& gateset = settings_.gateset;
    const int radius = split_radius(price_, gateset, epsilon, prefix_scale);
    return prefix_ball(price_, gateset.step == 1, radius);
}

int QuickSynth::rotation_solvable() const {
    if (settings_.rotation_solvable) return *settings_.rotation_solvable;
    return settings_.fast ? kRotationSolvableFast : kRotationSolvable;
}

// Stage 0: the word within epsilon of U that the z-rotation path found, and
// in `final` whether it is the answer. Otherwise the stages run and it is
// returned only if they find nothing cheaper: Clifford+sqrt(T) stopped by a
// cap before the count rule. Nullopt: no diagonal form within the tolerance,
// no candidate, or an answer accept() rejects. The walks go to `log`.
std::optional<Synthesised> QuickSynth::rotation_stage(const Unitary& U, double epsilon,
                                                      std::optional<double> deadline,
                                                      std::vector<std::pair<int, SearchLog>>& log,
                                                      bool& final) const {
    const Settings& s = settings_;
    const bool use_q = s.gateset.step == 1;
    // kEpsilonMargin: an exact target's residual is rounding, up to about 1e-16
    const std::optional<DiagonalForm> form =
        diagonal_form(U, s.rotation_tolerance * epsilon + kEpsilonMargin);
    if (!form) return std::nullopt;
    // d(C1 W C2, U) <= d(W, D) + residual: the chord distance is unitarily
    // invariant and obeys the triangle inequality. So W within `target` of D
    // spells a word within epsilon - kEpsilonMargin of U.
    const double target = epsilon - form->residual - kEpsilonMargin;
    if (!(target > 0)) return std::nullopt;
    RotationSearch::Settings settings;
    settings.solvable = rotation_solvable();
    settings.deadline = deadline;
    RotationSearch::Result found = RotationSearch(s.gateset, price_, settings).run(form->d, target);
    for (const SearchLog& entry : found.log) log.emplace_back(0, entry);
    // past the deadline no stage starts
    final = !use_q || found.stop == RotationStop::bound || found.stop == RotationStop::count ||
            found.stop == RotationStop::time_limit;
    std::optional<Synthesised> answer;
    if (found.best) {
        const CyclotomicMatrix& e = found.best->element;
        const std::string word =
            form->left +
            Descent(use_q ? price_.t : 1, price_.q).word(e.matrix, e.denominator_exponent) +
            form->right;
        answer = accept(U, epsilon, word, found.best->cost);
    }
    if (s.verbose)
        print_line("stage 0 rotation: cost " + (!found.best
                                                    ? std::string("None")
                                                    : std::to_string(found.best->cost.value) +
                                                          (!answer  ? ", rejected at epsilon"
                                                           : !final ? ", the stages search for less"
                                                                    : "")),
                   false);
    if (!answer) return std::nullopt;
    // The Clifford factors leave the cost unchanged.
    const Cost check = Pricer::for_gateset(s.gateset, price_).price(word_matrix(answer->word));
    if (check.value != answer->cost.value || check.n_t != answer->cost.n_t ||
        check.n_q != answer->cost.n_q)
        throw std::logic_error("stage 0: the word's price differs from the search's");
    // The search certifies the optimum around D; that is the epsilon-ball around
    // U only when D is U up to rounding. Clifford+sqrt(T) answers are never
    // certified, as in the stages.
    answer->certified = found.certified && form->residual <= kEpsilonMargin && !use_q;
    answer->stages = 0;
    answer->search = std::move(found.summary);
    return answer;
}

Synthesised QuickSynth::synthesise(const Unitary& U, double epsilon) const {
    const double start = wall_seconds();
    const Settings& s = settings_;
    Synthesised out = search(U, epsilon);
    if (s.gateset.step != 1) return out;
    if (s.time_limit && wall_seconds() > start + *s.time_limit) return out;
    Settings t_only = s;
    t_only.gateset = GateSet::clifford_t();
    t_only.verbose = false;
    std::optional<double> deadline;
    if (s.time_limit) deadline = start + *s.time_limit;
    std::optional<Synthesised> pure =
        QuickSynth(t_only).cheaper_than(U, epsilon, out.cost.value, deadline);
    if (!pure) return out;
    // respelt at the sqrt(T) price, as the stages spell (T as Q Q where 2q <= t)
    const CyclotomicMatrix e = word_matrix(pure->word);
    std::optional<Synthesised> answer =
        accept(U, epsilon, Descent(price_.t, price_.q).word(e.matrix, e.denominator_exponent), {});
    if (!answer) return out;
    answer->cost = Pricer::for_gateset(s.gateset, price_).price(word_matrix(answer->word));
    if (answer->cost.value >= out.cost.value) return out;
    if (s.verbose)
        print_line("Clifford+T comparison: cost " + std::to_string(answer->cost.value) +
                       ", cheaper than " + std::to_string(out.cost.value),
                   false);
    answer->stages = pure->stages;
    answer->search = std::move(pure->search);
    answer->search->source = Source::clifford_t;
    answer->log = std::move(out.log);
    return *answer;
}

std::optional<Synthesised> QuickSynth::cheaper_than(const Unitary& U, double epsilon, int64_t cost,
                                                    std::optional<double> deadline) const {
    const Settings& s = settings_;
    if (s.rotation_path) {
        std::vector<std::pair<int, SearchLog>> log;
        bool final = false;
        std::optional<Synthesised> found = rotation_stage(U, epsilon, deadline, log, final);
        if (found && found->cost.value < cost) return found;
        if (found) return std::nullopt;  // final on Clifford+T
    }
    const bool fallback = s.stopping_rule == StoppingRule::model && stopping_fallback(U);
    const std::vector<Stage> all = stages(epsilon, s.gateset, first_stage(epsilon, fallback));
    const double target = epsilon - kEpsilonMargin;
    const int threads = std::max(1, s.workers);
    for (size_t i = 1; i <= all.size(); ++i) {
        poll_interrupt();
        if (deadline && wall_seconds() > *deadline) return std::nullopt;
        ScanThenSplit::Settings stage = all[i - 1].settings;
        stage.scan.incumbent = cost;
        // the scan alone first: at most prices it reaches the incumbent and
        // the split's prefix ball is never built
        ScaleScan scan(target, s.gateset, price_, threads, stage.scan);
        const ScaleScan::Result alone = scan.run(U);
        // nothing below the incumbent at any scale
        if (alone.incumbent_reached && !(alone.best && alone.best->cost.value < cost))
            return std::nullopt;
        ScanThenSplit::Result found{alone.best, scan_summary(alone), {}};
        if (!alone.incumbent_reached && !alone.certified) {
            ScanThenSplit search(target, s.gateset, price_, threads,
                                 stage_ball(epsilon, stage.prefix_scale.value()), stage,
                                 stopping_ball(stage));
            found = search.run(U);
        }
        if (!found.best || found.best->cost.value >= cost) {
            if (alone.certified) return std::nullopt;  // the optimum is no cheaper
            continue;
        }
        const CyclotomicMatrix& e = found.best->element;
        std::optional<Synthesised> answer =
            accept(U, epsilon, Descent(1, price_.q).word(e.matrix, e.denominator_exponent),
                   found.best->cost);
        if (!answer) continue;
        answer->stages = static_cast<int>(i);
        answer->search = std::move(found.summary);
        return answer;
    }
    return std::nullopt;
}

Synthesised QuickSynth::search(const Unitary& U, double epsilon) const {
    const double start = wall_seconds();
    const Settings& s = settings_;
    const bool use_q = s.gateset.step == 1;
    // t_cost > 2 q_cost, overflow-free for prices >= 1
    if (use_q && price_.t - price_.q > price_.q)
        throw std::invalid_argument(
            "t_cost exceeds two sqrt(T) gates, so T is never used; set t_cost <= 2 * q_cost");
    const double target = epsilon - kEpsilonMargin;
    // targets outside the model use the fixed stopping rule in every stage
    const bool fallback = s.stopping_rule == StoppingRule::model && stopping_fallback(U);
    const std::vector<Stage> all = stages(epsilon, s.gateset, first_stage(epsilon, fallback));
    const size_t n = all.size();
    Synthesised out;
    // stage 0's word when the stages search for a cheaper one
    std::optional<Synthesised> bound;
    if (s.rotation_path) {
        std::optional<double> deadline;
        if (s.time_limit) deadline = start + *s.time_limit;
        bool final = false;
        bound = rotation_stage(U, epsilon, deadline, out.log, final);
        if (bound && final) {
            bound->log = std::move(out.log);
            return *bound;
        }
    }
    const auto stage_0 = [&] {
        bound->log = std::move(out.log);
        return *bound;
    };
    for (size_t i = 1; i <= n; ++i) {
        poll_interrupt();
        const Stage& stage = all[i - 1];
        if (s.time_limit && (i > 1 || bound) && wall_seconds() > start + *s.time_limit) {
            if (bound) return stage_0();
            throw SynthesisTimeout("no word within epsilon after " + std::to_string(i - 1) +
                                   " of " + std::to_string(n) + " stages and " +
                                   format("%.0f", wall_seconds() - start) + " s (time_limit " +
                                   s.time_limit_text + " s)");
        }
        auto ball = stage_ball(epsilon, stage.settings.prefix_scale.value());
        ScanThenSplit search(target, s.gateset, price_, std::max(1, s.workers), ball,
                             stage.settings, stopping_ball(stage.settings));
        ScaleScan::Progress progress;
        if (s.verbose)
            progress = [](const Enumerator::Result& r, double est) {
                std::vector<int64_t> costs;
                for (const Enumerator::Hit& h : r.hits) costs.push_back(h.cost.value);
                print_line(walk_line(r.K, r.det_class, r.points, std::move(costs), est, r.seconds),
                           true);
            };
        ScanThenSplit::Result found = search.run(U, progress);
        for (const SearchLog& entry : found.log) out.log.emplace_back(static_cast<int>(i), entry);
        if (s.verbose)
            print_line("stage " + std::to_string(i) + "/" + std::to_string(n) + " " +
                           stage.changes + ": cost " +
                           (found.best ? std::to_string(found.best->cost.value) : "None"),
                       false);
        if (!found.best) continue;
        const CyclotomicMatrix& e = found.best->element;
        // Clifford+T spells at t_cost 1, so no Q (T becomes Q Q only if
        // t_cost >= 2 q_cost).
        const std::string word =
            Descent(use_q ? price_.t : 1, price_.q).word(e.matrix, e.denominator_exponent);
        std::optional<Synthesised> answer = accept(U, epsilon, word, found.best->cost);
        if (!answer) continue;
        if (bound && bound->cost.value < answer->cost.value) return stage_0();
        answer->certified = found.summary.scan_certified && !use_q;
        answer->stages = static_cast<int>(i);
        answer->log = std::move(out.log);
        answer->search = std::move(found.summary);
        return *answer;
    }
    if (bound) return stage_0();
    throw std::runtime_error("no word within epsilon=" + format("%g", epsilon) + " after all " +
                             std::to_string(n) +
                             " escalation stages; this precision is beyond what the synthesiser "
                             "reaches");
}

void QuickSynth::prepare(double epsilon) const {
    const Settings& s = settings_;
    for (const Stage& stage : stages(epsilon, s.gateset, first_stage(epsilon, false))) {
        const int scale = stage.settings.prefix_scale.value();
        auto ball = stage_ball(epsilon, scale);
        ball->prefix_sets(split_scales(stage.settings, scale), stage.settings.max_prefixes);
        if (auto more = stopping_ball(stage.settings))
            more->prefix_sets(stopping_extra_scales(scale), stage.settings.max_prefixes);
    }
}

}  // namespace quicksynth

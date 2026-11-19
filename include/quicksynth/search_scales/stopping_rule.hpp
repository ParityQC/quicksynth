// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Model stopping rule: decides how far the scan and the split search past their first
// answer from a counting model instead of a fixed number of scales.
//
// Census N(K, c): gates modulo phase of least denominator scale K and cost c,
// in closed form over normal forms. For a Haar target, the gates of cell
// (K, c) within epsilon are modelled as independent Poisson with mean
// N(K, c) mu(epsilon), mu the epsilon-ball measure (ball_measure). Walking
// scale K finds all gates of least scale <= K, so group 0 holds every K <= s
// (s: first scale) and each later group a single scale.
//
// State (K, b): groups up to K walked, best cost b. A_g(c): gates of group g
// with cost <= c.
//   S_1(K, b) = sum_{c<b} 1 - exp(-mu A_{K+1}(c))   expected saving of group K+1
//   P(K, b)   = 1 - exp(-mu A_{K+1}(b - 1))         chance of a cheaper gate
// The scan walks group K+1 iff S_1(K, b) > rate sec(K+1). Where per-scale work
// grows too slowly against the overhead, a dynamic programme over the groups
// replaces this one-step test. c_int(K) is the proved stop.
//
// The split values composite scale K' with prefix scale k as
// S_split(K', k, b) - rate sec_split(K', k). The saving counts gates of scale
// K' - j first found at K', a fraction phi_k[j] (j = 0, 1, 2), over scales the
// scan did not walk. U = P R has K(P) + K(R) = K(U) + j; j = 2 only after a T
// syllable.
//
// QuickSynth uses the fixed stopping rule by default; with stopping_rule = model its
// first stage uses the model rule and later stages the fixed one.
#pragma once

#include "quicksynth/search_scales/scale_bound.hpp"
#include "quicksynth/price_candidates/price.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

namespace quicksynth {

// fixed: past the first hit, the scan walks up to extra_scales scales while
// fixed_worth_it, the split kFixedSplitScales. model: the tests above.
enum class StoppingRule { fixed, model };

// Rate rule threshold, magic states saved per second of compile time.
constexpr double kDefaultRate = 0.3;

// Fast setting: the rate rule at this threshold in the scan and the split, so no
// walk past the first answer is predicted to take over kScaleSaving x best / kFastRate
// seconds (about 1 s at cost 50).
constexpr double kFastRate = 3.0;

// Fast setting with finer levels on Clifford+sqrt(T): LevelSearch's rate (level
// gain over walk_seconds), fitted by replaying the fine-level sweep logs at cost
// parity with kFastRate on whole scales.
constexpr double kFineFastRate = 6.0;

// Finer levels without fast: LevelSearch's rate. Chosen on 40 Haar targets at 31
// precisions: mean cost 0.3 below whole scales at kDefaultRate, at 0.7x their CPU time.
constexpr double kFineDefaultRate = 0.2;

// Fixed stopping rule: predicted saving of a further scale, as a fraction of best cost.
constexpr double kScaleSaving = 0.06;

// Single-thread walk throughput: estimated points per second on
// Clifford+sqrt(T), leaves per second on Clifford+T (~1 leaf per point).
constexpr double kReferenceThroughputSqrtT = 8.3e6;
constexpr double kReferenceThroughputCliffordT = 2.1e7;

double reference_throughput(const GateSet& gateset);

// Composite scales the fixed-rule split walks from its first hit.
constexpr int kFixedSplitScales = 3;

// Rate rule of the fixed stopping rule: whether kScaleSaving x best_cost exceeds `rate` x the
// predicted time of a scale with `est` estimated points per determinant class.
bool fixed_worth_it(double est, int64_t best_cost, double rate, double throughput);

// Model-rule walk time, in seconds:
//   scan scale:      o + det_classes est / theta
//   split composite: det_classes M_k (a + est(K' - k) / theta_split)
// Least-squares fit (relative error) to quiet one-thread timings of
// tests/fixtures/calibration_targets.json, the three Clifford+sqrt(T) prices
// and Clifford+T, epsilon 1e-3..1e-8. Covers only walks made with an answer in
// hand; the first scale and per-call overhead (0.1-1.0 ms) are not modelled.
// Known error on Clifford+sqrt(T): large walks run at 7.7e6 points/s at 1e-7 vs
// 1.1e7 at 1e-3; split walks above 1e4 points run 2.2-2.9x faster than priced,
// those of 1e2-1e4 points take 1.1-1.3x their priced time.
struct TimeModel {
    double scan_throughput;   // theta, estimated points/s
    double scan_overhead;     // o, s per scale
    double split_overhead;    // a, s per prefix x class x composite
    double split_throughput;  // theta_split, estimated points/s
};
inline constexpr TimeModel kTimeModelSqrtT{9.61e6, 1.88e-3, 1.85e-4, 2.71e6};
inline constexpr TimeModel kTimeModelCliffordT{1.77e7, 1.78e-5, 5.28e-5, 1.94e7};

// Not constexpr: GateSet holds a vector.
inline const TimeModel& time_model(const GateSet& gateset) {
    return gateset.step == 1 ? kTimeModelSqrtT : kTimeModelCliffordT;
}

// Prefix scales the model-rule split adds to the fixed ones: those up to 3.
std::vector<int> stopping_extra_scales(int prefix_scale);

// Census N(K, c) for K = 0..kmax, c = 0..cost_cap(), as doubles (counts exceed
// 1e40; exact below 2^53 via Pascal's rule). A T syllable costs t_eff =
// t_syllable_cost; Clifford+T ignores q. Rows are complete: every gate of scale
// <= kmax costs at most cost_cap(). With `by_level` a row is a finer level
// l = 0..L kmax instead (SearchRegion::first_level): the gate's delta-adic exponent,
// the first level that holds it.
class Census {
  public:
    Census(const GateSet& gateset, const Price& price, int kmax, bool by_level = false);
    int kmax() const { return kmax_; }
    int rows() const { return rows_; }
    int64_t cost_cap() const { return cmax_; }
    int64_t t_eff() const { return t_; }
    double count(int row, int64_t c) const { return n_[static_cast<size_t>(row) * width() + c]; }

  private:
    size_t width() const { return static_cast<size_t>(cmax_) + 1; }
    bool sqrt_t_;
    int64_t t_, q_;
    int kmax_, rows_;
    int64_t cmax_;
    std::vector<double> n_;  // rows x (cmax + 1), row-major
};

// Returns phi_k[0..2] for the `prefixes` cheapest prefixes of scale k (sets
// are nested, so the count identifies the set). Sampled from normal forms in
// the cells of the 1e-4..1e-6 optima. Measured for first-stage sets at (1, 1),
// (1, 3), (2, 3) and multiples on Clifford+sqrt(T), and on Clifford+T; other
// sets get a default, k >= 4 uses k = 3.
std::array<double, 3> split_coverage(const GateSet& gateset, const Price& price, int k,
                                     int64_t prefixes);

// Grouped intensities and their sums for one (gate set, price, epsilon), over
// the scales a search can walk. Cached by stopping_model.
class StoppingModel {
  public:
    StoppingModel(const GateSet& gateset, const Price& price, double epsilon, int start);

    const GateSet& gateset() const { return gateset_; }
    const Price& price() const { return price_; }
    int start() const { return start_; }  // s, scale of group 0
    int groups() const { return groups_; }
    int64_t cost_cap() const { return cmax_; }  // costs >= this mean "none"
    double mu() const { return mu_; }
    int64_t c_int(int K) const { return integer_scale_bound(gateset_, price_, K); }

    // S_1(K, best); 0 past the last group.
    double saving(int K, int64_t best) const;

    // Whether to walk group K+1 from (K, best). `seconds(g)` is scale g's time,
    // infinite over the point budget. Uses the dynamic programme when
    // dp_needed (next two scales' times too close for the one-step test).
    bool walk(int K, int64_t best, double rate, const std::function<double(int)>& seconds) const;
    bool dp_needed(int K, const std::function<double(int)>& seconds) const;
    // Q(K, best): value of walking group K+1 then stopping optimally, minus
    // value of stopping now. Walk iff > 0.
    double dp_value(int K, int64_t best, double rate,
                    const std::function<double(int)>& seconds) const;

    // Expected saving of composite scale Kc with coverage `phi`, over scales
    // past `last_full`.
    double split_saving(int Kc, const std::array<double, 3>& phi, int64_t best,
                        int last_full) const;
    // Chance of a hit cheaper than `best`; best = cost_cap() + 1 gives the
    // chance of any hit.
    double split_p_cheaper(int Kc, const std::array<double, 3>& phi, int64_t best,
                           int last_full) const;

    // Lambda_g(<= c - 1) for group g, c = 0..cost_cap().
    const std::vector<double>& intensity_sums(int g) const { return cum_[static_cast<size_t>(g)]; }

  private:
    int64_t clamp(int64_t best) const;
    // Expected hits cheaper than c at composite scale Kc:
    // sum_j phi[j] Lambda_{Kc-j}(<= c - 1) over groups past last_full.
    double split_intensity(int Kc, const std::array<double, 3>& phi, size_t c, int last_full) const;
    GateSet gateset_;
    Price price_;
    int start_;
    int groups_;
    int64_t cmax_;
    double mu_;
    // cum_[g][c] = mu A_g(c - 1)
    std::vector<std::vector<double>> cum_;
    // next_saving_[g][b] = sum_{c<b} 1 - exp(-cum_[g][c + 1])
    std::vector<std::vector<double>> next_saving_;
    // DP tables, m = least cost in group g: none_below_[g][c] = P(m >= c) =
    // exp(-cum_[g][c]), least_at_[g][c] = P(m = c)
    std::vector<std::vector<double>> none_below_, least_at_;
};

// Largest model cost cap; tables hold cost_cap() + 1 doubles per scale. The cap
// is 98-134 at (1, 3) over 1e-3..1e-8. It is exceeded at 1e-8 by a sqrt(T)
// price above 372 (t = 1) or a Clifford+T price above 138 (511 and 204 at
// 1e-3); a sqrt(T) price of 1e9 would need terabytes.
constexpr int64_t kMaxModelCosts = 1 << 14;

// Whether the model's cost cap is at most kMaxModelCosts, for a scan starting at
// `start` or at epsilon's first scale. Otherwise QuickSynth uses the fixed
// stopping rule and skips the model's prefix sets.
bool has_stopping_model(const GateSet& gateset, const Price& price, int start);
bool has_stopping_model(const GateSet& gateset, const Price& price, double epsilon);

// Returns the cached model, built on first use (LRU of 32). `start` is the
// scan's first scale; Clifford+T ignores q. Null if the cost cap would exceed
// kMaxModelCosts. Thread-safe.
std::shared_ptr<const StoppingModel> stopping_model(const GateSet& gateset, const Price& price,
                                                    double epsilon, int start);

// --- Finer levels (level_search.hpp) -------------------------------------------------

// Seconds of one walk with `est` estimated points and a cached set-up, on
// Clifford+sqrt(T): log-log between decade knots measured on one thread
// (21 us at 0.1 points, 3.5 ms at 1e4, 123 ms at 1e6), flat below the first
// and slope 0.75 above the last.
double walk_seconds(double est);

// LLL set-up per walk direction, on top of walk_seconds: a fresh reduction, and
// a finer level's basis derived from the next coarser one (next_sublattice_basis).
constexpr double kFreshSetupSeconds = 2e-4;
constexpr double kDerivedSetupSeconds = 6e-5;

// Magic states per predicted second for scan and split levels.
struct LevelRates {
    double scan, split;
};

// Expected saving, in cost units, of walking further finer levels with best cost
// `best` in hand. Levels nest, so walking level a walks every level up to a.
class LevelGain {
  public:
    virtual ~LevelGain() = default;
    // Scan levels (cov, a]; the split has walked (cov, scov] at its coverage
    // (scov <= cov: none).
    virtual double scan(int cov, int a, int scov, int64_t best) const = 0;
    // Split levels (lo, y], lo the deeper of both searches' last levels.
    virtual double split(int lo, int y, int64_t best) const = 0;
    // Scan and split rates for the search's rate, fitted with the gain.
    virtual LevelRates rates(double rate) const = 0;
};

// Census gain: gates of level l and cost c within epsilon as Poisson with mean
// N(l, c) mu(epsilon) (Census by_level). A range saves
//   sum_{c<best} 1 - exp(-w mu A(c)),
// A(c) its gates of cost <= c, w = 1 for the scan and kLevelSplitCoverage for
// the split. Rates: scan at the search's rate, split at the rate over
// kCensusSplitRateDivisor. Fitted by replaying the 20-target fine-level sweep logs.
constexpr double kLevelSplitCoverage = 0.25;
constexpr double kCensusSplitRateDivisor = 4.0;
class CensusGain final : public LevelGain {
  public:
    // Levels up to L (start + model span), start the scan's first scale; requires
    // has_stopping_model(gateset, price, start).
    CensusGain(const GateSet& gateset, const Price& price, double epsilon, int start);
    double scan(int cov, int a, int scov, int64_t best) const override;
    double split(int lo, int y, int64_t best) const override;
    LevelRates rates(double rate) const override { return {rate, rate / kCensusSplitRateDivisor}; }

  private:
    // mu x gates of levels (lo, hi] and cost <= c
    double intensity(int lo, int hi, size_t c) const;
    template <class Intensity>
    double saving(int64_t best, const Intensity& at) const;
    std::shared_ptr<const std::vector<std::vector<double>>> counts_;  // [l][c], cumulative
    double mu_;
};

// Decay gain, for targets the census does not describe (stopping_fallback, or no
// model for the price): level o past the first answer's saves
// best kDecayFirstSaving kDecayPerLevel^o; the split covers kDecaySplitCoverage.
// Rates: the search's rate over kDecayScanRateDivisor and kDecaySplitRateDivisor.
constexpr double kDecayFirstSaving = 0.026;
constexpr double kDecayPerLevel = 0.8;
constexpr double kDecaySplitCoverage = 0.6;
constexpr double kDecayScanRateDivisor = 2.0;
constexpr double kDecaySplitRateDivisor = 6.0;
class DecayGain final : public LevelGain {
  public:
    explicit DecayGain(int first_answer_level) : first_(first_answer_level) {}
    double scan(int cov, int a, int scov, int64_t best) const override;
    double split(int lo, int y, int64_t best) const override;
    LevelRates rates(double rate) const override {
        return {rate / kDecayScanRateDivisor, rate / kDecaySplitRateDivisor};
    }

  private:
    // saving fraction of levels up to l
    double reached(int l) const;
    int first_;
};

// CensusGain, or DecayGain where stopping_fallback(V) or the price has no model.
std::unique_ptr<LevelGain> level_gain(const GateSet& gateset, const Price& price, double epsilon,
                                      const Unitary& V, int start, int first_answer_level);

// Whether V is outside the model: fixed by a nontrivial U -> C U C^-1 or
// U -> C U^T C^-1 (C Clifford), so its epsilon-ball holds multi-gate orbits.
// Detected as an axis n on an octahedral mirror plane (a zero component or two
// equal magnitudes). Covers diagonal V, symmetric V (n_y = 0), and their
// Clifford conjugates.
bool stopping_fallback(const Unitary& V);

}  // namespace quicksynth

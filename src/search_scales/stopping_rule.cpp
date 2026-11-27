// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/search_scales/stopping_rule.hpp"

#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/support/lru.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace quicksynth {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Groups past the first scale. Far beyond any point budget (work grows 256x per
// scale on Clifford+sqrt(T), 16x on Clifford+T), plus room for the split's
// prefix scale and +2 offset.
int model_span(const GateSet& gateset) { return gateset.step == 1 ? 16 : 28; }

// Census cost cap to scale kmax: gates of scale <= kmax have l <= 4 kmax
// (2 kmax on Clifford+T), which bounds E.
int64_t census_cost_cap(const GateSet& gateset, const Price& price, int kmax) {
    const int64_t t = t_syllable_cost(gateset, price);
    if (gateset.step != 1) return t * (2 * kmax + 4);
    return std::max<int64_t>(t * (2 * kmax + 4), (price.q * (4 * kmax + 8) + 2) / 3) + 2;
}

// The one-step test is optimal while sec(K+2) >= this x sec(K+1): the largest
// census growth A_{K+1}(<= c) / A_K(<= c), against 256x (16x) the work.
double one_step_ratio(const GateSet& gateset) { return gateset.step == 1 ? 16.2 : 4.0; }

// Number of syllable-type sequences of length n with n_T T's, first type f and
// last type l (0: T, 1: Q): C(n - 2, r_Q), read from `row` (Pascal row n - 2).
double sequences(int n, int n_t, int f, int l, const std::vector<double>& row) {
    const int n_q = n - n_t;
    if (n == 1) return (f == l) && ((f == 0) == (n_t == 1)) ? 1.0 : 0.0;
    const int r_t = n_t - (f == 0) - (l == 0);
    const int r_q = n_q - (f == 1) - (l == 1);
    return r_t >= 0 && r_q >= 0 ? row[static_cast<size_t>(r_q)] : 0.0;
}

// Delta-adic exponent a syllable of `type` (0: T, 1: Q) adds to a normal form.
int syllable_exponent(bool sqrt_t, int type) { return sqrt_t ? (type == 0 ? 2 : 3) : 1; }

// Delta-adic exponent l of a normal form a_1 b_1 ... a_n b_n kappa (n >= 1):
// E = 2 n_T + 3 n_Q (n_T on Clifford+T), minus overlaps for a_1 = 1 and
// kappa in <S, X>. One syllable with both is monomial.
int exponent(bool sqrt_t, int n, int n_t, int n_q, int f, int l, bool a1_one, bool kappa_in) {
    if (n == 1 && a1_one && kappa_in) return 0;
    if (sqrt_t) {
        const auto overlap = [](int type) { return syllable_exponent(true, type); };
        return 4 + 2 * n_t + 3 * n_q - (a1_one ? overlap(f) : 0) - (kappa_in ? overlap(l) : 0);
    }
    return 2 + n_t - a1_one - kappa_in;
}

}  // namespace

double reference_throughput(const GateSet& gateset) {
    return gateset.step == 1 ? kReferenceThroughputSqrtT : kReferenceThroughputCliffordT;
}

bool fixed_worth_it(double est, int64_t best_cost, double rate, double throughput) {
    const double det_classes = GateSet::kDeterminantClasses;
    const double best = static_cast<double>(best_cost);
    const double seconds = det_classes * est / throughput;
    return kScaleSaving * best > rate * seconds;
}

std::vector<int> stopping_extra_scales(int prefix_scale) {
    std::vector<int> out;
    for (int k = std::max(2, prefix_scale + 1); k <= 3; ++k) out.push_back(k);
    return out;
}

Census::Census(const GateSet& gateset, const Price& price, int kmax, bool by_level,
               std::optional<int> first_syllable)
    : sqrt_t_(gateset.step == 1),
      t_(t_syllable_cost(gateset, price)),
      q_(sqrt_t_ ? price.q : 0),
      kmax_(kmax) {
    if (price.t < 1 || (sqrt_t_ && price.q < 1) || kmax < 1)
        throw std::invalid_argument("census: a positive price and kmax");
    const int places = sqrt_t_ ? 4 : 2;
    rows_ = (by_level ? places * kmax : kmax) + 1;
    cmax_ = census_cost_cap(gateset, price, kmax);
    n_.assign(static_cast<size_t>(rows_) * width(), 0.0);
    const auto add = [&](int row, int64_t c, double v) {
        if (v == 0.0) return;
        if (row < rows_) n_[static_cast<size_t>(row) * width() + static_cast<size_t>(c)] += v;
    };
    // 24 Cliffords: 8 in <S, X> at scale 0, 16 at scale 1
    if (!first_syllable) {
        add(0, 0, 8);
        add(by_level ? places : 1, 0, 16);
    }
    const int types = sqrt_t_ ? 2 : 1;
    // Normal forms of cost <= cmax by syllable count n. `row` is Pascal row n - 2
    // up to n_Q = q_max, exact below 2^53.
    const int q_max = sqrt_t_ ? static_cast<int>(cmax_ / q_) : 0;
    const int n_max = static_cast<int>(cmax_ / (sqrt_t_ ? std::min(t_, q_) : t_));
    std::vector<double> row(static_cast<size_t>(q_max) + 1, 0.0);
    row[0] = 1.0;
    for (int n = 1; n <= n_max; ++n) {
        if (n >= 3)
            for (int r = std::min(n - 2, q_max); r > 0; --r)
                row[static_cast<size_t>(r)] += row[static_cast<size_t>(r) - 1];
        for (int n_q = 0; n_q <= std::min(n, q_max); ++n_q) {
            const int n_t = n - n_q;
            const int64_t c = t_ * n_t + q_ * n_q;
            if (c > cmax_) continue;
            // 2^(n-1) later connectors; Q and XQ per Q syllable
            const int doublings = n - 1 + n_q;
            for (int f = 0; f < types; ++f)
                for (int l = 0; l < types; ++l) {
                    if (first_syllable && f != *first_syllable) continue;
                    const double s = sequences(n, n_t, f, l, row);
                    if (s == 0.0) continue;
                    // a_1 = 1 (one way) or SH, HS^-1 (two); kappa in <S, X> (8) or not (16)
                    for (const bool a1_one : {true, false})
                        for (const bool kappa_in : {true, false}) {
                            if (first_syllable && a1_one) continue;
                            const int e = exponent(sqrt_t_, n, n_t, n_q, f, l, a1_one, kappa_in);
                            const int row = by_level ? e : (e + places - 1) / places;
                            const double ways = s * (a1_one ? 1 : 2) * (kappa_in ? 8 : 16);
                            add(row, c, std::ldexp(ways, doublings));
                        }
                }
        }
    }
}

std::array<double, 3> split_coverage(const GateSet& gateset, const Price& price, int k,
                                     int64_t prefixes) {
    if (gateset.step == 1) {
        // Keyed by (t, q, k, prefixes) with (t, q) reduced; multiples order
        // prefixes alike. Study sets: radius 6 (k = 2), 7 (k = 3). First-stage
        // sets from larger balls: k = 3 radius 10 at (2, 3); k = 2 below 1e-7
        // radius 7 at (1, 3), 10 at (2, 3).
        using Set = std::tuple<int64_t, int64_t, int, int64_t>;
        // clang-format off
        static const std::map<Set, std::array<double, 3>> kMeasured{
            {{1, 1, 2, 240}, {0.4838, 0.5162, 0.0}},
            {{1, 1, 3, 600}, {0.1721, 0.5973, 0.0127}},
            {{1, 3, 2, 128}, {0.4349, 0.5310, 0.0}},
            {{1, 3, 2, 176}, {0.4368, 0.5613, 0.0}},
            {{1, 3, 3, 600}, {0.2637, 0.6184, 0.0126}},
            {{2, 3, 2, 48}, {0.4813, 0.4475, 0.0}},
            {{2, 3, 2, 240}, {0.6446, 0.3554, 0.0}},
            {{2, 3, 3, 32}, {0.0, 0.1969, 0.0179}},
            {{2, 3, 3, 600}, {0.2734, 0.6227, 0.0065}},
        };
        // clang-format on
        const int64_t g = std::gcd(price.t, price.q);
        const auto it = kMeasured.find(Set{price.t / g, price.q / g, k, prefixes});
        if (it != kMeasured.end()) return it->second;
    } else if (k >= 2 && k <= 4 && prefixes >= (int64_t{24} << (2 * (k - 2)))) {
        // Clifford+T deep optima (cost 2K - 3) split after 2k - 1 syllables at
        // offset 0. The 24, 96, 384 cheapest prefixes of scales 2, 3, 4
        // (radius 6, 6, 8) contain all such prefixes.
        return {1.0, 0.0, 0.0};
    }
    static const std::array<std::array<double, 3>, 2> kFallback{
        {{0.45, 0.5, 0.0}, {0.2, 0.6, 0.01}}};
    return kFallback[k <= 2 ? 0 : 1];
}

StoppingModel::StoppingModel(const GateSet& gateset, const Price& price, double epsilon, int start)
    : gateset_(gateset), price_(price), start_(start) {
    const Census census(gateset, price, start + model_span(gateset));
    groups_ = census.kmax() - start + 1;
    cmax_ = census.cost_cap();
    mu_ = ball_measure(epsilon);
    const size_t costs = static_cast<size_t>(cmax_) + 1;
    cum_.assign(static_cast<size_t>(groups_), std::vector<double>(costs + 1, 0.0));
    for (int g = 0; g < groups_; ++g) {
        std::vector<double>& cum = cum_[static_cast<size_t>(g)];
        for (size_t c = 0; c < costs; ++c) {
            double n = 0.0;
            if (g == 0)
                for (int K = 0; K <= start; ++K) n += census.count(K, static_cast<int64_t>(c));
            else
                n = census.count(start + g, static_cast<int64_t>(c));
            cum[c + 1] = cum[c] + n * mu_;
        }
    }
    next_saving_.assign(cum_.size(), std::vector<double>(costs + 1, 0.0));
    for (size_t g = 0; g < cum_.size(); ++g)
        for (size_t b = 0; b < costs; ++b)
            next_saving_[g][b + 1] = next_saving_[g][b] - std::expm1(-cum_[g][b + 1]);
    none_below_.assign(static_cast<size_t>(groups_), std::vector<double>(costs + 1));
    least_at_.assign(static_cast<size_t>(groups_), std::vector<double>(costs, 0.0));
    for (size_t g = 0; g < static_cast<size_t>(groups_); ++g) {
        const std::vector<double>& cum = cum_[g];
        for (size_t c = 0; c <= costs; ++c) none_below_[g][c] = std::exp(-cum[c]);
        for (size_t c = 0; c < costs; ++c)
            least_at_[g][c] = std::exp(-cum[c]) * -std::expm1(-(cum[c + 1] - cum[c]));
    }
}

int64_t StoppingModel::clamp(int64_t best) const { return std::clamp<int64_t>(best, 0, cmax_ + 1); }

double StoppingModel::saving(int K, int64_t best) const {
    const int j = std::max(K - start_ + 1, 0);
    if (j >= groups_) return 0.0;
    return next_saving_[static_cast<size_t>(j)][static_cast<size_t>(clamp(best))];
}

bool StoppingModel::dp_needed(int K, const std::function<double(int)>& seconds) const {
    if (K + 2 - start_ >= groups_) return false;
    return seconds(K + 2) < one_step_ratio(gateset_) * seconds(K + 1);
}

bool StoppingModel::walk(int K, int64_t best, double rate,
                         const std::function<double(int)>& seconds) const {
    if (dp_needed(K, seconds)) return dp_value(K, best, rate, seconds) > 0.0;
    return saving(K, best) > rate * seconds(K + 1);
}

double StoppingModel::dp_value(int K, int64_t best, double rate,
                               const std::function<double(int)>& seconds) const {
    // Backward induction, V(i, b) = max(0, Q(i, b)), V = 0 past the last group:
    //   Q(i, b) = -rate sec(i+1) + sum_{c<b} P(m = c) (b - c + V(i+1, c))
    //             + P(m >= b) V(i+1, b)
    // m: least cost in group i+1. Q = -inf where the next scale is over budget.
    const int first = std::max(K - start_, -1);
    // Start below the first group over budget, where V = 0.
    int over = first;
    while (over <= groups_ - 2 && std::isfinite(seconds(start_ + over + 1))) ++over;
    if (over == first) return -kInf;
    const size_t costs = static_cast<size_t>(cmax_) + 1;
    const size_t b = static_cast<size_t>(clamp(best));
    std::vector<double> after(costs + 1, 0.0), q(costs + 1);
    for (int i = over - 1;; --i) {
        const std::vector<double>& none_below = none_below_[static_cast<size_t>(i + 1)];
        const std::vector<double>& least_at = least_at_[static_cast<size_t>(i + 1)];
        const double sec = seconds(start_ + i + 1);
        double a1 = 0.0, a2 = 0.0, a3 = 0.0;
        for (size_t v = 0; v <= costs; ++v) {
            q[v] = -rate * sec + static_cast<double>(v) * a1 - a2 + a3 + none_below[v] * after[v];
            if (v == costs) break;
            const double pm = least_at[v];
            a1 += pm;
            a2 += pm * static_cast<double>(v);
            a3 += pm * after[v];
        }
        if (i == first) return q[b];
        for (size_t v = 0; v <= costs; ++v) after[v] = std::max(0.0, q[v]);
    }
}

double StoppingModel::split_intensity(int Kc, const std::array<double, 3>& phi, size_t c,
                                      int last_full) const {
    double lambda = 0.0;
    for (int j = 0; j < 3; ++j) {
        const int g = Kc - j - start_;
        if (Kc - j > last_full && g >= 0 && g < groups_)
            lambda += phi[static_cast<size_t>(j)] * cum_[static_cast<size_t>(g)][c];
    }
    return lambda;
}

double StoppingModel::split_saving(int Kc, const std::array<double, 3>& phi, int64_t best,
                                   int last_full) const {
    const size_t b = static_cast<size_t>(clamp(best));
    double total = 0.0;
    for (size_t c = 0; c < b; ++c) total -= std::expm1(-split_intensity(Kc, phi, c + 1, last_full));
    return total;
}

double StoppingModel::split_p_cheaper(int Kc, const std::array<double, 3>& phi, int64_t best,
                                      int last_full) const {
    return -std::expm1(-split_intensity(Kc, phi, static_cast<size_t>(clamp(best)), last_full));
}

bool has_stopping_model(const GateSet& gateset, const Price& price, int start) {
    return census_cost_cap(gateset, price, start + model_span(gateset)) <= kMaxModelCosts;
}

bool has_stopping_model(const GateSet& gateset, const Price& price, double epsilon) {
    return has_stopping_model(gateset, price, SearchRegion(epsilon, gateset).first_scale());
}

std::shared_ptr<const StoppingModel> stopping_model(const GateSet& gateset, const Price& price,
                                                    double epsilon, int start) {
    Price priced = price;
    // q is unused on Clifford+T
    if (gateset.step != 1) priced.q = 1;
    using Key = std::tuple<int, int64_t, int64_t, double, int>;
    static std::mutex mutex;
    static Lru<Key, std::shared_ptr<const StoppingModel>> models(32);
    if (!has_stopping_model(gateset, price, start)) return nullptr;
    const Key key{gateset.step, priced.t, priced.q, epsilon, start};
    std::lock_guard<std::mutex> lock(mutex);
    if (auto* hit = models.find(key)) return *hit;
    auto model = std::make_shared<const StoppingModel>(gateset, priced, epsilon, start);
    models.put(key, model);
    return model;
}

double walk_seconds(double est) {
    // (estimated points, seconds) at decade knots
    static constexpr std::array<std::array<double, 2>, 8> kKnots{{{0.1, 21e-6},
                                                                  {1, 30e-6},
                                                                  {10, 56e-6},
                                                                  {1e2, 163e-6},
                                                                  {1e3, 685e-6},
                                                                  {1e4, 3.5e-3},
                                                                  {1e5, 18.5e-3},
                                                                  {1e6, 123e-3}}};
    constexpr double kSlopeAbove = 0.75;
    if (!(est > kKnots[0][0])) return kKnots[0][1];
    for (size_t i = 1; i < kKnots.size(); ++i) {
        const auto& [x0, y0] = kKnots[i - 1];
        const auto& [x1, y1] = kKnots[i];
        if (est <= x1) return y0 * std::pow(y1 / y0, std::log(est / x0) / std::log(x1 / x0));
    }
    const auto& [x, y] = kKnots.back();
    return y * std::pow(est / x, kSlopeAbove);
}

namespace {

using LevelCounts = std::vector<std::vector<double>>;

// Census by_level to scale `kmax`, cumulative over levels and costs: [l][c] counts
// the gates of level <= l and cost <= c (with `first_syllable`, the remainders that
// start with a connector and that syllable type). Cached (LRU of 24); thread-safe.
std::shared_ptr<const LevelCounts> level_counts(const GateSet& gateset, const Price& price,
                                                int kmax,
                                                std::optional<int> first_syllable = std::nullopt) {
    using Key = std::tuple<int, int64_t, int64_t, int, int>;
    static std::mutex mutex;
    static Lru<Key, std::shared_ptr<const LevelCounts>> cache(24);
    const Key key{gateset.step, price.t, gateset.step == 1 ? price.q : 1, kmax,
                  first_syllable.value_or(-1)};
    std::lock_guard<std::mutex> lock(mutex);
    if (auto* hit = cache.find(key)) return *hit;
    const Census census(gateset, price, kmax, true, first_syllable);
    const size_t costs = static_cast<size_t>(census.cost_cap()) + 1;
    auto counts = std::make_shared<LevelCounts>(static_cast<size_t>(census.rows()),
                                                std::vector<double>(costs, 0.0));
    for (int l = 0; l < census.rows(); ++l) {
        std::vector<double>& row = (*counts)[static_cast<size_t>(l)];
        double below = 0.0;  // this level, cost <= c
        for (size_t c = 0; c < costs; ++c) {
            below += census.count(l, static_cast<int64_t>(c));
            row[c] = below + (l > 0 ? (*counts)[static_cast<size_t>(l) - 1][c] : 0.0);
        }
    }
    cache.put(key, counts);
    return counts;
}

// Gates of levels (lo, hi] and cost <= c in cumulative `counts`.
double count_range(const LevelCounts& counts, int lo, int hi, size_t c) {
    hi = std::min(hi, static_cast<int>(counts.size()) - 1);
    if (lo >= hi) return 0.0;
    const double below = lo >= 0 ? counts[static_cast<size_t>(lo)][c] : 0.0;
    return counts[static_cast<size_t>(hi)][c] - below;
}

}  // namespace

CensusGain::CensusGain(const GateSet& gateset, const Price& price, double epsilon, int start,
                       const PrefixSet* split)
    : counts_(level_counts(gateset, price, start + model_span(gateset))),
      mu_(ball_measure(epsilon)) {
    if (!split) return;
    const int kmax = start + model_span(gateset);
    for (size_t type = 0; type < 2; ++type)
        remainders_[type] = level_counts(gateset, price, kmax, static_cast<int>(type));
    const bool sqrt_t = gateset.step == 1;
    offset_ = (sqrt_t ? 4 : 2) * split->scale;
    std::map<std::tuple<int64_t, int, bool, bool>, int> shapes;
    for (const PrefixForm& f : split->forms) {
        const auto weight = [&](int type) { return syllable_exponent(sqrt_t, type); };
        const int shift =
            f.n_t * weight(0) + f.n_q * weight(1) - (f.bare_first ? weight(*f.bare_first) : 0);
        const std::array<bool, 2> first{!f.extended[0], !f.extended[1]};
        const auto key = std::make_tuple(f.cost, shift, first[0], first[1]);
        const int shape = shapes.emplace(key, static_cast<int>(shapes.size())).first->second;
        prefixes_.push_back({f.cost, shift, first, shape});
    }
}

double CensusGain::intensity(int lo, int hi, size_t c) const {
    return mu_ * count_range(*counts_, lo, hi, c);
}

double CensusGain::reached(const Prefix& p, int lo, int hi, size_t c) const {
    if (static_cast<int64_t>(c) < p.cost) return 0.0;
    const size_t rest = c - static_cast<size_t>(p.cost);
    double n = 0.0;
    for (size_t type = 0; type < 2; ++type)
        if (p.first[type]) n += count_range(*remainders_[type], lo, hi, rest);
    return mu_ * n;
}

template <class Intensity>
double CensusGain::saving(int64_t best, const Intensity& at) const {
    // rows are complete: past the cost cap every gate of the range is counted
    const size_t last = counts_->front().size() - 1;
    double total = 0.0;
    for (int64_t c = 0; c < best; ++c)
        total -= std::expm1(-at(std::min(static_cast<size_t>(c), last)));
    return total;
}

double CensusGain::scan(int cov, int a, const std::vector<int>& depths, int64_t best) const {
    // prefixes of one shape at one depth reach the same gates: (shape, depth) ->
    // (a prefix, their number)
    std::map<std::pair<int, int>, std::pair<const Prefix*, int>> alike;
    for (size_t i = 0; i < depths.size(); ++i) {
        auto& [prefix, count] = alike[{prefixes_[i].shape, depths[i]}];
        prefix = &prefixes_[i];
        ++count;
    }
    return saving(best, [&](size_t c) {
        double v = intensity(cov, a, c);
        // less the gates P R the prefixes reached: level(R) <= depth - L k
        for (const auto& [key, group] : alike) {
            const Prefix& p = *group.first;
            v -= group.second *
                 reached(p, cov - p.shift, std::min(a - p.shift, key.second - offset_), c);
        }
        return std::max(v, 0.0);
    });
}

double CensusGain::split(int cov, const std::vector<PrefixWalk>& walks, int64_t best) const {
    return saving(best, [&](size_t c) {
        double v = 0.0;
        // remainders the walk newly reaches whose U lies past the scan's levels
        for (const PrefixWalk& w : walks) {
            const Prefix& p = prefixes_[w.prefix];
            v += reached(p, std::max(w.from - offset_, cov - p.shift), w.to - offset_, c);
        }
        return v;
    });
}

double DecayGain::reached(int l) const {
    const double o = std::max(0, l - first_);
    return kDecayFirstSaving * (1.0 - std::pow(kDecayPerLevel, o)) / (1.0 - kDecayPerLevel);
}

double DecayGain::scan(int cov, int a, const std::vector<int>& depths, int64_t best) const {
    double gain = reached(a) - reached(cov);
    // the split walks every prefix together, so their depths are equal
    const int scov = depths.empty() ? cov : depths.front();
    if (scov > cov) gain -= kDecaySplitCoverage * (reached(std::min(a, scov)) - reached(cov));
    return static_cast<double>(best) * gain;
}

double DecayGain::split(int cov, const std::vector<PrefixWalk>& walks, int64_t best) const {
    if (walks.empty()) return 0.0;
    // the walks share one level; each prefix saves its share of it
    const PrefixWalk& w = walks.front();
    const double share = static_cast<double>(walks.size()) / static_cast<double>(prefixes_);
    return static_cast<double>(best) * kDecaySplitCoverage *
           (reached(w.to) - reached(std::max(cov, w.from))) * share;
}

std::unique_ptr<LevelGain> level_gain(const GateSet& gateset, const Price& price, double epsilon,
                                      const Unitary& V, int start, int first_answer_level,
                                      const PrefixSet* split) {
    if (stopping_fallback(V) || !has_stopping_model(gateset, price, start))
        return std::make_unique<DecayGain>(first_answer_level, split ? split->prefixes.size() : 0);
    return std::make_unique<CensusGain>(gateset, price, epsilon, start, split);
}

bool stopping_fallback(const Unitary& V) {
    // V = cos(theta/2) - i s n.sigma up to phase: Im a = -s n_z, Re b = s n_y,
    // Im b = -s n_x. Mirror planes: n_i = 0, |n_i| = |n_j|; s = 0 lies on all.
    constexpr double kTolerance = 1e-12;
    const Vec4 u = su2_vector(V);
    const double x = std::abs(u[3]), y = std::abs(u[2]), z = std::abs(u[1]);
    const auto same = [](double a, double b) { return std::abs(a - b) <= kTolerance; };
    return std::min({x, y, z}) <= kTolerance || same(x, y) || same(y, z) || same(x, z);
}

}  // namespace quicksynth

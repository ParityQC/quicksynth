// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/cost_ball.hpp"

#include "quicksynth/data_structures/clifford.hpp"
#include "quicksynth/data_structures/prefix_sets.hpp"
#include "quicksynth/support/host.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace quicksynth {

namespace {

bool is_zero(const Z16& z) {
    for (int i = 0; i < kZetaDegree; ++i)
        if (z[i] != 0) return false;
    return true;
}

}  // namespace

std::vector<SyllableGenerator> syllable_generators(const Price& price, bool use_q) {
    // levels are indexed by int cost
    if (price.t > std::numeric_limits<int>::max() || price.q > std::numeric_limits<int>::max())
        throw std::invalid_argument("cost ball: a price must be below 2^31");
    std::vector<SyllableGenerator> out;
    const auto add = [&](const std::string& word, int64_t weight) {
        const CyclotomicMatrix e = word_matrix(word);
        out.push_back({word, static_cast<int>(weight), e.matrix, e.denominator_exponent});
    };
    for (const char* connector : kConnectors) {
        const std::string c = connector;
        add(c + "T", price.t);
        if (use_q)
            for (const char* sign : {"", "X"}) add(c + sign + "Q", price.q);
    }
    return out;
}

CostBall::CostBall(const Price& price, bool use_q) : price_(price), use_q_(use_q) {
    if (price.t < 1 || price.q < 1) throw std::invalid_argument("gate prices must be positive");
    // price.t > 2 price.q without overflow
    if (use_q && price.t - price.q > price.q)
        throw std::invalid_argument("a price needs t_cost <= 2 * q_cost");
    generators_ = syllable_generators(price, use_q);
    // a * b with b monomial is monomial iff connector a is
    for (const auto& g : generators_) {
        const M2& m = g.matrix;
        const bool diagonal = is_zero(m[0][1]) && is_zero(m[1][0]);
        const bool antidiagonal = is_zero(m[0][0]) && is_zero(m[1][1]);
        connector_in_K_.push_back(diagonal || antidiagonal);
    }
    Level root;
    root.cost = 0;
    root.parent_index.push_back(-1);
    root.generator.push_back(-1);
    root.sqrt2_exponent.push_back(0);
    levels_.push_back(std::move(root));
    radius_ = 0;
}

int CostBall::max_weight() const {
    int w = 0;
    for (const auto& g : generators_) w = std::max(w, g.weight);
    return w;
}

size_t CostBall::size() const {
    size_t n = 0;
    for (const auto& l : levels_) n += l.size();
    return n;
}

std::vector<int8_t> CostBall::chain(int cost, int index) const {
    std::vector<int8_t> chain;
    int c = cost, i = index;
    while (true) {
        const Level& l = levels_.at(c);
        const int g = l.generator[i];
        if (g < 0) break;
        chain.push_back(static_cast<int8_t>(g));
        i = l.parent_index[i];
        c -= generators_[g].weight;
    }
    return chain;
}

M2 CostBall::numerator(int cost, int index, int& k) const {
    const std::vector<int8_t> path = chain(cost, index);
    M2 m{};
    m[0][0][0] = 1;
    m[1][1][0] = 1;
    k = 0;
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        const SyllableGenerator& g = generators_[*it];
        m = matrix_mul(m, g.matrix);
        k += g.sqrt2_exponent;
        reduce_sqrt2(m, k);
    }
    return m;
}

const PrefixSets& CostBall::prefix_sets(const std::vector<int>& scales,
                                        std::optional<int64_t> cap) const {
    const auto key = std::make_tuple(radius_, scales, cap);
    {
        std::lock_guard<std::mutex> lock(prefix_mutex_);
        auto it = prefix_cache_.find(key);
        if (it != prefix_cache_.end()) return it->second;
    }
    // Collect unlocked: the interrupt hook takes the GIL, and a GIL holder may
    // be waiting for this lock.
    PrefixSets sets;
    for (PrefixSet& s : collect_prefixes(*this, scales, cap))
        sets.push_back(std::make_shared<const PrefixSet>(std::move(s)));
    std::lock_guard<std::mutex> lock(prefix_mutex_);
    // keeps a concurrent insert; map references are stable
    return prefix_cache_.emplace(key, std::move(sets)).first->second;
}

int64_t CostBall::scale_count(int k) const {
    int64_t n = 0;
    for (const Level& l : levels_)
        n += std::count(l.sqrt2_exponent.begin(), l.sqrt2_exponent.end(), k);
    return n;
}

std::string CostBall::word(int cost, int index) const {
    const std::vector<int8_t> path = chain(cost, index);
    std::string w;
    for (auto it = path.rbegin(); it != path.rend(); ++it) w += generators_[*it].word;
    return w;
}

void CostBall::grow(int c) {
    // The identity takes every generator; other cosets only connectors outside
    // K (normal-form successors). Each parent numerator is computed once.
    Level level;
    level.cost = c;
    std::vector<int> weights;
    for (const auto& g : generators_)
        if (std::find(weights.begin(), weights.end(), g.weight) == weights.end())
            weights.push_back(g.weight);
    std::sort(weights.begin(), weights.end());
    for (const int w : weights) {
        const int pc = c - w;
        if (pc < 0) continue;
        const Level& parent = levels_[pc];
        for (size_t pi = 0; pi < parent.size(); ++pi) {
            const bool root = parent.generator[pi] < 0;
            int k0;
            const M2 base = numerator(pc, static_cast<int>(pi), k0);
            for (size_t gi = 0; gi < generators_.size(); ++gi) {
                const SyllableGenerator& g = generators_[gi];
                if (g.weight != w || (!root && connector_in_K_[gi])) continue;
                M2 m = matrix_mul(base, g.matrix);
                int k = k0 + g.sqrt2_exponent;
                reduce_sqrt2(m, k);
                level.parent_index.push_back(static_cast<int32_t>(pi));
                level.generator.push_back(static_cast<int8_t>(gi));
                level.sqrt2_exponent.push_back(static_cast<int16_t>(k));
            }
        }
    }
    levels_.push_back(std::move(level));
}

void CostBall::build(int radius) {
    for (int c = radius_ + 1; c <= radius; ++c) {
        poll_interrupt();
        grow(c);
        radius_ = c;
    }
}

}  // namespace quicksynth

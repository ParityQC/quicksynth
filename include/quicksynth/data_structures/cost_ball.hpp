// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The weighted syllable graph, grown level by level.
//
// A node is a right coset U * Clifford (up to phase) over Z[zeta16, 1/sqrt2].
// Children are U * G for syllables G = a * b: a connector Clifford a, then a
// T-type (weight t_cost) or sqrt(T)-type (weight q_cost) syllable b. Level c
// holds every coset of cost exactly c.
//
// The graph grows along the normal-form tree: each non-identity coset has one
// parent, and children use only connectors outside K = <S, X> (non-monomial
// matrices). Children are never duplicates, so no keys are needed.
//
// Memory per coset: 4 (parent) + 1 (generator) + 2 (sqrt2 exponent) bytes.
// Numerators are recomputed from the parent chain.
#pragma once

#include "quicksynth/data_structures/cost.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace quicksynth {

struct SyllableGenerator {
    std::string word;
    int weight;
    M2 matrix{};
    int sqrt2_exponent;
};

struct Level {
    int cost = 0;
    std::vector<int32_t> parent_index;
    std::vector<int8_t> generator;
    std::vector<int16_t> sqrt2_exponent;
    size_t size() const { return parent_index.size(); }
};

struct PrefixSet;  // prefix_sets.hpp
using PrefixSets = std::vector<std::shared_ptr<const PrefixSet>>;

// Returns the syllables C T (weight price.t) and C Q, C X Q (weight price.q),
// C over kConnectors in order; Level::generator indexes this list.
// use_q = false omits the Q syllables.
std::vector<SyllableGenerator> syllable_generators(const Price& price, bool use_q);

class CostBall {
  public:
    // Throws std::invalid_argument if use_q and price.t > 2 price.q (the normal
    // form needs t <= 2q).
    CostBall(const Price& price, bool use_q = true);

    // Cached collect_prefixes (prefix_sets.hpp) for the current radius.
    // Thread-safe.
    const PrefixSets& prefix_sets(const std::vector<int>& scales, std::optional<int64_t> cap) const;
    // Grows the ball to `radius`, polling the interrupt check before each level.
    void build(int radius);
    int radius() const { return radius_; }
    int max_weight() const;
    const Level& level(int c) const { return levels_.at(c); }
    size_t size() const;
    std::string word(int cost, int index) const;
    M2 numerator(int cost, int index, int& k) const;
    // Number of cosets of denominator scale k, over all levels.
    int64_t scale_count(int k) const;

  private:
    // Normal-form generators of a coset, last syllable first.
    std::vector<int8_t> chain(int cost, int index) const;
    void grow(int c);
    Price price_;
    bool use_q_;
    std::vector<SyllableGenerator> generators_;
    std::vector<bool> connector_in_K_;  // connector is in K (monomial)
    std::vector<Level> levels_;
    int radius_ = -1;
    mutable std::mutex prefix_mutex_;
    mutable std::map<std::tuple<int, std::vector<int>, std::optional<int64_t>>, PrefixSets>
        prefix_cache_;
};

}  // namespace quicksynth

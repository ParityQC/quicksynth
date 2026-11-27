// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/prefix_sets.hpp"

#include "quicksynth/support/host.hpp"

#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

using BallPosition = std::pair<int, int>;  // (cost, index) in the cost ball

PrefixSet prefix_set(const CostBall& ball, int scale, std::vector<CyclotomicMatrix> prefixes,
                     const std::vector<BallPosition>& positions) {
    PrefixSet set;
    set.scale = scale;
    set.prefixes = std::move(prefixes);
    for (const CyclotomicMatrix& p : set.prefixes)
        set.inverses.push_back(dagger(matrix_unitary(p)));
    std::map<BallPosition, size_t> index;
    for (size_t i = 0; i < positions.size(); ++i) index[positions[i]] = i;
    set.forms.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
        const auto [cost, at] = positions[i];
        const std::vector<CostBall::Syllable> syllables = ball.syllables(cost, at);
        PrefixForm& form = set.forms[i];
        form.cost = cost;
        for (const CostBall::Syllable& s : syllables) ++(s.type == 0 ? form.n_t : form.n_q);
        if (!syllables.empty() && syllables.front().bare) form.bare_first = syllables.front().type;
        const std::optional<BallPosition> parent = ball.parent(cost, at);
        if (parent && index.count(*parent))
            set.forms[index.at(*parent)].extended[static_cast<size_t>(syllables.back().type)] =
                true;
    }
    return set;
}

}  // namespace

std::vector<int> prefix_scales(int scale) {
    std::vector<int> scales{scale};
    for (int k = scale - 1; k > 0; --k) scales.push_back(k);
    return scales;
}

std::vector<PrefixSet> collect_prefixes(const CostBall& ball, const std::vector<int>& scales,
                                        std::optional<int64_t> cap) {
    const std::set<int> wanted(scales.begin(), scales.end());
    std::vector<std::vector<CyclotomicMatrix>> found(wanted.size());
    std::vector<std::vector<BallPosition>> positions(wanted.size());
    const auto slot = [&](int scale) {
        return static_cast<size_t>(std::distance(wanted.begin(), wanted.find(scale)));
    };
    for (int c = 0; c <= ball.radius(); ++c) {
        poll_interrupt();
        const Level& level = ball.level(c);
        for (size_t i = 0; i < level.size(); ++i) {
            const int scale = level.sqrt2_exponent[i];
            if (wanted.count(scale) == 0) continue;
            CyclotomicMatrix e;
            e.matrix = ball.numerator(c, static_cast<int>(i), e.denominator_exponent);
            if (e.denominator_exponent != scale)
                throw std::logic_error("prefixes: a coset's scale is stale");
            found[slot(scale)].push_back(e);
            positions[slot(scale)].emplace_back(c, static_cast<int>(i));
        }
        if (cap && std::all_of(found.begin(), found.end(), [&](const auto& f) {
                return static_cast<int64_t>(f.size()) >= *cap;
            }))
            break;
    }
    std::vector<PrefixSet> out;
    for (int scale : scales) {
        std::vector<CyclotomicMatrix> prefixes = found[slot(scale)];
        std::vector<BallPosition> at = positions[slot(scale)];
        if (cap && static_cast<int64_t>(prefixes.size()) > *cap) {
            prefixes.resize(static_cast<size_t>(*cap));
            at.resize(static_cast<size_t>(*cap));
        }
        out.push_back(prefix_set(ball, scale, std::move(prefixes), at));
    }
    return out;
}

}  // namespace quicksynth

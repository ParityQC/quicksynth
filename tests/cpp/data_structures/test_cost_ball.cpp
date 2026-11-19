// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Cost ball: levels, words and syllable graph.

#include "fixtures.hpp"
#include "oracles.hpp"

#include "quicksynth/data_structures/clifford.hpp"
#include "quicksynth/data_structures/cost_ball.hpp"
#include "quicksynth/data_structures/cyclotomic_matrix.hpp"
#include "quicksynth/price_candidates/price.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth::test {

namespace {

// Level sizes at the catalysis price (paper).
const std::vector<size_t> kExactLevels{1, 3, 6, 18, 48, 120, 312, 816, 2112, 5472};

std::unique_ptr<CostBall> built(int64_t t, int64_t q, bool use_q, int radius) {
    auto ball = std::make_unique<CostBall>(Price{t, q}, use_q);
    ball->build(radius);
    return ball;
}

std::vector<size_t> levels(const CostBall& ball) {
    std::vector<size_t> out;
    for (int c = 0; c <= ball.radius(); ++c) out.push_back(ball.level(c).size());
    return out;
}

// Brute-force least cost per coset up to `radius`.
std::map<CosetKey, int64_t> brute_force(int64_t t, int64_t q, bool use_q, int64_t radius) {
    std::vector<std::pair<std::string, int64_t>> letters{{"T", t}};
    if (use_q) letters.push_back({"Q", q});
    std::map<CosetKey, int64_t> dist{{coset_key(word_matrix("")), 0}};
    std::vector<std::pair<int64_t, std::string>> frontier{{0, ""}};
    std::vector<std::string> cliffords;
    for (const auto& [rotation, word] : s_h_spellings()) cliffords.push_back(word);
    for (int64_t c = 0; c <= radius; ++c) {
        std::vector<std::string> here;
        for (const auto& [cost, word] : frontier)
            if (cost == c) here.push_back(word);
        for (const std::string& u : here)
            for (const std::string& cl : cliffords)
                for (const auto& [g, w] : letters) {
                    if (c + w > radius) continue;
                    const std::string v = u + cl + g;
                    const CosetKey k = coset_key(word_matrix(v));
                    const auto it = dist.find(k);
                    if (it == dist.end() || it->second > c + w) {
                        dist[k] = c + w;
                        frontier.push_back({c + w, v});
                    }
                }
    }
    return dist;
}

}  // namespace

TEST_CASE("the T-only ball is the Clifford+T syllable graph") {
    CHECK(levels(*built(1, 3, false, 6)) == std::vector<size_t>{1, 3, 6, 12, 24, 48, 96});
}

TEST_CASE("the weighted ball is the brute force's coset by coset") {
    // prices on both sides of t = q and on Clifford+T; no coset twice
    struct Case {
        int64_t t, q;
        bool use_q;
        int radius;
    };
    for (const Case& k :
         {Case{1, 3, true, 5}, Case{1, 1, true, 4}, Case{2, 3, true, 8}, Case{1, 3, false, 8}}) {
        INFO(k.t, ",", k.q, " use_q ", k.use_q);
        const std::map<CosetKey, int64_t> dist = brute_force(k.t, k.q, k.use_q, k.radius);
        std::map<int64_t, size_t> brute;
        for (const auto& [key, c] : dist) ++brute[c];
        const auto ball_held = built(k.t, k.q, k.use_q, k.radius);
        const CostBall& ball = *ball_held;
        for (int c = 0; c <= k.radius; ++c) {
            const Level& level = ball.level(c);
            CHECK(level.size() == brute[c]);
            std::set<CosetKey> keys;
            for (size_t i = 0; i < level.size(); ++i) {
                const CyclotomicMatrix e = word_matrix(ball.word(c, static_cast<int>(i)));
                const CosetKey key = coset_key(e);
                keys.insert(key);
                CHECK(dist.at(key) == c);
            }
            CHECK(keys.size() == level.size());
        }
    }
}

TEST_CASE("the descent prices every coset of the ball at its level") {
    // with the cheapest word's magic counts
    for (const auto& [price, radius] :
         std::vector<std::pair<Price, int>>{{{1, 3}, 8}, {{1, 1}, 5}, {{2, 3}, 9}, {{2, 5}, 11}}) {
        const auto ball_held = built(price.t, price.q, true, radius);
        const CostBall& ball = *ball_held;
        const Pricer pricer = Pricer::descent(price.t, price.q);
        for (int c = 0; c <= radius; ++c)
            for (size_t i = 0; i < ball.level(c).size(); ++i) {
                const std::string word = ball.word(c, static_cast<int>(i));
                INFO(word);
                const Cost got = pricer.price(word_matrix(word));
                const auto [n_t, n_q] = magic_counts(word);
                CHECK(got.value == c);
                CHECK(got.n_t == n_t);
                CHECK(got.n_q == n_q);
            }
    }
}

TEST_CASE("the ball's words are the reference's") {
    for (const json& b : golden("balls")["balls"]) {
        const int64_t t = b["price"][0], q = b["price"][1];
        const int radius = b["radius"];
        INFO(t, ",", q);
        const auto ball_held = built(t, q, true, radius);
        const CostBall& ball = *ball_held;
        CHECK(json(levels(ball)) == b["levels"]);
        std::vector<std::string> words;
        for (int c = 0; c <= radius; ++c)
            for (size_t i = 0; i < ball.level(c).size(); ++i)
                words.push_back(ball.word(c, static_cast<int>(i)));
        CHECK(json(words) == b["words"]);
    }
}

TEST_CASE("the level sizes are the exact counts and satisfy the recurrence") {
    const auto ball_held = built(1, 3, true, 9);
    const CostBall& ball = *ball_held;
    CHECK(levels(ball) == kExactLevels);
    for (size_t n = 4; n < kExactLevels.size(); ++n)
        CHECK(kExactLevels[n] == 2 * kExactLevels[n - 1] + 4 * kExactLevels[n - 3]);
    size_t total = 0;
    for (size_t n : kExactLevels) total += n;
    CHECK(ball.size() == total);
    CHECK(ball.max_weight() == 3);
}

TEST_CASE("other prices reproduce the known growths") {
    // (2, 3): Bloch lambda-exponent; (1, 1): non-Clifford count
    CHECK(levels(*built(2, 3, true, 10)) ==
          std::vector<size_t>{1, 0, 3, 6, 6, 24, 36, 72, 168, 288, 624});
    CHECK(levels(*built(1, 1, true, 5)) == std::vector<size_t>{1, 9, 54, 324, 1944, 11664});
}

TEST_CASE("the words have their level's cost") {
    const auto ball_held = built(1, 3, true, 9);
    const CostBall& ball = *ball_held;
    Rng rng(0);
    for (int c = 1; c <= 9; ++c) {
        const Level& level = ball.level(c);
        for (int n = 0; n < 20; ++n) {
            const size_t i = static_cast<size_t>(rng.below(static_cast<int64_t>(level.size())));
            const std::string word = ball.word(c, static_cast<int>(i));
            INFO(word);
            const auto [n_t, n_q] = magic_counts(word);
            CHECK(n_t + 3 * n_q == c);
        }
    }
}

TEST_CASE("every coset of cost >= 3 has one of two neighbour profiles") {
    // profiles over 3 T-edges and 6 sqrt(T)-edges give a_n = 2 a_{n-1} + 4 a_{n-3}
    const auto ball_held = built(1, 3, true, 8);
    const CostBall& ball = *ball_held;
    std::map<CosetKey, int> cost_of;
    for (int c = 0; c <= 8; ++c)
        for (size_t i = 0; i < ball.level(c).size(); ++i)
            cost_of[coset_key(word_matrix(ball.word(c, static_cast<int>(i))))] = c;
    using Profile = std::pair<std::vector<int>, std::vector<int>>;
    const Profile type_a{{-1, 1, 1}, {2, 2, 3, 3, 3, 3}}, type_b{{0, 1, 1}, {-3, -2, 3, 3, 3, 3}};
    const std::vector<SyllableGenerator> generators = syllable_generators({1, 3}, true);
    for (int c = 3; c < 6; ++c) {
        std::map<Profile, size_t> counts;
        for (size_t i = 0; i < ball.level(c).size(); ++i) {
            const std::string w = ball.word(c, static_cast<int>(i));
            Profile p;
            for (const SyllableGenerator& g : generators) {
                const int d = cost_of.at(coset_key(word_matrix(w + g.word))) - c;
                (g.weight == 1 ? p.first : p.second).push_back(d);
            }
            std::sort(p.first.begin(), p.first.end());
            std::sort(p.second.begin(), p.second.end());
            ++counts[p];
        }
        CHECK(counts.size() == 2);
        CHECK(counts[type_a] == 2 * ball.level(c - 1).size());
        if (c >= 4) CHECK(counts[type_b] == 4 * ball.level(c - 3).size());
    }
}

TEST_CASE("the syllables of each price are the reference's") {
    for (const json& g : golden("words")["generators"]) {
        const int64_t t = g["price"][0], q = g["price"][1];
        INFO(t, ",", q);
        json got = json::array();
        for (const SyllableGenerator& s : syllable_generators({t, q}, g["use_q"])) {
            // sqrt(T)-only subset, in full-list order
            if (!g["use_t"].get<bool>() && s.word.find('T') != std::string::npos) continue;
            got.push_back({s.word, s.weight, row_of(CyclotomicMatrix{s.matrix, s.sqrt2_exponent})});
        }
        CHECK(got == g["generators"]);
    }
}

TEST_CASE("the cosets of each denominator scale are the reference's") {
    for (const json& c : golden("words")["scale_counts"]) {
        const int radius = c["radius"];
        const auto ball_held = built(c["price"][0], c["price"][1], c["use_q"], radius);
        const CostBall& ball = *ball_held;
        json counts = json::object();
        for (int k = 0; k < 2 * radius + 2; ++k)
            if (const int64_t n = ball.scale_count(k)) counts[std::to_string(k)] = n;
        CHECK(counts == c["counts"]);
    }
}

TEST_CASE("the normal form refuses a price it cannot represent") {
    CHECK_THROWS_AS(CostBall({7, 3}), std::invalid_argument);  // t > 2q is invalid
}

}  // namespace quicksynth::test

// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/data_structures/clifford.hpp"

#include <algorithm>
#include <stdexcept>

namespace quicksynth {

CliffordRotation rotation_product(const CliffordRotation& a, const CliffordRotation& b) {
    CliffordRotation out{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) out[i][j] += a[i][k] * b[k][j];
    return out;
}

CliffordRotation rotation_transpose(const CliffordRotation& a) {
    CliffordRotation out{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) out[i][j] = a[j][i];
    return out;
}

std::optional<CliffordRotation> clifford_rotation(char letter) {
    switch (letter) {
        case 'H': return CliffordRotation{{{0, 0, 1}, {0, -1, 0}, {1, 0, 0}}};
        case 'S': return CliffordRotation{{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}}};
        case 's': return CliffordRotation{{{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}}};
        case 'X': return CliffordRotation{{{1, 0, 0}, {0, -1, 0}, {0, 0, -1}}};
        case 'Y': return CliffordRotation{{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}}};
        case 'Z': return CliffordRotation{{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}};
        default: return std::nullopt;
    }
}

CliffordRotation rotation_of(const std::string& word) {
    CliffordRotation p = kIdentityRotation;
    for (char letter : word)
        if (const std::optional<CliffordRotation> g = clifford_rotation(letter))
            p = rotation_product(p, *g);
    return p;
}

const SpellingTables& spelling_tables() {
    static const SpellingTables t = [] {
        SpellingTables out;
        for (const char* a : {"", "H", "SH"}) out.cosets.emplace_back(a, rotation_of(a));
        for (bool b : {false, true}) {
            CliffordRotation m = b ? rotation_of("X") : kIdentityRotation;
            for (int i = 0; i < 4; ++i) {
                out.sx.emplace_back(m, b);
                m = rotation_product(m, rotation_of("S"));
            }
        }
        out.words.emplace_back(kIdentityRotation, "");
        std::vector<std::pair<std::string, CliffordRotation>> frontier{{"", kIdentityRotation}};
        while (!frontier.empty()) {
            std::vector<std::pair<std::string, CliffordRotation>> next;
            for (const auto& [w, m] : frontier)
                for (char letter : std::string("HSsXZY")) {
                    const CliffordRotation n = rotation_product(m, *clifford_rotation(letter));
                    const bool known = std::any_of(out.words.begin(), out.words.end(),
                                                   [&](const auto& e) { return e.first == n; });
                    if (known) continue;
                    out.words.emplace_back(n, w + letter);
                    next.emplace_back(w + letter, n);
                }
            frontier = std::move(next);
        }
        return out;
    }();
    return t;
}

CosetSplit split_coset(const CliffordRotation& p) {
    const SpellingTables& t = spelling_tables();
    for (const auto& [representative, a] : t.cosets) {
        const CliffordRotation d = rotation_product(rotation_transpose(a), p);
        const auto it =
            std::find_if(t.sx.begin(), t.sx.end(), [&](const auto& e) { return e.first == d; });
        if (it != t.sx.end()) return {representative, d, it->second};
    }
    throw std::logic_error("split_coset: a Clifford outside every coset of <S, X>");
}

const std::string& clifford_word(const CliffordRotation& p) {
    const auto& words = spelling_tables().words;
    const auto it =
        std::find_if(words.begin(), words.end(), [&](const auto& e) { return e.first == p; });
    if (it == words.end()) throw std::logic_error("clifford_word: a Clifford without a word");
    return it->second;
}

const std::vector<std::pair<CliffordRotation, std::string>>& s_h_spellings() {
    static const std::vector<std::pair<CliffordRotation, std::string>> out = [] {
        const CliffordRotation H = rotation_of("H");
        const CliffordRotation S = rotation_of("S");
        CliffordRotation Sp[4] = {kIdentityRotation};
        for (int p = 1; p < 4; ++p) Sp[p] = rotation_product(Sp[p - 1], S);
        const auto s = [](int n) { return std::string(static_cast<size_t>(n), 'S'); };
        // layer h appends H S^k to each spelling of layer h - 1, k innermost
        std::vector<std::pair<CliffordRotation, std::string>> layer;
        for (int l = 0; l < 4; ++l) layer.emplace_back(Sp[l], s(l));
        std::vector<std::pair<CliffordRotation, std::string>> v = layer;
        for (int h = 1; h < 4; ++h) {
            std::vector<std::pair<CliffordRotation, std::string>> next;
            for (const auto& [p, w] : layer)
                for (int k = 0; k < 4; ++k)
                    next.emplace_back(rotation_product(rotation_product(p, H), Sp[k]),
                                      w + "H" + s(k));
            v.insert(v.end(), next.begin(), next.end());
            layer = std::move(next);
        }
        // the first spelling of each Clifford
        std::vector<std::pair<CliffordRotation, std::string>> first;
        for (auto& spelling : v) {
            const CliffordRotation& p = spelling.first;
            if (std::none_of(first.begin(), first.end(),
                             [&](const auto& f) { return f.first == p; }))
                first.push_back(std::move(spelling));
        }
        if (first.size() != 24) throw std::logic_error("s_h_spellings: not 24 Cliffords");
        return first;
    }();
    return out;
}

}  // namespace quicksynth

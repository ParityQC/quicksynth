// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/exact_synthesis/normal_form.hpp"

#include "quicksynth/data_structures/int_mp.hpp"
#include "quicksynth/data_structures/clifford.hpp"

#include <optional>
#include <stdexcept>
#include <utility>

namespace quicksynth {

namespace {

bool is_magic(char letter) {
    return letter == 'T' || letter == 't' || letter == 'Q' || letter == 'q';
}

char swapcase(char letter) {
    return letter >= 'a' ? static_cast<char>(letter - 'a' + 'A')
                         : static_cast<char>(letter - 'A' + 'a');
}

// Non-overlapping replace, left to right, like Python's str.replace.
std::string replace_all(const std::string& s, const std::string& from, const std::string& to) {
    std::string out;
    size_t pos = 0;
    for (size_t hit; (hit = s.find(from, pos)) != std::string::npos; pos = hit + from.size()) {
        out.append(s, pos, hit - pos);
        out += to;
    }
    out.append(s, pos, std::string::npos);
    return out;
}

}  // namespace

std::string merge_letters(std::string word) {
    static const std::pair<const char*, const char*> kMerges[] = {
        {"SSSS", ""}, {"HH", ""},  {"XX", ""},  {"YY", ""},
        {"ZZ", ""},   {"TT", "S"}, {"QQ", "T"}, {"SS", "Z"}};
    static const std::pair<const char*, const char*> kReorderings[] = {
        {"ST", "TS"}, {"SQ", "QS"}, {"SZ", "ZS"}, {"QT", "TQ"}, {"ZT", "TZ"}, {"ZQ", "QZ"}};
    for (;;) {
        std::string next = word;
        for (const auto& [from, to] : kMerges) next = replace_all(next, from, to);
        for (const auto& [from, to] : kReorderings) next = replace_all(next, from, to);
        if (next == word) return word;
        word = std::move(next);
    }
}

std::string t_as_qq(const std::string& word) { return replace_all(word, "T", "QQ"); }

std::string normal_form(const std::string& word) {
    std::string out;
    CliffordRotation p = kIdentityRotation;  // the Clifford since the last magic letter
    for (char letter : word) {
        if (!is_magic(letter)) {
            const std::optional<CliffordRotation> g = clifford_rotation(letter);
            if (!g) throw unknown_letter(word);
            p = rotation_product(p, *g);
            continue;
        }
        // P = A D: A is written, D is carried past the letter, which an X in D inverts
        const CosetSplit split = split_coset(p);
        out += split.representative;
        out += split.holds_x ? swapcase(letter) : letter;
        p = split.remainder;
    }
    out += clifford_word(p);
    if (!equal_up_to_phase(word_matrix_mp(out), word_matrix_mp(word)))
        throw std::logic_error("normal_form changed the gate: " + word + " -> " + out);
    return out;
}

}  // namespace quicksynth

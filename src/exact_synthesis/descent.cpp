// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/exact_synthesis/descent.hpp"

#include "quicksynth/exact_synthesis/normal_form.hpp"
#include "quicksynth/data_structures/bloch.hpp"
#include "quicksynth/data_structures/int_mp.hpp"
#include "quicksynth/data_structures/clifford.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace quicksynth {
namespace {

// canonicalise()'s column orders (new -> old); Clifford words in kConnectors.
constexpr int kColumnOrders[3][3] = {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}};

// Matrix with reordered columns; order indexes kColumnOrders and kConnectors.
template <class D>
struct Canonical {
    LambdaMatrix<D> m;
    int order = 0;
};

// Puts the two columns not divisible by lambda (an odd constant digit) first.
template <class D>
Canonical<D> canonicalise(const LambdaMatrix<D>& m) {
    bool odd[3] = {false, false, false};
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r)
            if (!is_even(at(m, r, c)[0])) odd[c] = true;
    Canonical<D> out;
    if (odd[0] && odd[1] && !odd[2])
        out.order = 0;
    else if (odd[0] && odd[2] && !odd[1])
        out.order = 1;
    else if (odd[1] && odd[2] && !odd[0])
        out.order = 2;
    else
        throw std::invalid_argument("descent: not exactly two columns off the lambda ideal");
    const int* columns = kColumnOrders[out.order];
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) at(out.m, r, c) = at(m, r, columns[c]);
    return out;
}

template <class D>
LambdaMatrix<D> next_syllable(const LambdaMatrix<D>& m, int s) {
    LambdaMatrix<D> out;
    for (int r = 0; r < 3; ++r) {
        const Lambda<D>& x = at(m, r, 0);
        const Lambda<D>& y = at(m, r, 1);
        at(out, r, 0) = times_u(add(x, times_theta(times_sign(y, s))));
        at(out, r, 1) = times_sign(times_u(sub(y, times_theta(times_sign(x, s)))), s);
        at(out, r, 2) = times_sign(times_lambda_pow(at(m, r, 2), 3), s);
    }
    return out;
}

// Applies a Q syllable of sign s to a canonical image; none if it ascends.
template <class D>
std::optional<BlochImage<D>> descend(const BlochImage<D>& image, int s) {
    BlochImage<D> out;
    out.m = next_syllable(image.m, s);
    const int64_t d = strip_lambda(out.m);
    if (d % 4 != 2) return std::nullopt;
    out.denominator_exponent = image.denominator_exponent + 3 - d;
    return out;
}

// Normal-form syllable: column order of its Clifford prefix (kConnectors),
// sign (Q for +1, X Q X for -1), and T (Q Q on one axis) or a single Q.
struct Syllable {
    int order = 0;
    int sign = 1;
    bool is_t = false;
};

// Descends an image with no common lambda factor, passing each syllable to
// visit; returns the remaining exponent-0 image, a Clifford.
template <class D, class Visit>
BlochImage<D> walk_normal_form(BlochImage<D> image, Visit&& visit) {
    const int64_t limit = 64 + 4 * image.denominator_exponent;
    for (int64_t steps = 0; image.denominator_exponent > 0; ++steps) {
        if (steps > limit) throw std::invalid_argument("descent: the exponent does not fall");
        // new syllable on any axis: reorder the columns
        const Canonical<D> c = canonicalise(image.m);
        const BlochImage<D> canonical{c.m, image.denominator_exponent};
        Syllable syllable{c.order, 1, false};
        std::optional<BlochImage<D>> one = descend(canonical, 1);
        if (!one) {
            syllable.sign = -1;
            one = descend(canonical, -1);
        }
        if (!one) throw std::invalid_argument("descent: no descending syllable");
        // Q^-1 = X Q X: descend(m, -1) is descend(m X, +1); carry the trailing X
        if (syllable.sign == -1) apply_x(one->m);
        // same axis and sign again: Q Q = T, or X Q X X Q X = T^-1
        std::optional<BlochImage<D>> two = descend(*one, syllable.sign);
        if (two) {
            if (syllable.sign == -1) apply_x(two->m);
            syllable.is_t = true;
            image = std::move(*two);
        } else {
            image = std::move(*one);
        }
        visit(syllable);
    }
    return image;
}

// Normal-form weight of the image.
template <class D>
Cost normal_form_cost(BlochImage<D> image, int64_t t_cost, int64_t q_cost) {
    int64_t n_t = 0, n_q = 0;
    walk_normal_form(std::move(image), [&](const Syllable& syllable) {
        if (syllable.is_t)
            ++n_t;
        else
            ++n_q;
    });
    if (2 * q_cost <= t_cost) {  // Q Q costs no more than T
        n_q += 2 * n_t;
        n_t = 0;
    }
    return Price{t_cost, q_cost}.cost(n_t, n_q);
}

// --- normal-form words -------------------------------------------------------------

// Inverse up to phase: letters reversed, T -> X T X (one T, unlike S^3 T),
// Q -> X Q X, S -> S S S.
std::string get_inverse(const std::string& sequence) {
    std::string out;
    for (char letter : sequence) {
        switch (letter) {
            case 'T': out = "XTX" + out; break;
            case 'Q': out = "XQX" + out; break;
            case 'S': out = "SSS" + out; break;
            case 'X':
            case 'Y':
            case 'Z':
            case 'H': out = std::string(1, letter) + out; break;
            default: break;
        }
    }
    return out;
}

// Descent::cost for a CyclotomicMatrixMp.
Cost cost_of(const CyclotomicMatrixMp& element, int64_t t_cost, int64_t q_cost) {
    return normal_form_cost(bloch_image(element), t_cost, q_cost);
}

// Word for the element up to global phase, every T spelt T.
std::string exact_synthesis(const CyclotomicMatrixMp& element) {
    // per syllable: Clifford prefix, then Q or X Q X, twice for a T
    std::string sequence;
    const BlochImage<int_mp> leftover =
        walk_normal_form(bloch_image(element), [&](const Syllable& syllable) {
            const char* q = syllable.sign == -1 ? "XQX" : "Q";
            sequence += kConnectors[syllable.order];
            sequence += q;
            if (syllable.is_t) sequence += q;
        });
    // leftover Clifford as a signed permutation
    CliffordRotation c{};
    for (int r = 0; r < 3; ++r)
        for (int col = 0; col < 3; ++col) {
            const Lambda<int_mp>& x = at(leftover.m, r, col);
            if (!x[1].zero() || !x[2].zero() || !x[3].zero() ||
                !(x[0].zero() || x[0] == 1 || x[0] == -1))
                throw std::invalid_argument("descent: the remainder is not a Clifford");
            c[r][col] = static_cast<int>(x[0].to_int64());
        }
    const auto& spellings = s_h_spellings();
    const auto it = std::find_if(spellings.begin(), spellings.end(),
                                 [&](const auto& sp) { return sp.first == c; });
    if (it == spellings.end())
        throw std::invalid_argument("descent: the remainder is not a Clifford");
    std::string w = merge_letters(it->second + get_inverse(merge_letters(sequence)));
    // verify: the word reproduces the element up to phase
    CyclotomicMatrixMp built;
    if (!word_matrix_mp(w, built) || !equal_up_to_phase(element, built))
        throw std::logic_error("descent: the word does not reproduce the element");
    return w;
}

// exact_synthesis's word, with T spelt Q Q when 2 q_cost <= t_cost.
std::string price_word(const CyclotomicMatrixMp& element, int64_t t_cost, int64_t q_cost) {
    std::string w = exact_synthesis(element);
    if (2 * q_cost <= t_cost) w = t_as_qq(w);
    return w;
}

}  // namespace

Cost Descent::cost(BlochImage128 image) const {
    image.denominator_exponent -= strip_lambda(image.m);
    return normal_form_cost(std::move(image), t_, q_);
}

Cost Descent::cost(const M2& u, int k) const { return cost_of(widen(u, k), t_, q_); }

Cost Descent::cost(const CyclotomicMatrixDigits& u) const {
    return cost_of(from_digits(u), t_, q_);
}

std::string Descent::word(const M2& u, int k) const { return price_word(widen(u, k), t_, q_); }

std::string Descent::word(const CyclotomicMatrixDigits& u) const {
    return price_word(from_digits(u), t_, q_);
}

}  // namespace quicksynth

// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "fixtures.hpp"

#include "quicksynth/data_structures/int_mp.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace quicksynth::test {

namespace {

// nlohmann_json DOM parser that keeps integers wider than 64 bits as digit
// strings instead of rounding them to double.
struct WideIntegers : nlohmann::detail::json_sax_dom_parser<json> {
    using json_sax_dom_parser::json_sax_dom_parser;
    bool number_float(double value, const std::string& text) {
        if (text.find_first_of(".eE") == std::string::npos) {
            std::string digits = text;
            return string(digits);
        }
        return json_sax_dom_parser::number_float(value, text);
    }
};

// Quotes Python's bare Infinity, -Infinity and NaN outside strings.
std::string quote_bare_words(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    bool in_string = false;
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (in_string) {
            out += c;
            if (c == '\\') {
                out += text[++i];
            } else if (c == '"') {
                in_string = false;
            }
            continue;
        }
        if (c == '"') in_string = true;
        bool replaced = false;
        for (const char* word : {"-Infinity", "Infinity", "NaN"}) {
            const std::string w(word);
            if (text.compare(i, w.size(), w) == 0) {
                out += '"' + w + '"';
                i += w.size() - 1;
                replaced = true;
                break;
            }
        }
        if (!replaced) out += c;
    }
    return out;
}

json parse(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::stringstream text;
    text << in.rdbuf();
    json out;
    WideIntegers sax(out);
    json::sax_parse(quote_bare_words(text.str()), &sax);
    return out;
}

int_mp integer(const std::string& digits, int base) {
    int_mp x;
    if (!int_mp::parse(digits.c_str(), base, x))
        throw std::invalid_argument("not an integer in base " + std::to_string(base) + ": " +
                                    digits);
    return x;
}

std::string convert_base(const std::string& digits, int from, int to) {
    return integer(digits, from).str(to);
}

}  // namespace

const json& repository_json(const std::string& relative) {
    static std::mutex mutex;
    static std::map<std::string, std::unique_ptr<json>> documents;
    std::lock_guard<std::mutex> lock(mutex);
    std::unique_ptr<json>& slot = documents[relative];
    if (!slot)
        slot = std::make_unique<json>(parse(std::string(QUICKSYNTH_REPOSITORY) + "/" + relative));
    return *slot;
}

const json& golden(const std::string& group) {
    return repository_json("tests/fixtures/golden_" + group + ".json");
}

double number(const json& value) {
    if (value.is_string()) {
        const std::string& s = value.get_ref<const std::string&>();
        if (s == "Infinity") return std::numeric_limits<double>::infinity();
        if (s == "-Infinity") return -std::numeric_limits<double>::infinity();
        if (s == "NaN") return std::numeric_limits<double>::quiet_NaN();
        return std::stod(s);
    }
    return value.get<double>();
}

std::string decimal(const json& value) {
    if (value.is_string()) return value.get<std::string>();
    if (value.is_number_unsigned()) return std::to_string(value.get<uint64_t>());
    return std::to_string(value.get<int64_t>());
}

std::optional<int64_t> optional_of(const json& value) {
    if (value.is_null()) return std::nullopt;
    return value.get<int64_t>();
}

Unitary decode_target(const json& target) {
    Unitary out{};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            out[4 * i + 2 * j] = number(target["re"][i][j]);
            out[4 * i + 2 * j + 1] = number(target["im"][i][j]);
        }
    return out;
}

std::vector<Unitary> decode_targets(const json& targets) {
    std::vector<Unitary> out;
    for (const json& t : targets) out.push_back(decode_target(t));
    return out;
}

const std::vector<Unitary>& calibration_targets() {
    static const std::vector<Unitary> targets =
        decode_targets(repository_json("tests/fixtures/calibration_targets.json"));
    return targets;
}

GateSet ring_gateset(const std::string& ring) {
    if (ring == "sqrtT") return GateSet::sqrt_t();
    if (ring == "T") return GateSet::clifford_t();
    throw std::invalid_argument("unknown ring " + ring);
}

CyclotomicMatrix matrix_of_row(const json& row) {
    CyclotomicMatrix e;
    e.denominator_exponent = row.at(0).get<int>();
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int t = 0; t < kZetaDegree; ++t)
                e.matrix[a][b][t] = row.at(1 + 16 * a + 8 * b + t);
    return e;
}

json row_of(const CyclotomicMatrix& element) {
    json row = json::array({element.denominator_exponent});
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int t = 0; t < kZetaDegree; ++t) row.push_back(element.matrix[a][b][t]);
    return row;
}

CyclotomicMatrixDigits digits_of_row(const json& row) {
    CyclotomicMatrixDigits e;
    e.denominator_exponent = row.at(0).get<int>();
    for (size_t i = 0; i < 32; ++i) e.hex[i] = convert_base(decimal(row.at(1 + i)), 10, 16);
    return e;
}

std::vector<std::string> decimal_row(const CyclotomicMatrixDigits& element) {
    std::vector<std::string> out{std::to_string(element.denominator_exponent)};
    for (const std::string& h : element.hex) out.push_back(convert_base(h, 16, 10));
    return out;
}

CyclotomicMatrixDigits digits_of(const CyclotomicMatrix& element) {
    return digits_of_row(row_of(element));
}

CyclotomicMatrixDigits widened(const CyclotomicMatrix& element, int bits) {
    CyclotomicMatrixDigits out = digits_of(element);
    out.denominator_exponent += 2 * bits;
    for (std::string& h : out.hex) h = shift_left(integer(h, 16), bits).str(16);
    return out;
}

}  // namespace quicksynth::test

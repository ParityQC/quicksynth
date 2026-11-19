// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// int_mp against values computed with Python ints, and its conversions at limb
// and rounding boundaries.

#include "quicksynth/data_structures/int_mp.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace quicksynth::test {

namespace {

int_mp hex(const char* text) {
    int_mp x;
    REQUIRE(int_mp::parse(text, 16, x));
    return x;
}

int_mp power_of_two(int bits) { return shift_left(int_mp(1), bits); }

// a, b, s; then a + b, a - b, a b, a s, a / 2 truncated, and a in decimal.
struct Case {
    const char* a;
    const char* b;
    int64_t s;
    const char* sum;
    const char* difference;
    const char* product;
    const char* scaled;
    const char* half;
    const char* decimal;
};

// Python ints at limb boundaries and of random sizes up to 200 bits.
const Case kCases[] = {
    {"0", "7fffffffffffffff", 0, "7fffffffffffffff", "-7fffffffffffffff", "0", "0", "0", "0"},
    {"1", "-80000000000000000000000000000000", INT64_MIN, "-7fffffffffffffffffffffffffffffff",
     "80000000000000000000000000000001", "-80000000000000000000000000000000", "-8000000000000000",
     "0", "1"},
    {"-1", "-1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a61", -3,
     "-1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a62",
     "1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a60",
     "1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a61", "3", "0", "-1"},
    {"7fffffffffffffff", "0", 1, "7fffffffffffffff", "7fffffffffffffff", "0", "7fffffffffffffff",
     "3fffffffffffffff", "9223372036854775807"},
    {"-8000000000000000", "10000000000000000", -1, "8000000000000000", "-18000000000000000",
     "-80000000000000000000000000000000", "8000000000000000", "-4000000000000000",
     "-9223372036854775808"},
    {"8000000000000000", "-ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 7,
     "-ffffffffffffffffffffffffffffffffffffffffffffffff7fffffffffffffff",
     "10000000000000000000000000000000000000000000000007fffffffffffffff",
     "-7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff8000000000000000",
     "38000000000000000", "4000000000000000", "9223372036854775808"},
    {"ffffffffffffffff", "24d3ac94af0f21ddb6", 1, "25d3ac94af0f21ddb5", "-23d3ac94af0f21ddb7",
     "24d3ac94af0f21dd912c536b50f0de224a", "ffffffffffffffff", "7fffffffffffffff",
     "18446744073709551615"},
    {"10000000000000000", "-8000000000000000", INT64_MIN, "8000000000000000", "18000000000000000",
     "-80000000000000000000000000000000", "-80000000000000000000000000000000", "8000000000000000",
     "18446744073709551616"},
    {"-10000000000000000", "100000000000000000000000000000000", 0,
     "ffffffffffffffff0000000000000000", "-100000000000000010000000000000000",
     "-1000000000000000000000000000000000000000000000000", "0", "-8000000000000000",
     "-18446744073709551616"},
    {"7fffffffffffffffffffffffffffffff", "-224be40ad23f0824128b2f330c5c7fd0a6a3a4506513270e",
     INT64_MIN, "-224be40ad23f0823928b2f330c5c7fd0a6a3a4506513270f",
     "224be40ad23f0824928b2f330c5c7fd0a6a3a4506513270d",
     "-1125f205691f841209459799862e3fe83105ee1d604a8b62ed74d0ccf3a3802f595c5baf9aecd8f2",
     "-3fffffffffffffffffffffffffffffff8000000000000000", "3fffffffffffffffffffffffffffffff",
     "170141183460469231731687303715884105727"},
    {"-80000000000000000000000000000000", "1", -1, "-7fffffffffffffffffffffffffffffff",
     "-80000000000000000000000000000001", "-80000000000000000000000000000000",
     "80000000000000000000000000000000", "-40000000000000000000000000000000",
     "-170141183460469231731687303715884105728"},
    {"100000000000000000000000000000000", "-10000000000000000", INT64_MIN,
     "ffffffffffffffff0000000000000000", "100000000000000010000000000000000",
     "-1000000000000000000000000000000000000000000000000",
     "-800000000000000000000000000000000000000000000000", "80000000000000000000000000000000",
     "340282366920938463463374607431768211456"},
    {"ffffffffffffffffffffffffffffffff",
     "10000000000000000000000000000000000000000000000000000000000000000", -3,
     "100000000000000000000000000000000ffffffffffffffffffffffffffffffff",
     "-ffffffffffffffffffffffffffffffff00000000000000000000000000000001",
     "ffffffffffffffffffffffffffffffff0000000000000000000000000000000000000000000000000000000000000000",
     "-2fffffffffffffffffffffffffffffffd", "7fffffffffffffffffffffffffffffff",
     "340282366920938463463374607431768211455"},
    {"100000000000000000000000000000000ffffffffffffffff", "-25a09f76b5a170b338",
     9223372036854775807, "ffffffffffffffffffffffffffffffdb5f60894a5e8f4cc7",
     "100000000000000000000000000000026a09f76b5a170b337",
     "-25a09f76b5a170b3380000000000000025a09f76b5a170b3125f60894a5e8f4cc8",
     "7fffffffffffffff00000000000000007ffffffffffffffe8000000000000001",
     "800000000000000000000000000000007fffffffffffffff",
     "6277101735386680763835789423207666416120802188537744064511"},
    {"-ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", "8000000000000000", 1,
     "-ffffffffffffffffffffffffffffffffffffffffffffffff7fffffffffffffff",
     "-10000000000000000000000000000000000000000000000007fffffffffffffff",
     "-7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff8000000000000000",
     "-ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
     "-7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
     "-115792089237316195423570985008687907853269984665640564039457584007913129639935"},
    {"10000000000000000000000000000000000000000000000000000000000000000",
     "ffffffffffffffffffffffffffffffff", 0,
     "100000000000000000000000000000000ffffffffffffffffffffffffffffffff",
     "ffffffffffffffffffffffffffffffff00000000000000000000000000000001",
     "ffffffffffffffffffffffffffffffff0000000000000000000000000000000000000000000000000000000000000000",
     "0", "8000000000000000000000000000000000000000000000000000000000000000",
     "115792089237316195423570985008687907853269984665640564039457584007913129639936"},
    {"1fd5863c3eb0469ec21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ae225bfaff1eaaf8b0a1",
     "3a0ed904759531985d", INT64_MIN,
     "1fd5863c3eb0469ec21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ae5c6ad40394402a48fe",
     "1fd5863c3eb0469ec21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ade84d21faa915c71844",
     "73839138a3dab83b6072a49e12435356cf6e88d60dcd43939ace570cd1b8426960a7049b394cef940acee99a44908c27d",
     "-feac31e1f58234f610d49bd3b79a197feb9ecbf223b035b41f67b7b7253d7112dfd7f8f557c58508000000000000000",
     "feac31e1f58234f610d49bd3b79a197feb9ecbf223b035b41f67b7b7253d7112dfd7f8f557c5850",
     "265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699044001"},
    {"-1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a61", "-1", INT64_MIN,
     "-1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a62",
     "-1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a60",
     "1aa3b2c5319d5e494c9a977611d99b7b5cb34b967d4a2c6aecef68933be1fc93d3a1a61",
     "d51d96298ceaf24a64d4bbb08eccdbdae59a5cb3ea516357677b4499df0fe49e9d0d308000000000000000",
     "-d51d96298ceaf24a64d4bbb08eccdbdae59a5cb3ea516357677b4499df0fe49e9d0d30",
     "-3234476509624757991344647769100216810857203198904625400933895331391691459636928060001"},
    {"-224be40ad23f0824128b2f330c5c7fd0a6a3a4506513270e", "7fffffffffffffffffffffffffffffff",
     9223372036854775807, "-224be40ad23f0823928b2f330c5c7fd0a6a3a4506513270f",
     "-224be40ad23f0824928b2f330c5c7fd0a6a3a4506513270d",
     "-1125f205691f841209459799862e3fe83105ee1d604a8b62ed74d0ccf3a3802f595c5baf9aecd8f2",
     "-1125f205691f8411e6f9b38eb3ef37c440c6a2f5262d13b6595c5baf9aecd8f2",
     "-1125f205691f841209459799862e3fe85351d22832899387",
     "-840946467821730214414669612260512535109722854322804565774"},
    {"3a0ed904759531985d",
     "1fd5863c3eb0469ec21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ae225bfaff1eaaf8b0a1", 1,
     "1fd5863c3eb0469ec21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ae5c6ad40394402a48fe",
     "-1fd5863c3eb0469ec21a937a76f3432ffd73d97e447606b683ecf6f6e4a7ade84d21faa915c71844",
     "73839138a3dab83b6072a49e12435356cf6e88d60dcd43939ace570cd1b8426960a7049b394cef940acee99a44908c27d",
     "3a0ed904759531985d", "1d076c823aca98cc2e", "1070981047564691937373"},
    {"-11e20b8f6b0d549b6f03675a1600a35a", "ef9ebdacc0cb1e29c", -1,
     "-11e20b8f6b0d548c75178c8e094ec0be", "-11e20b8f6b0d54aa68ef422622b285f6",
     "-10bd1f8dfb55e3a23bfa5f30ed95783d5e6dc902321d2fed8", "11e20b8f6b0d549b6f03675a1600a35a",
     "-8f105c7b586aa4db781b3ad0b0051ad", "-23770569487636332224155707499473838938"},
    {"24d3ac94af0f21ddb6", "ffffffffffffffff", 0, "25d3ac94af0f21ddb5", "23d3ac94af0f21ddb7",
     "24d3ac94af0f21dd912c536b50f0de224a", "0", "1269d64a578790eedb", "679335516171134950838"},
    {"-25a09f76b5a170b338", "100000000000000000000000000000000ffffffffffffffff", INT64_MIN,
     "ffffffffffffffffffffffffffffffdb5f60894a5e8f4cc7",
     "-100000000000000000000000000000026a09f76b5a170b337",
     "-25a09f76b5a170b3380000000000000025a09f76b5a170b3125f60894a5e8f4cc8",
     "12d04fbb5ad0b8599c0000000000000000", "-12d04fbb5ad0b8599c", "-694103630817088549688"},
    {"ef9ebdacc0cb1e29c", "-11e20b8f6b0d549b6f03675a1600a35a", 9223372036854775807,
     "-11e20b8f6b0d548c75178c8e094ec0be", "11e20b8f6b0d54aa68ef422622b285f6",
     "-10bd1f8dfb55e3a23bfa5f30ed95783d5e6dc902321d2fed8", "77cf5ed660658f13f06142533f34e1d64",
     "77cf5ed660658f14e", "276263145136825950876"},
};

}  // namespace

TEST_CASE("int_mp arithmetic matches Python ints") {
    for (const Case& c : kCases) {
        CAPTURE(c.a);
        CAPTURE(c.b);
        const int_mp a = hex(c.a), b = hex(c.b);
        CHECK((a + b).str(16) == c.sum);
        CHECK((a - b).str(16) == c.difference);
        CHECK((a * b).str(16) == c.product);
        CHECK((b * a).str(16) == c.product);
        CHECK((a * c.s).str(16) == c.scaled);
        CHECK((c.s * a).str(16) == c.scaled);
        CHECK(half_of(a).str(16) == c.half);
        CHECK(a.str(10) == c.decimal);
        CHECK(a.str(16) == c.a);
        CHECK((-(-a)) == a);
        CHECK((a - b) == -(b - a));
        int_mp x = a;
        x += b;
        CHECK(x == hex(c.sum));
        x -= b;
        CHECK(x == a);
        x -= b;
        CHECK(x == hex(c.difference));
    }
}

TEST_CASE("int_mp in-place operations may alias their operand") {
    for (const Case& c : kCases) {
        int_mp x = hex(c.a);
        x += x;
        CHECK(x == hex(c.a) * 2);
        x -= x;
        CHECK(x.zero());
        CHECK(x.sign() == 0);
        CHECK(x.str(16) == "0");
    }
}

TEST_CASE("int_mp construction and narrowing at the 64- and 128-bit limits") {
    const int64_t lo = std::numeric_limits<int64_t>::min();
    const int64_t hi = std::numeric_limits<int64_t>::max();
    CHECK(int_mp(lo).str(16) == "-8000000000000000");
    CHECK(int_mp(hi).str(10) == "9223372036854775807");
    CHECK(int_mp(lo).fits_int64());
    CHECK(int_mp(hi).fits_int64());
    CHECK(int_mp(lo).to_int64() == lo);
    CHECK(int_mp(hi).to_int64() == hi);
    CHECK_FALSE((int_mp(hi) + int_mp(1)).fits_int64());
    CHECK_FALSE((int_mp(lo) - int_mp(1)).fits_int64());
    CHECK(int_mp(lo) == lo);
    CHECK_FALSE(power_of_two(64) == 0);
    CHECK(int_mp(int8_t{-5}) == -5);
    const __int128 big = static_cast<__int128>(~(static_cast<unsigned __int128>(1) << 127));
    CHECK(int_mp(big).str(16) == "7fffffffffffffffffffffffffffffff");
    CHECK(int_mp(-big - 1).str(16) == "-80000000000000000000000000000000");
    CHECK(int_mp(-big - 1) == -power_of_two(127));
    CHECK(int_mp(__int128{0}).zero());
    CHECK(int_mp(big).to_int128() == big);
    CHECK(int_mp(-big).to_int128() == -big);
    CHECK(int_mp(lo).to_int128() == lo);
    CHECK((power_of_two(100) + int_mp(3)).to_int128() == (static_cast<__int128>(1) << 100) + 3);
}

TEST_CASE("int_mp sign parity and bit length") {
    CHECK(int_mp().zero());
    CHECK(int_mp().even());
    CHECK(int_mp().bits() == 0);
    CHECK(int_mp(-3).sign() == -1);
    CHECK(int_mp(3).sign() == 1);
    CHECK_FALSE(int_mp(-3).even());
    CHECK(int_mp(-4).even());
    CHECK(int_mp(-1).bits() == 1);
    CHECK(power_of_two(63).bits() == 64);
    CHECK(power_of_two(64).bits() == 65);
    CHECK((power_of_two(128) - int_mp(1)).bits() == 128);
    CHECK(power_of_two(5000).bits() == 5001);
    CHECK((power_of_two(64) + int_mp(1)).even() == false);
}

TEST_CASE("int_mp halving truncates toward zero") {
    CHECK(half_of(int_mp(-1)).zero());
    CHECK(half_of(int_mp(-1)).sign() == 0);
    CHECK(half_of(int_mp(-3)) == -1);
    CHECK(half_of(int_mp(7)) == 3);
    CHECK(half_of(-(power_of_two(64) + int_mp(1))) == -power_of_two(63));
    CHECK(half_of(power_of_two(64)) == power_of_two(63));
    CHECK(half_of(power_of_two(200)) == power_of_two(199));
}

TEST_CASE("int_mp shifts cross limb boundaries") {
    CHECK(shift_left(int_mp(1), 0) == 1);
    CHECK(shift_left(int_mp(-3), 63).str(16) == "-18000000000000000");
    CHECK(shift_left(int_mp(-3), 64).str(16) == "-30000000000000000");
    CHECK(shift_left(hex("ffffffffffffffff"), 1).str(16) == "1fffffffffffffffe");
    CHECK(shift_left(int_mp(), 100).zero());
}

TEST_CASE("int_mp parses signed digits and skips white space") {
    int_mp x;
    CHECK(int_mp::parse("  -1F", 16, x));
    CHECK(x == -31);
    CHECK(int_mp::parse("1 0", 10, x));
    CHECK(x == 10);
    CHECK(int_mp::parse("-0", 16, x));
    CHECK(x.zero());
    CHECK(x.sign() == 0);
    CHECK(int_mp::parse("0000000000000000000000000000000000000012", 10, x));
    CHECK(x == 12);
    x = int_mp(5);
    CHECK_FALSE(int_mp::parse("", 16, x));
    CHECK_FALSE(int_mp::parse("-", 16, x));
    CHECK_FALSE(int_mp::parse("- 1", 16, x));
    CHECK_FALSE(int_mp::parse("+1", 16, x));
    CHECK_FALSE(int_mp::parse("0x1", 16, x));
    CHECK_FALSE(int_mp::parse("1g", 16, x));
    CHECK_FALSE(int_mp::parse("a", 10, x));
    CHECK(x == 5);
}

TEST_CASE("int_mp holds thousands of bits") {
    // 3^3000: 1432 decimal digits, 231080957811190927... ...769798453765660001
    int_mp x(1);
    for (int i = 0; i < 3000; ++i) x = x * 3;
    const std::string text = x.str(10);
    CHECK(text.size() == 1432);
    CHECK(text.substr(0, 18) == "231080957811190927");
    CHECK(text.substr(text.size() - 18) == "769798453765660001");
    int_mp y(1);
    for (int i = 0; i < 1500; ++i) y = y * int_mp(3);
    CHECK(y * y == x);
    CHECK(hex(x.str(16).c_str()) == x);
    int_mp z;
    REQUIRE(int_mp::parse(text.c_str(), 10, z));
    CHECK(z == x);
    CHECK(half_of(x * 2) == x);
    CHECK((x - y * y).zero());
    CHECK((-x + x).sign() == 0);
}

TEST_CASE("int_mp to_double rounds to 53 bits then to a double") {
    CHECK(int_mp().to_double(10) == 0.0);
    CHECK(int_mp(-5).to_double(-1) == -2.5);
    const int_mp two53 = power_of_two(53);
    // ties to even, and a set bit below the half rounds up
    CHECK((two53 + int_mp(1)).to_double(0) == 9007199254740992.0);
    CHECK((two53 + int_mp(3)).to_double(0) == 9007199254740996.0);
    CHECK((shift_left(two53 + int_mp(1), 1) + int_mp(1)).to_double(-1) == 9007199254740994.0);
    CHECK((-(two53 + int_mp(3))).to_double(0) == -9007199254740996.0);
    // a carry to the next power of two
    CHECK((power_of_two(54) - int_mp(1)).to_double(0) == 18014398509481984.0);
    CHECK((power_of_two(64) - int_mp(1)).to_double(-64) == 1.0);
    CHECK(power_of_two(1023).to_double(0) == std::ldexp(1.0, 1023));
    CHECK(power_of_two(1024).to_double(0) == std::numeric_limits<double>::infinity());
    CHECK(power_of_two(5000).to_double(-4990) == 1024.0);
    // below the normal range the 53-bit value is rounded again: 2^54 + 5 rounds
    // to (2^52 + 1) 4, whose last bit is then a tie to even
    CHECK((power_of_two(54) + int_mp(5)).to_double(-1077) == std::ldexp(1.0, -1023));
}

TEST_CASE("int_mp truncated_frexp truncates to 53 bits") {
    int e = -1;
    CHECK(int_mp().truncated_frexp(e) == 0.0);
    CHECK(e == 0);
    CHECK(int_mp(-3).truncated_frexp(e) == -0.75);
    CHECK(e == 2);
    CHECK((power_of_two(64) - int_mp(1)).truncated_frexp(e) == 1.0 - std::ldexp(1.0, -53));
    CHECK(e == 64);
    CHECK(power_of_two(2000).truncated_frexp(e) == 0.5);
    CHECK(e == 2001);
}

TEST_CASE("int_mp from_double truncates toward zero") {
    CHECK(int_mp::from_double(0.0).zero());
    CHECK(int_mp::from_double(-0.75).zero());
    CHECK(int_mp::from_double(-2.5) == -2);
    CHECK(int_mp::from_double(9007199254740993.0) == 9007199254740992);
    CHECK(int_mp::from_double(std::ldexp(1.0, 100)) == power_of_two(100));
    CHECK(int_mp::from_double(-std::ldexp(3.0, 70)) == -(power_of_two(71) + power_of_two(70)));
    CHECK(int_mp::from_double(std::ldexp(1.0, 63)).str(16) == "8000000000000000");
}

namespace {

// a, b, then a / b truncated and floored, from Python ints. The divisions by
// three limbs take Knuth's add-back step.
const char* const kDivisions[][4] = {
    {"7", "2", "3", "3"},
    {"-7", "2", "-3", "-4"},
    {"7", "-2", "-3", "-4"},
    {"-7", "-2", "3", "3"},
    {"-8000000000000000", "-1", "8000000000000000", "8000000000000000"},
    {"5", "10000000000000000", "0", "0"},
    {"-5", "10000000000000000", "0", "-1"},
    {"100000000000000000000000000000000000000000000003039", "3b9aca07",
     "44b82f98895147f23df9f377d4795fb1436c44a6bd0", "44b82f98895147f23df9f377d4795fb1436c44a6bd0"},
    {"-100000000000000000000000000000000000000000000003039", "3b9aca07",
     "-44b82f98895147f23df9f377d4795fb1436c44a6bd0",
     "-44b82f98895147f23df9f377d4795fb1436c44a6bd1"},
    {"2fffffffffffffffe000000000000000000000000ffffffff", "180000000000000008000000000000001",
     "1fffffffffffffffd", "1fffffffffffffffd"},
    {"-7fffffffffffffff800000000000000000000000000000007fffffffffffffff",
     "7fffffffffffffff8000000000000000534291e584299dd7", "-ffffffffffffffff", "-10000000000000000"},
    {"800000000000000080000000000000018000000000000000",
     "-80000000000000008000000000000001d173efcd694038e6", "0", "-1"},
    {"212f052e2af4170eeca89b3c9cb6981b6af609a4000b810322e0df2d20dac705661479b1b",
     "a6ebd69fe29d76d4330f1446beab0c11fdf", "32e4735bc8fbbcbde5c0994164d8399f767c45",
     "32e4735bc8fbbcbde5c0994164d8399f767c45"},
    {"-212f052e2af4170eeca89b3c9cb6981b6af609a4000b810322e0df2d20dac705661479b1b",
     "a6ebd69fe29d76d4330f1446beab0c11fdf", "-32e4735bc8fbbcbde5c0994164d8399f767c45",
     "-32e4735bc8fbbcbde5c0994164d8399f767c45"},
    {"212f052e2af4170eeca89b3c9cb6981b6af609a4000b810322e0df2d20dac705661479b1c",
     "-a6ebd69fe29d76d4330f1446beab0c11fdf", "-32e4735bc8fbbcbde5c0994164d8399f767c45",
     "-32e4735bc8fbbcbde5c0994164d8399f767c46"},
};

}  // namespace

TEST_CASE("int_mp division truncates and floors as Python ints") {
    for (const auto& d : kDivisions) {
        INFO(d[0], " / ", d[1]);
        const int_mp a = hex(d[0]), b = hex(d[1]);
        CHECK((a / b).str(16) == d[2]);
        CHECK(floor_div(a, b).str(16) == d[3]);
    }
}

TEST_CASE("int_mp compares by value") {
    std::vector<int_mp> values;  // increasing
    for (int bits : {130, 64}) values.push_back(-power_of_two(bits));
    for (int v : {-1, 0, 1}) values.push_back(int_mp(v));
    for (int bits : {63, 64, 130}) values.push_back(power_of_two(bits));
    for (size_t i = 0; i < values.size(); ++i)
        for (size_t j = 0; j < values.size(); ++j) {
            INFO(i, " ", j);
            CHECK(compare(values[i], values[j]) == (i > j) - (i < j));
            CHECK((values[i] < values[j]) == (i < j));
            CHECK((values[i] > values[j]) == (i > j));
        }
}

}  // namespace quicksynth::test

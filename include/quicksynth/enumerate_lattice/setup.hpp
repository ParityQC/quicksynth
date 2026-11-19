// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// The set-up of a search direction, all a walk needs besides its scale, and the
// cache that builds each one once: the region's Loewner ellipsoid around the
// direction, LLL-reduced, and the triangular factors of the walk's forms.
#pragma once

#include "quicksynth/enumerate_lattice/enumerator.hpp"
#include "quicksynth/enumerate_lattice/lll.hpp"
#include "quicksynth/enumerate_lattice/search_region.hpp"
#include "quicksynth/data_structures/gate_set.hpp"
#include "quicksynth/data_structures/su2.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace quicksynth {

// A search direction's exact bits, the key of the caches.
using DirectionBits = std::array<uint64_t, 4>;
using SubBits = std::array<uint64_t, 5>;  // direction and sublattice

struct DirectionBitsHash {
    template <size_t N>
    size_t operator()(const std::array<uint64_t, N>& v) const {
        uint64_t h = 0xcbf29ce484222325ULL;
        for (uint64_t x : v) {
            h ^= x;
            h *= 0x100000001b3ULL;
            h ^= h >> 29;
        }
        return static_cast<size_t>(h);
    }
};

// A map with a bound on its size: an insertion that finds it holding more than
// `capacity` entries clears it first. Unlike evicting single entries, this makes
// what it holds depend on the order of the insertions alone.
//
// Insertions come in batches, planned key by key and committed at the end. That
// leaves the map as inserting the keys one at a time would: an entry from before
// the batch counts until the batch's first clear, and of the batch's own keys
// those inserted after its last clear are kept. Until commit() the map itself is
// untouched.
template <class Key, class Value, class Hash>
class ClearingMap {
  public:
    explicit ClearingMap(int64_t capacity) : capacity_(capacity) {}

    class Batch {
      public:
        explicit Batch(ClearingMap& map) : map_(map), size_(map.entries_.size()) {}

        // The value stored before the batch, unless the batch has cleared it.
        const Value* stored(const Key& key) const {
            if (cleared_) return nullptr;
            const auto it = map_.entries_.find(key);
            return it == map_.entries_.end() ? nullptr : &it->second;
        }
        // Whether the batch inserted `key` after its last clear.
        bool inserted(const Key& key) const { return kept_.count(key) > 0; }
        // Inserts a key neither stored() nor inserted(), clearing first if full.
        void insert(const Key& key) {
            if (static_cast<int64_t>(size_) > map_.capacity_) {
                cleared_ = true;
                size_ = 0;
                kept_.clear();
            }
            ++size_;
            kept_.insert(key);
        }
        // Writes the batch into the map: value(key) for every key kept.
        template <class ValueOf>
        void commit(const ValueOf& value) {
            if (cleared_) map_.entries_.clear();
            for (const Key& key : kept_) map_.entries_.emplace(key, value(key));
        }

      private:
        ClearingMap& map_;
        size_t size_;  // the map's size as the insertions so far leave it
        bool cleared_ = false;
        std::unordered_set<Key, Hash> kept_;
    };

    // A new batch on this map, kept as a local of the function that plans it.
    // The plan lives in the batch alone, so a batch never committed (say because
    // computing the values threw) disappears with that function, and the map
    // needs no cleanup.
    Batch batch() { return Batch(*this); }

  private:
    int64_t capacity_;
    std::unordered_map<Key, Value, Hash> entries_;
};

// An LLL-reduced basis of a direction's lattice at one level, the sublattice delta^sub
// Z[zeta]^2: the n x n matrix B whose column j is basis vector j in power-basis
// coefficients, with to_ball B LLL-reduced for the direction's ellipsoid, and B^-1.
struct SublatticeBasis {
    int sub = 0;
    // B: a point's coordinates in this basis to its power-basis coefficients, the walk's
    // transform; its columns are points of delta^sub Z[zeta]^2
    std::vector<int64_t> to_power_basis;
    // B^-1: power-basis coefficients to coordinates in this basis, as for the walk's centre.
    // Its entries k / 2^sub are exact doubles while |k| < 2^53 (next_sublattice_basis()
    // checks).
    std::vector<double> from_power_basis;
};

// Set-up of one direction, built at its reference scale. Valid at every scale
// of the same parity, where the basis differs by a power of two. Holds all a
// walk needs besides its scale (walk_of()). The matrices are the upper
// triangular factors R of the walk's forms at the reference scale, n x n and
// zero past their rank; the conjugate forms keep their Q for the offsets.
struct Setup {
    int reference_scale = 0;           // K0, the scale it was built at
    Vec4 direction{};                  // the target direction v
    int lattice_dimension = 0;         // n
    double ellipsoid_radius_sq = 0.0;  // P, of the walk's ellipsoid
    double conjugate_radius_sq = 0.0;  // P - 1, of the conjugate block
    std::shared_ptr<const SublatticeBasis> basis;
    std::vector<double> ellipsoid_upper;      // the walk's ellipsoid
    std::vector<double> reduced_center;       // n, the ellipsoid's centre in this basis
    std::vector<double> cap_upper;            // the cap block, rows 0..3
    QrFactors conjugate_block;                // R rows 0..n-5
    std::vector<QrFactors> conjugate_places;  // per place where several, R rows 0..3
};

// The basis of the next sublattice, delta^(sub+1) Z[zeta]^2, from `parent`, the basis of
// delta^sub Z[zeta]^2 for the same direction: parent's columns combined into ones divisible
// by delta once more (parent B Kernel), then LLL-reduced against to_ball, the direction's
// ellipsoid's (LoewnerEllipsoid). Empty when an entry would pass the walk's limits.
std::optional<SublatticeBasis> next_sublattice_basis(const GateSet& gateset,
                                                     const SublatticeBasis& parent,
                                                     const std::vector<double>& to_ball);

// The set-ups of the directions met so far, their sublattice bases and the
// embeddings of the scales, so that a later scale or problem in the same
// direction pays only for its walk. A basis serves every scale, a set-up
// every scale of its parity. Set-ups and bases are bounded by `capacity`
// each (ClearingMap), so a batch leaves them as its problems taken one at a time
// would.
class SetupCache {
  public:
    explicit SetupCache(int64_t capacity) : setups_(capacity), bases_(capacity) {}

    // The set-up of each problem at the scale K = `denominator_exponent`, built
    // where missing.
    std::vector<std::shared_ptr<const Setup>> for_problems(
        const SearchRegion& region, int threads, const std::vector<Enumerator::Problem>& problems,
        int denominator_exponent, int sub);
    // Embedding at the scale K = `denominator_exponent`.
    const std::vector<double>& embedding(const GateSet& gateset, int denominator_exponent);

  private:
    // What a set-up is valid for: one direction and determinant class, every
    // scale of one parity.
    struct Key {
        SubBits direction;
        int det_class;
        int scale_parity;
        bool operator==(const Key& o) const {
            return direction == o.direction && det_class == o.det_class &&
                   scale_parity == o.scale_parity;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            return DirectionBitsHash()(k.direction) ^
                   (static_cast<size_t>(k.det_class) * 0x9e3779b97f4a7c15ULL) ^
                   (static_cast<size_t>(k.scale_parity) << 1);
        }
    };

    // The set-ups of `directions`, none of them cached, at the scale K =
    // `denominator_exponent`.
    std::vector<std::shared_ptr<const Setup>> build(const SearchRegion& region, int threads,
                                                    const std::vector<Vec4>& directions,
                                                    int denominator_exponent, int sub);

    ClearingMap<Key, std::shared_ptr<const Setup>, KeyHash> setups_;
    ClearingMap<SubBits, std::shared_ptr<const SublatticeBasis>, DirectionBitsHash> bases_;
    std::map<int, std::vector<double>> embeddings_;  // by scale
};

}  // namespace quicksynth

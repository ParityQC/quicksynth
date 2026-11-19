// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// LRU cache bounded by count or weight. Not thread-safe; the owner locks.
#pragma once

#include <cstdint>
#include <list>
#include <utility>

namespace quicksynth {

// Keeps the most recent entries with total weight <= `capacity` (default 1 each).
template <class Key, class Value>
class Lru {
  public:
    explicit Lru(int64_t capacity) : capacity_(capacity) {}

    Value* find(const Key& key) {
        for (auto it = items_.begin(); it != items_.end(); ++it)
            if (it->key == key) {
                items_.splice(items_.begin(), items_, it);
                return &items_.front().value;
            }
        return nullptr;
    }

    // Evicts the oldest entries until `weight` more fits or the cache is empty.
    void make_room(int64_t weight) {
        while (!items_.empty() && weight_ + weight > capacity_) {
            weight_ -= items_.back().weight;
            items_.pop_back();
        }
    }

    // Inserts `value` as most recent, after make_room(weight).
    Value& put(Key key, Value value, int64_t weight = 1) {
        make_room(weight);
        items_.push_front({std::move(key), std::move(value), weight});
        weight_ += weight;
        return items_.front().value;
    }

  private:
    struct Entry {
        Key key;
        Value value;
        int64_t weight;
    };
    int64_t capacity_;
    int64_t weight_ = 0;
    std::list<Entry> items_;
};

}  // namespace quicksynth

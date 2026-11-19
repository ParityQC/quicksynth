// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

#include "quicksynth/policy.hpp"

#include "quicksynth/search_scales/scan.hpp"
#include "quicksynth/support/host.hpp"

#include <optional>
#include <string>

namespace quicksynth {

std::vector<Stage> stages(double epsilon, const GateSet& gateset,
                          const ScanThenSplit::Settings& first) {
    const int base = prefix_scale_for(epsilon, gateset);
    const int top = gateset.step != 1 && base >= 4 ? base + 1 : 4;
    ScanThenSplit::Settings later = first;
    later.scan.stopping_rule = StoppingRule::fixed;
    later.scan.point_budget = kPointBudget;
    later.max_prefixes = std::nullopt;
    std::vector<Stage> out{{first, "{}"}};
    // a change only if stage 1 caps the prefixes
    const std::string uncapped = first.max_prefixes ? "'max_prefixes': None, " : "";
    for (int k = base + 1; k <= top; ++k) {
        Stage stage{later, "{" + uncapped + "'prefix_scale': " + std::to_string(k) +
                               ", 'split_budget': " + float_repr(3e8) + "}"};
        stage.settings.prefix_scale = k;
        stage.settings.split_budget = 3e8;
        out.push_back(stage);
    }
    // Clifford+T below kDeeperEpsilon: the last stage splits one prefix scale
    // past stage 2, at kDeeperLastBudget
    const bool deeper = gateset.step != 1 && base > 4;
    const double last_budget = deeper ? kDeeperLastBudget : 3e9;
    Stage last{later, "{" + uncapped +
                          (deeper ? "'prefix_scale': " + std::to_string(top + 1) + ", " : "") +
                          "'split_budget': " + float_repr(last_budget) +
                          ", 'point_budget': " + float_repr(1e9) + "}"};
    if (deeper) last.settings.prefix_scale = top + 1;
    last.settings.split_budget = last_budget;
    last.settings.scan.point_budget = 1e9;
    out.push_back(last);
    return out;
}

}  // namespace quicksynth

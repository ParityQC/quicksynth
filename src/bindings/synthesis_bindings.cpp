// ===----------------------------------------------------------------------===//
//
// Part of QuickSynth, under the MIT License.
// See <repo-root>/LICENSE.txt for license information.
// Copyright (c) 2026 Parity Quantum Computing GmbH
// SPDX-License-Identifier: MIT
//
// ===----------------------------------------------------------------------===//

// Bindings for synthesise.hpp and policy.hpp. Long calls release the GIL.
#include "quicksynth/bindings/python_conversions.hpp"

#include "quicksynth/policy.hpp"
#include "quicksynth/synthesise.hpp"
#include "quicksynth/synthesise_rotation/rotation_search.hpp"

#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <optional>
#include <string>
#include <utility>

namespace quicksynth::bindings {

namespace {

// One scan, split or rotation walk, keyed as in quicksynth.SynthesisResult.log.
nb::dict log_entry(int stage, const SearchLog& e) {
    nb::dict d;
    d["search"] = e.rotation ? "rotation" : e.prefix_scale ? "split" : "scan";
    d["scale"] = e.scale;
    if (e.sub) d["sub"] = e.sub;
    if (e.prefix_scale)
        d["prefix_scale"] = *e.prefix_scale;
    else
        d["parity"] = e.det_class;
    d["seconds"] = e.seconds;
    d["estimated_points"] = e.estimated;
    d["points"] = e.points;
    d["best_cost"] = e.best_cost;
    d["stage"] = stage;
    if (e.rotation) {
        d["candidates"] = e.rotation->candidates;
        d["solvable"] = e.rotation->solvable;
        d["undecided"] = e.rotation->undecided;
        d["capped"] = e.rotation->capped;
        d["pricings"] = e.rotation->pricings;
    }
    return d;
}

nb::object summary_dict(const std::optional<SearchSummary>& summary) {
    if (!summary) return nb::none();
    const SearchSummary& s = *summary;
    nb::dict d;
    d["scales"] = s.scales;
    d["first_hit_scale"] = s.first_hit_scale;
    d["last_full_scale"] = s.last_full_scale;
    d["lower_bound"] = s.lower_bound;
    d["estimated_points"] = s.estimated_points;
    d["scan_cost"] = s.scan_cost;
    d["scan_certified"] = s.scan_certified;
    d["split_ran"] = s.split_ran;
    d["source"] = s.source == Source::split        ? "split"
                  : s.source == Source::rotation   ? "rotation"
                  : s.source == Source::clifford_t ? "clifford_t"
                                                   : "scan";
    return d;
}

// QuickSynth on a (2, 2) unitary, GIL released.
Synthesised synthesise(const QuickSynth& q, const Matrix& target, double epsilon) {
    const Unitary U = unitary_of(target);
    nb::gil_scoped_release release;
    return q.synthesise(U, epsilon);
}

void prepare(const QuickSynth& q, double epsilon) {
    nb::gil_scoped_release release;
    q.prepare(epsilon);
}

void bind_results(nb::module_& mod) {
    nb::class_<Synthesised>(mod, "Synthesised", "A checked word (synthesise.hpp).")
        .def_ro("word", &Synthesised::word)
        .def_prop_ro("cost", [](const Synthesised& r) { return r.cost.value; })
        .def_prop_ro("n_t", [](const Synthesised& r) { return r.cost.n_t; })
        .def_prop_ro("n_q", [](const Synthesised& r) { return r.cost.n_q; })
        .def_ro("distance", &Synthesised::distance)
        .def_ro("certified", &Synthesised::certified)
        .def_ro("stages", &Synthesised::stages)
        .def_prop_ro(
            "log",
            [](const Synthesised& r) {
                nb::list out;
                for (const auto& [stage, entry] : r.log) out.append(log_entry(stage, entry));
                return out;
            },
            "QuickSynth's walks, as SynthesisResult.log documents.")
        .def_prop_ro(
            "search", [](const Synthesised& r) { return summary_dict(r.search); },
            "QuickSynth's answering stage, as SynthesisResult.search documents.");
}

void bind_entry_points(nb::module_& mod) {
    nb::class_<QuickSynth>(mod, "QuickSynth", "QuickSynth's one-call synthesis.")
        .def(
            "__init__",
            [](QuickSynth* self, const std::string& gateset, int64_t t_cost, int64_t q_cost,
               std::optional<int> extra_scales, std::optional<double> point_budget, int workers,
               std::optional<double> time_limit, std::string time_limit_text, bool verbose,
               const std::string& stopping_rule, bool fast, std::optional<double> rate,
               bool partial_walk, std::optional<bool> fine_levels, bool level_search,
               bool skip_setup_levels, std::optional<int> extra_levels,
               std::optional<int> split_levels, bool rotation_path, double rotation_tolerance,
               std::optional<int> rotation_solvable) {
                QuickSynth::Settings s;
                s.gateset = GateSet::from_name(gateset);
                s.price = {t_cost, q_cost};
                s.extra_scales = extra_scales;
                s.point_budget = point_budget;
                s.workers = workers;
                s.time_limit = time_limit;
                s.time_limit_text = std::move(time_limit_text);
                s.verbose = verbose;
                s.stopping_rule = stopping_rule_named(stopping_rule);
                s.fast = fast;
                s.rate = rate;
                s.partial_walk = partial_walk;
                s.fine_levels = fine_levels;
                s.level_search = level_search;
                s.skip_setup_levels = skip_setup_levels;
                s.extra_levels = extra_levels;
                s.split_levels = split_levels;
                s.rotation_path = rotation_path;
                s.rotation_tolerance = rotation_tolerance;
                s.rotation_solvable = rotation_solvable;
                new (self) QuickSynth(std::move(s));
            },
            nb::arg("gateset"), nb::arg("t_cost"), nb::arg("q_cost"),
            nb::arg("extra_scales").none() = nb::none(),
            nb::arg("point_budget").none() = nb::none(), nb::arg("workers") = 1,
            nb::arg("time_limit").none() = nb::none(), nb::arg("time_limit_text") = "",
            nb::arg("verbose") = false, nb::arg("stopping_rule") = "fixed", nb::arg("fast") = false,
            nb::arg("rate").none() = nb::none(), nb::arg("partial_walk") = false,
            nb::arg("fine_levels").none() = nb::none(), nb::arg("level_search") = true,
            nb::arg("skip_setup_levels") = true, nb::arg("extra_levels").none() = nb::none(),
            nb::arg("split_levels").none() = nb::none(), nb::arg("rotation_path") = true,
            nb::arg("rotation_tolerance") = 1e-2, nb::arg("rotation_solvable").none() = nb::none(),
            "time_limit_text is the limit as the TimeoutError's message gives it.")
        .def("synthesise", &synthesise, nb::arg("target"), nb::arg("epsilon"),
             "The cheapest word found within epsilon of a (2, 2) unitary.")
        .def("prepare", &prepare, nb::arg("epsilon"),
             "Build every stage's prefix ball and collect its prefix sets.");
}

}  // namespace

void bind_synthesis(nb::module_& mod) {
    // SynthesisTimeout -> builtin TimeoutError
    nb::register_exception_translator([](const std::exception_ptr& p, void*) {
        try {
            std::rethrow_exception(p);
        } catch (const SynthesisTimeout& e) {
            PyErr_SetString(PyExc_TimeoutError, e.what());
        }
    });
    mod.attr("EPSILON_MARGIN") = kEpsilonMargin;
    mod.attr("SPLIT_BELOW") = kSplitBelow;
    mod.attr("RATE") = kDefaultRate;
    mod.attr("FAST_RATE") = kFastRate;
    mod.attr("FINE_FAST_RATE") = kFineFastRate;
    mod.attr("FINE_RATE") = kFineDefaultRate;
    mod.attr("PREFIX_RADIUS") = kPrefixRadius;
    mod.attr("PREFIX_COUNT") = kMaxPrefixes;
    mod.attr("ROTATION_SOLVABLE") = kRotationSolvable;
    mod.attr("ROTATION_SOLVABLE_FAST") = kRotationSolvableFast;
    mod.def(
        "reference_throughput",
        [](const std::string& gateset) {
            return reference_throughput(GateSet::from_name(gateset));
        },
        nb::arg("gateset"), "The fixed stopping rule's reference throughput.");
    bind_results(mod);
    bind_entry_points(mod);
}

}  // namespace quicksynth::bindings

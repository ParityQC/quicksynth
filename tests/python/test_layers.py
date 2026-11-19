#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Checks the C++ pipeline layout (docs/architecture.md, Layout).

One directory per pipeline step under include/quicksynth/ and src/. Files
include only their own and later steps; only the bindings include Python. The
compiler cannot enforce this. C++ tests mirror the steps and CMake builds every
source and test.
"""

from __future__ import annotations

import re
from pathlib import Path

from conftest import ROOT

INCLUDE = ROOT / "include"
SRC = ROOT / "src"

# Bindings, entry (""), the steps in call order, then the foundations. A file
# may include its own directory and any later one.
ORDER = (
    "bindings",
    "",
    "synthesise_rotation",
    "solve_norm_equation",
    "search_scales",
    "enumerate_lattice",
    "price_candidates",
    "exact_synthesis",
    "data_structures",
    "support",
)

_INCLUDE = re.compile(r'^[ \t]*#[ \t]*include[ \t]*([<"])([^>"]+)[>"]', re.MULTILINE)
_PYTHON = re.compile(r"(nanobind|pybind11)/|Python\.h$")

Files = dict[str, tuple[list[str], list[str]]]


def _files() -> Files:
    """Maps each file of include/ and src/ to its quoted and angle includes.

    Headers are keyed by include path (quicksynth/...), sources by repository
    path (src/...).
    """
    files: Files = {}
    for root, base in ((INCLUDE, INCLUDE), (SRC, ROOT)):
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            quoted: list[str] = []
            system: list[str] = []
            if path.suffix in {".hpp", ".cpp"}:
                for kind, name in _INCLUDE.findall(path.read_text()):
                    (quoted if kind == '"' else system).append(name)
            files[path.relative_to(base).as_posix()] = (quoted, system)
    return files


def _step(name: str) -> str:
    """Step of quicksynth/[<step>/]<file>.hpp or src/[<step>/]<file>.cpp."""
    parts = name.split("/")
    return parts[1] if len(parts) == 3 else ""


def _resolve(files: Files, name: str, include: str) -> str | None:
    """Resolves a quoted include of `name` by include path."""
    if include.startswith("quicksynth/") and include in files:
        return include
    return None


def test_every_file_sits_in_a_step() -> None:
    files = _files()
    headers = [n for n in files if n.startswith("quicksynth/")]
    sources = [n for n in files if n.startswith("src/")]
    assert headers and sources
    assert len(headers) + len(sources) == len(files), "include/ holds only quicksynth/"
    misplaced = [
        n
        for n in files
        if not (
            (n in headers and n.endswith(".hpp") and n.count("/") <= 2)
            or (n in sources and n.endswith(".cpp") and n.count("/") <= 2)
        )
        or _step(n) not in ORDER
    ]
    assert misplaced == []
    # every step has headers and sources
    assert {_step(n) for n in headers} == set(ORDER)
    assert {_step(n) for n in sources} == set(ORDER)


def test_each_source_is_in_the_step_of_its_header() -> None:
    files = _files()
    steps: dict[str, set[str]] = {}
    for name in files:
        steps.setdefault(Path(name).stem, set()).add(_step(name))
    assert {stem: s for stem, s in steps.items() if len(s) > 1} == {}


def test_includes_name_headers_by_their_include_path() -> None:
    files = _files()
    unknown = [
        (name, inc)
        for name, (quoted, _) in files.items()
        for inc in quoted
        if _resolve(files, name, inc) is None
    ]
    assert unknown == []


def test_includes_follow_the_pipeline() -> None:
    """Files include only their own and later steps; bindings include any."""
    files = _files()
    backward = {
        (name, inc)
        for name, (quoted, _) in files.items()
        for inc in quoted
        if inc in files and ORDER.index(_step(inc)) < ORDER.index(_step(name))
    }
    assert backward == set(), "includes against the pipeline order"


def test_only_the_bindings_include_python() -> None:
    files = _files()
    python = [
        name
        for name, (_, system) in files.items()
        if _step(name) != "bindings" and any(_PYTHON.match(s) for s in system)
    ]
    assert python == []


def test_cmake_builds_every_source() -> None:
    """Sources are listed by path from src/."""
    text = (ROOT / "CMakeLists.txt").read_text()
    code = "\n".join(line.split("#", 1)[0] for line in text.splitlines())
    listed = set(re.findall(r"[\w/-]+\.cpp", code))
    sources = {n.removeprefix("src/") for n in _files() if n.endswith(".cpp")}
    assert sources == listed


CPP_TESTS = ROOT / "tests" / "cpp"


def test_the_cpp_tests_mirror_the_steps() -> None:
    """tests/cpp/[<step>/]test_<x>.cpp tests include/quicksynth/[<step>/]<x>.hpp."""
    tests = sorted(CPP_TESTS.rglob("test_*.cpp"))
    assert tests
    misplaced = [
        path.relative_to(CPP_TESTS).as_posix()
        for path in tests
        if not (
            INCLUDE
            / "quicksynth"
            / path.parent.relative_to(CPP_TESTS)
            / f"{path.stem[len('test_') :]}.hpp"
        ).exists()
    ]
    assert misplaced == []
    directories = {
        path.parent.relative_to(CPP_TESTS).as_posix().removeprefix(".")
        for path in CPP_TESTS.rglob("*.cpp")
    }
    assert directories <= set(ORDER)


def test_every_cpp_test_is_built() -> None:
    text = (CPP_TESTS / "CMakeLists.txt").read_text()
    code = "\n".join(line.split("#", 1)[0] for line in text.splitlines())
    listed = set(re.findall(r"[\w/-]+\.cpp", code))
    sources = {p.relative_to(CPP_TESTS).as_posix() for p in CPP_TESTS.rglob("*.cpp")}
    assert sources == listed


def test_each_cpp_test_case_is_its_own_ctest_filter() -> None:
    """Test case names have no commas: doctest splits --test-case filters at them,
    so such a test would select nothing and pass."""
    names = [
        name
        for path in CPP_TESTS.rglob("*.cpp")
        for name in re.findall(r'TEST_CASE\("([^"]*)"', path.read_text())
    ]
    assert names
    assert [n for n in names if "," in n or "\\" in n] == []
    assert len(set(names)) == len(names)

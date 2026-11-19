#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

cmake_minimum_required(VERSION 3.18)
project(quicksynth LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT MSVC)
  add_compile_options(-O3)
endif()

# Windows: the core's arithmetic is written in __int128, which MSVC's compiler
# lacks and clang-cl, Visual Studio's LLVM toolset (-T ClangCL), has. Its
# conversions between __int128 and double are calls into compiler-rt's
# builtins, which the Microsoft linker does not find by itself.
if(MSVC)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
      "The C++ core needs __int128: on Windows build it with clang-cl "
      "(Visual Studio generator: -T ClangCL; Ninja: CMAKE_CXX_COMPILER=clang-cl).")
  endif()
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} /clang:--rtlib=compiler-rt /clang:-print-libgcc-file-name
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE CLANG_RT_BUILTINS)
  # An LLVM may name a library it does not have, or keep it under the other of
  # its two layouts: then take the newest one in the LLVM installation itself,
  # <llvm>/lib/clang/<version>/lib/{windows,x86_64-pc-windows-msvc}.
  if(NOT EXISTS "${CLANG_RT_BUILTINS}")
    get_filename_component(_LLVM_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_LLVM_ROOT "${_LLVM_BIN}" DIRECTORY)
    file(GLOB _CLANG_RT_FOUND
      "${_LLVM_ROOT}/lib/clang/*/lib/windows/clang_rt.builtins-x86_64.lib"
      "${_LLVM_ROOT}/lib/clang/*/lib/x86_64-pc-windows-msvc/clang_rt.builtins.lib")
    if(_CLANG_RT_FOUND)
      list(SORT _CLANG_RT_FOUND COMPARE NATURAL ORDER DESCENDING)
      message(STATUS "clang-cl names '${CLANG_RT_BUILTINS}', which does not exist")
      list(GET _CLANG_RT_FOUND 0 CLANG_RT_BUILTINS)
    endif()
  endif()
  if(NOT EXISTS "${CLANG_RT_BUILTINS}")
    message(FATAL_ERROR "clang-cl names no compiler-rt builtins library ('${CLANG_RT_BUILTINS}'), "
      "and ${_LLVM_ROOT}/lib/clang/*/lib holds none")
  endif()
  message(STATUS "compiler-rt builtins: ${CLANG_RT_BUILTINS}")
  add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
endif()
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Python + nanobind (scikit-build-core passes the interpreter; a manual build
# finds it through the venv).
find_package(Python 3.12 COMPONENTS Interpreter Development.Module REQUIRED)
if(NOT nanobind_ROOT)
  execute_process(
    COMMAND ${Python_EXECUTABLE} -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
endif()
find_package(nanobind CONFIG REQUIRED)

# The native library proper: plain C++, no Python, unit-testable on its own.
# The headers are in include/quicksynth/, the sources in src/, both with one
# directory per step of a call: one source list per step, in the order a call
# runs them, then the foundations. A file includes only its own step, the
# steps after it and the foundations (docs/architecture.md);
# tests/python/test_layers.py checks this, since every header is on the one
# include path.
set(_QUICKSYNTH_SYNTHESISE policy.cpp synthesise.cpp)
set(_QUICKSYNTH_SYNTHESISE_ROTATION
  synthesise_rotation/diagonal_form.cpp synthesise_rotation/rotation_body.cpp
  synthesise_rotation/rotation_search.cpp)
set(_QUICKSYNTH_SOLVE_NORM_EQUATION
  solve_norm_equation/integer_factor.cpp solve_norm_equation/prime_field.cpp
  solve_norm_equation/prime_ideals.cpp solve_norm_equation/relative_norm.cpp)
set(_QUICKSYNTH_SEARCH_SCALES
  search_scales/level_search.cpp search_scales/scale_bound.cpp search_scales/scan.cpp
  search_scales/scan_then_split.cpp search_scales/split.cpp search_scales/stopping_rule.cpp)
set(_QUICKSYNTH_ENUMERATE_LATTICE
  enumerate_lattice/enumerator.cpp enumerate_lattice/lll.cpp enumerate_lattice/search_region.cpp
  enumerate_lattice/setup.cpp enumerate_lattice/walk.cpp)
set(_QUICKSYNTH_PRICE_CANDIDATES price_candidates/price.cpp)
set(_QUICKSYNTH_EXACT_SYNTHESIS exact_synthesis/descent.cpp exact_synthesis/normal_form.cpp)
set(_QUICKSYNTH_DATA_STRUCTURES
  data_structures/bloch.cpp data_structures/clifford.cpp data_structures/cost_ball.cpp
  data_structures/cyclotomic_matrix.cpp data_structures/distance.cpp
  data_structures/double_double.cpp data_structures/gate_set.cpp
  data_structures/int_mp.cpp data_structures/number_field.cpp data_structures/prefix_sets.cpp
  data_structures/su2.cpp)
set(_QUICKSYNTH_SUPPORT support/diagnostics.cpp support/host.cpp support/parallel.cpp)
set(_QUICKSYNTH_SOURCES
  ${_QUICKSYNTH_SYNTHESISE} ${_QUICKSYNTH_SYNTHESISE_ROTATION} ${_QUICKSYNTH_SOLVE_NORM_EQUATION}
  ${_QUICKSYNTH_SEARCH_SCALES}
  ${_QUICKSYNTH_ENUMERATE_LATTICE}
  ${_QUICKSYNTH_PRICE_CANDIDATES} ${_QUICKSYNTH_EXACT_SYNTHESIS} ${_QUICKSYNTH_DATA_STRUCTURES}
  ${_QUICKSYNTH_SUPPORT})
list(TRANSFORM _QUICKSYNTH_SOURCES PREPEND src/)
add_library(quicksynth_core STATIC ${_QUICKSYNTH_SOURCES})
find_package(Threads REQUIRED)
target_link_libraries(quicksynth_core PUBLIC Threads::Threads)
if(MSVC)
  target_link_libraries(quicksynth_core PRIVATE ${CLANG_RT_BUILTINS})
endif()
target_include_directories(quicksynth_core PUBLIC include)
set_target_properties(quicksynth_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# The Python module quicksynth._core, the step before the entry: its header is
# in include/quicksynth/bindings/, its sources in src/bindings/.
set(_QUICKSYNTH_BINDINGS
  bindings/module.cpp bindings/python_conversions.cpp bindings/exact_bindings.cpp
  bindings/synthesis_bindings.cpp)
list(TRANSFORM _QUICKSYNTH_BINDINGS PREPEND src/)
nanobind_add_module(_core NB_STATIC ${_QUICKSYNTH_BINDINGS})
target_link_libraries(_core PRIVATE quicksynth_core)

# The lattice inputs (search_region.hpp, su2.hpp) and the LLL's Gram-Schmidt (lll.cpp)
# are fixed sequences of correctly rounded operations; a fused multiply-add
# would change them. The
# module compiles the library's inline functions too and the linker keeps one
# copy of each, so both targets need the same rounding. clang-cl takes the
# option through /clang:, and its /W4 is clang's -Wall -Wextra.
if(MSVC)
  set(QUICKSYNTH_CXX_OPTIONS /clang:-ffp-contract=off /W4)
else()
  set(QUICKSYNTH_CXX_OPTIONS -ffp-contract=off -Wall -Wextra)
endif()
foreach(_target quicksynth_core _core)
  target_compile_options(${_target} PRIVATE ${QUICKSYNTH_CXX_OPTIONS})
endforeach()

# The lattice walk's functions each start on a 64-byte cache line. Unaligned,
# every edit anywhere in walk.cpp moved its hot loops within their cache lines
# and shifted single-threaded timings by 1-2%, whatever the edit did; aligned,
# they stay put, so that a benchmark measures the change itself. It costs
# nothing measurable. Aligning the loops as well (-falign-loops=32) cost 10% on
# Clifford+T. The path is from src/, as everywhere in this file.
if(MSVC)
  set(_QUICKSYNTH_WALK_ALIGNMENT /clang:-falign-functions=64)
else()
  set(_QUICKSYNTH_WALK_ALIGNMENT -falign-functions=64)
endif()
set(_QUICKSYNTH_WALK enumerate_lattice/walk.cpp)
set_source_files_properties(src/${_QUICKSYNTH_WALK} PROPERTIES
  COMPILE_OPTIONS "${_QUICKSYNTH_WALK_ALIGNMENT}")

install(TARGETS _core LIBRARY DESTINATION quicksynth)

# nanobind (BSD-3-Clause), and the tsl::robin_map it includes (MIT), are
# compiled into _core: every wheel carries their notices in its metadata, in
# licenses/bundled, from the nanobind the build used. An installed nanobind
# keeps its LICENSE in its dist-info and robin_map's notice only atop its
# headers; a source checkout of nanobind has both files.
if(SKBUILD)
  get_filename_component(_NB_PACKAGE "${nanobind_DIR}" DIRECTORY)
  get_filename_component(_NB_SITE "${_NB_PACKAGE}" DIRECTORY)
  file(GLOB _NB_LICENSE
    "${_NB_PACKAGE}/LICENSE"
    "${_NB_SITE}/nanobind-${nanobind_VERSION}.dist-info/licenses/LICENSE"
    "${_NB_SITE}/nanobind-${nanobind_VERSION}.dist-info/LICENSE")
  if(NOT _NB_LICENSE)
    message(FATAL_ERROR "nanobind's LICENSE is neither in ${_NB_PACKAGE} nor in its dist-info")
  endif()
  list(GET _NB_LICENSE 0 _NB_LICENSE)
  set(_RM_DIR "${_NB_PACKAGE}/ext/robin_map")
  if(EXISTS "${_RM_DIR}/LICENSE")
    file(READ "${_RM_DIR}/LICENSE" _RM_LICENSE)
  else()
    # the comment that opens robin_map.h, without its " * " margin
    file(READ "${_RM_DIR}/include/tsl/robin_map.h" _RM_HEADER)
    string(FIND "${_RM_HEADER}" "*/" _RM_END)
    string(SUBSTRING "${_RM_HEADER}" 0 ${_RM_END} _RM_LICENSE)
    string(REGEX REPLACE "^/\\*\\*\n" "" _RM_LICENSE "${_RM_LICENSE}")
    string(REGEX REPLACE "(^|\n) \\* ?" "\\1" _RM_LICENSE "${_RM_LICENSE}")
    string(STRIP "${_RM_LICENSE}" _RM_LICENSE)
    string(APPEND _RM_LICENSE "\n")
  endif()
  if(NOT _RM_LICENSE MATCHES "^MIT License\n.*Permission is hereby granted")
    message(FATAL_ERROR "no MIT notice for tsl::robin_map in ${_RM_DIR}")
  endif()
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/licenses/tsl-robin-map/LICENSE" "${_RM_LICENSE}")
  install(FILES "${_NB_LICENSE}" DESTINATION "${SKBUILD_METADATA_DIR}/licenses/bundled/nanobind")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/licenses/tsl-robin-map/LICENSE"
    DESTINATION "${SKBUILD_METADATA_DIR}/licenses/bundled/tsl-robin-map")
endif()

# The C++ tests (tests/cpp), registered with ctest; `make test` builds them in
# a tree of their own. The wheel leaves them out.
option(QUICKSYNTH_TESTS "Build the C++ tests of tests/cpp" OFF)
if(QUICKSYNTH_TESTS)
  enable_testing()
  add_subdirectory(tests/cpp)
endif()

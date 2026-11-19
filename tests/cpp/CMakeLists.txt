#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

# doctest and nlohmann_json, pinned by URL and hash. Only headers and doctest's
# ctest script are used; SOURCE_SUBDIR skips their CMake files, which CMake 4
# rejects.
include(FetchContent)
FetchContent_Declare(doctest
  URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.11.tar.gz
  URL_HASH SHA256=632ed2c05a7f53fa961381497bf8069093f0d6628c5f26286161fbd32a560186
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  SOURCE_SUBDIR headers-only)
FetchContent_Declare(nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
  URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  SOURCE_SUBDIR headers-only)
FetchContent_MakeAvailable(doctest nlohmann_json)

# One file per concept; tests/cpp/<step>/ mirrors include/quicksynth/<step>/.
set(_QUICKSYNTH_TESTS
  main.cpp fixtures.cpp oracles.cpp
  test_policy.cpp test_synthesise.cpp
  synthesise_rotation/test_diagonal_form.cpp synthesise_rotation/test_rotation_body.cpp
  synthesise_rotation/test_rotation_search.cpp
  solve_norm_equation/fp_poly_oracle.cpp
  solve_norm_equation/test_integer_factor.cpp solve_norm_equation/test_prime_field.cpp
  solve_norm_equation/test_prime_ideals.cpp solve_norm_equation/test_relative_norm.cpp
  search_scales/replay.cpp search_scales/test_level_search.cpp search_scales/test_scan.cpp
  search_scales/test_scan_then_split.cpp
  search_scales/test_split.cpp search_scales/test_stopping_rule.cpp
  enumerate_lattice/test_enumerator.cpp enumerate_lattice/test_lll.cpp
  enumerate_lattice/test_search_region.cpp enumerate_lattice/test_setup.cpp
  price_candidates/test_price.cpp
  exact_synthesis/test_descent.cpp exact_synthesis/test_normal_form.cpp
  data_structures/test_bloch.cpp data_structures/test_clifford.cpp
  data_structures/test_cost_ball.cpp data_structures/test_cyclotomic_matrix.cpp
  data_structures/test_distance.cpp data_structures/test_double_double.cpp
  data_structures/test_gate_set.cpp
  data_structures/test_int_mp.cpp data_structures/test_number_field.cpp
  data_structures/test_prefix_sets.cpp data_structures/test_su2.cpp)
add_executable(quicksynth_tests ${_QUICKSYNTH_TESTS})
target_link_libraries(quicksynth_tests PRIVATE quicksynth_core)
target_include_directories(quicksynth_tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${doctest_SOURCE_DIR} ${nlohmann_json_SOURCE_DIR}/include)
# Fixtures are read in place.
target_compile_definitions(quicksynth_tests PRIVATE
  QUICKSYNTH_REPOSITORY="${PROJECT_SOURCE_DIR}")
# Same rounding options as the library: the tests compile its inline functions
# and the linker keeps one copy of each.
target_compile_options(quicksynth_tests PRIVATE ${QUICKSYNTH_CXX_OPTIONS})

# One ctest test per test case, labelled by suite (`ctest -LE slow`).
include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)
doctest_discover_tests(quicksynth_tests ADD_LABELS 1)

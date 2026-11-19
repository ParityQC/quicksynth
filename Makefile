# Copyright (c) 2026 Parity Quantum Computing GmbH
# SPDX-License-Identifier: MIT

.DEFAULT_GOAL := help
UV ?= uv

.PHONY: help sync test test-fast test-cpp test-sanitizers lint format typecheck check demo clean cpp wheel devcontainer

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-16s\033[0m %s\n", $$1, $$2}'

sync: ## Install the project and all dev dependencies
	$(UV) sync --all-extras

# The C++ tests build in a tree of their own; CTEST_ARGS="-LE slow" leaves out
# the slow ones, and CMAKE_ARGS reaches the configure step.
CPP_TESTS ?= build/tests
CMAKE_ARGS ?=
CTEST_ARGS ?=
CTEST_JOBS ?= 8
test-cpp: ## Build and run the C++ tests (tests/cpp) with ctest
	$(UV) run cmake -S . -B $(CPP_TESTS) -G Ninja -DQUICKSYNTH_TESTS=ON $(CMAKE_ARGS)
	$(UV) run cmake --build $(CPP_TESTS) --target quicksynth_tests
	$(UV) run ctest --test-dir $(CPP_TESTS) --output-on-failure -j $(CTEST_JOBS) $(CTEST_ARGS)

# UBSan stops at its first finding, so that ctest sees it as a failure.
SANITIZE_FLAGS ?= -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer
test-sanitizers: ## Run the C++ tests without the slow ones under ASan and UBSan
	$(MAKE) test-cpp CPP_TESTS=build/sanitizers CTEST_ARGS="-LE slow" \
		CMAKE_ARGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_FLAGS='$(SANITIZE_FLAGS)'"

PYTEST_ARGS ?=
test: test-cpp ## Run the full test suite, C++ and Python
	$(UV) run pytest $(PYTEST_ARGS)

test-fast: ## Run both test suites without the slow tests
	$(MAKE) test-cpp CTEST_ARGS="-LE slow"
	$(UV) run pytest -m "not slow" -n auto $(PYTEST_ARGS)

lint: ## Check formatting and lint rules
	$(UV) run ruff format --check .
	$(UV) run ruff check .

format: ## Auto-format and auto-fix
	$(UV) run ruff format .
	$(UV) run ruff check --fix .

typecheck: ## Run mypy
	$(UV) run mypy

check: lint typecheck test-fast ## What CI runs: lint, mypy, both suites without the slow tests

cpp: ## Rebuild the C++ core (quicksynth._core) after editing include/ or src/
	$(UV) sync --all-extras --reinstall-package quicksynth

wheel: ## Build sdist and wheel into dist/
	$(UV) build

demo: ## Synthesise a few targets (examples/quickstart.py)
	$(UV) run python examples/quickstart.py

devcontainer: ## Create .devcontainer/local/ from the templates (keeps existing files)
	mkdir -p .devcontainer/local
	@for f in devcontainer.json Dockerfile; do \
		if [ -e .devcontainer/local/$$f ]; then echo "kept .devcontainer/local/$$f"; \
		else cp .devcontainer/local.$$f.example .devcontainer/local/$$f && echo "created .devcontainer/local/$$f"; fi; \
	done

clean: ## Remove caches and build artefacts
	rm -rf build dist .pytest_cache .mypy_cache .ruff_cache htmlcov .coverage
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

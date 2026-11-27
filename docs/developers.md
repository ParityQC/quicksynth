<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Developer guide

How the C++ library and the Python package fit together, how to build and test
them, and how to contribute and release. New contributors start with
[Contributing](../CONTRIBUTING.md).

- [Contributing](../CONTRIBUTING.md): setup, the checks to run before a pull
  request, and CI.
- [Architecture](architecture.md): the steps of a call, the source layout, the
  boundary between Python and C++, and the constants.
- [Development](development.md): prerequisites, building, wheels, the
  documentation site, the tests and the golden fixtures.
- [Fixture records](../tests/fixtures/README.md): how the golden fixtures were
  recorded and every change that re-recorded one.
- [Dev container](../.devcontainer/README.md): the shared and local container
  configs and the per-checkout Python environment.
- [Releasing](releasing.md): the steps of a release to PyPI.

```{toctree}
:maxdepth: 1
:hidden:

architecture
development
Fixture records <../tests/fixtures/README>
Contributing <../CONTRIBUTING>
Dev container <../.devcontainer/README>
releasing
```

<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# Releasing

A release is a tag `vX.Y.Z` on `main`. The tag starts
`.github/workflows/wheels.yml`, which builds the wheels and the sdist and
uploads them to PyPI. Zenodo archives the GitHub release and gives it a DOI.

## Before the first release

1. Remove the classifier `"Private :: Do Not Upload"` from `pyproject.toml`.
   PyPI rejects packages that carry it.
2. On PyPI, add a trusted publisher for the project `quicksynth`: this
   repository, the workflow `wheels.yml` and the environment `pypi`.
3. On GitHub, create the environment `pypi` in the repository settings. A
   required reviewer there makes each upload wait for approval.
4. In `.github/workflows/wheels.yml`, remove `false &&` from the publish job's
   `if:`. Once the Windows wheel build has passed, remove `continue-on-error`
   from the wheels job.
5. In Zenodo's GitHub settings, switch on this repository.
6. In `README.md`, under "From PyPI", replace the two sentences that say
   QuickSynth is not on PyPI yet with "A prebuilt wheel installs QuickSynth
   without a compiler:".

## Each release

1. Set `version` in `pyproject.toml` and in `CITATION.cff` to the same value,
   and add `date-released` to `CITATION.cff`. Run `uv lock` to update the
   package's entry in `uv.lock`. `quicksynth.__version__` reads the version
   from the installed package.
2. Once the paper has an arXiv identifier, add it to `preferred-citation` in
   `CITATION.cff`.
3. Merge these changes to `main` with `ci.yml` green. Then start `wheels.yml`
   by hand (Actions, Wheels, Run workflow) and wait until it passes on all five
   platforms.
4. Tag the release and push the tag:

   ```bash
   git tag -a v1.0.0 -m "QuickSynth 1.0.0"
   git push origin v1.0.0
   ```

   The publish job checks that the tag is `v` followed by the version in
   `pyproject.toml`. It uploads to PyPI only when every build is green.
5. Create the GitHub release from the tag. Zenodo archives it and mints a DOI.
   Add the DOI to `CITATION.cff` (`doi:`) in a follow-up commit.

## Checking the published wheels

- PyPI lists 15 wheels (CPython 3.12, 3.13 and 3.14 on manylinux x86_64 and
  aarch64, macOS arm64 and x86_64, and Windows x64) and the sdist.
- On each available platform, in a new virtual environment, run
  `pip install quicksynth` and then `quicksynth pi/4 1e-2`.
- Each wheel holds the license texts of nanobind and tsl::robin_map, which are
  compiled into the extension module, under
  `quicksynth-<version>.dist-info/licenses/bundled/`.

#  ===----------------------------------------------------------------------===//
#
#  Part of QuickSynth, under the MIT License.
#  See <repo-root>/LICENSE.txt for license information.
#  Copyright (c) 2026 Parity Quantum Computing GmbH
#  SPDX-License-Identifier: MIT
#
#  ===----------------------------------------------------------------------===//

"""Sphinx configuration of the documentation site (``make docs``).

The source directory is the repository root, so the Markdown files link to
each other by their repository paths, as on GitHub; ``include_patterns`` picks
the pages. The API reference imports ``python/quicksynth`` with the compiled
core mocked: the site builds without a C++ compiler.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

project = "QuickSynth"
author = "Parity Quantum Computing GmbH"
copyright = "2026, Parity Quantum Computing GmbH"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_design",
]

root_doc = "docs/index"
include_patterns = [
    "README.md",
    "RESULTS.md",
    "CONTRIBUTING.md",
    ".devcontainer/README.md",
    "docs/**",
]
exclude_patterns = ["docs/_static/**", "docs/_extra/**"]

# GitHub's heading anchors, for links such as development.md#wheels
myst_heading_anchors = 3
myst_enable_extensions = ["colon_fence", "deflist"]

autodoc_mock_imports = ["quicksynth._core"]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True}

html_theme = "pydata_sphinx_theme"
html_title = "QuickSynth"
html_static_path = ["_static"]
html_css_files = ["brand.css", "showcase/showcase.css"]

# GitHub Pages serves index.html at the root; it forwards to docs/index.html
html_extra_path = ["_extra"]
html_show_sourcelink = False
html_theme_options = {
    "logo": {
        "image_light": "_static/logo-light.svg",
        "image_dark": "_static/logo-dark.svg",
        "text": "QuickSynth",
        "alt_text": "QuickSynth",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/ParityQC/quicksynth",
            "icon": "fa-brands fa-github",
        },
    ],
    "navbar_align": "left",
    "show_toc_level": 2,
    "footer_start": ["copyright"],
    "footer_end": [],
}
html_context = {"default_mode": "auto"}


def _orphan(app, docname, source):
    # RESULTS.md: reached from the README's links, outside the navigation
    if docname == "RESULTS":
        source[0] = "---\norphan: true\n---\n" + source[0]


def setup(app):
    app.connect("source-read", _orphan)

# -*- coding: utf-8 -*-
# Sphinx configuration for pyCSAMT v2 documentation.
#
# Build with:
#   cd docs && make html
#
# Requirements (docs/requirements-docs.txt):
#   sphinx>=7.2
#   pydata-sphinx-theme>=0.15
#   numpydoc>=1.7
#   myst-parser>=2.0
#   sphinx-copybutton>=0.5
#   sphinx-design>=0.5
#   sphinx-gallery>=0.16
#   sphinx-autodoc-typehints>=2.0

from __future__ import annotations

import os
import sys

# Keep package imports side-effect-light during documentation builds.
os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
os.environ.setdefault(
    "MPLCONFIGDIR",
    os.path.abspath("../_build/matplotlib"),
)

# -- Path setup ----------------------------------------------------------------
# Make the package importable from the source tree without installation.
sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -------------------------------------------------------
project = "pyCSAMT"
author = "Kouadio K. Laurent"
copyright = "2022-2026, earthai-tech"

try:
    from pycsamt import __version__ as release
    version = ".".join(release.split(".")[:2])   # e.g. "2.0"
except Exception:
    version = release = "2.0"

switcher_version_match = os.environ.get("READTHEDOCS_VERSION", f"v{version.split('.')[0]}")
if switcher_version_match in {"latest", "stable"}:
    switcher_version_match = "v2"

# -- General configuration -----------------------------------------------------
extensions = [
    # Core Sphinx
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.todo",
    # Third-party
    "myst_parser",          # parse .md files (e.g. CHANGELOG, README)
    "sphinx_copybutton",    # copy-button on code blocks
    "sphinx_design",        # grid / tab / card directives
]

templates_path    = ["_templates"]
exclude_patterns  = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix     = {".rst": "restructuredtext", ".md": "markdown"}
master_doc        = "index"

# -- autodoc / autosummary -----------------------------------------------------
autosummary_generate      = True          # stubs not yet written; generate on build
autodoc_default_options   = {
    "members":          True,
    "undoc-members":    False,
    "show-inheritance": True,
    "member-order":     "bysource",
}
autodoc_typehints          = "description"
autodoc_typehints_format   = "short"

# -- Napoleon (NumPy-style docstrings) -----------------------------------------
napoleon_google_docstring       = False
napoleon_numpy_docstring        = True
napoleon_include_init_with_doc  = False
napoleon_use_param              = True
napoleon_use_rtype              = True
napoleon_use_ivar               = True

# -- numpydoc ------------------------------------------------------------------
numpydoc_show_class_members        = False
numpydoc_show_inherited_class_members = False
numpydoc_class_members_toctree     = False
numpydoc_validation_checks         = set()

# -- intersphinx ---------------------------------------------------------------
intersphinx_mapping = {
    "python":     ("https://docs.python.org/3", None),
    "numpy":      ("https://numpy.org/doc/stable", None),
    "scipy":      ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "pandas":     ("https://pandas.pydata.org/docs", None),
    "sklearn":    ("https://scikit-learn.org/stable", None),
}

# -- MathJax -------------------------------------------------------------------
mathjax3_config = {
    "tex": {"tags": "ams"},
}

# -- todo ----------------------------------------------------------------------
todo_include_todos = True          # flip to False before public release

# -- MyST (Markdown) -----------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "tasklist",
]
myst_heading_anchors = 3

# -- HTML output ---------------------------------------------------------------
html_theme = "pydata_sphinx_theme"

html_theme_options = {
    # Logo
    "logo": {
        "image_light": "_static/pycsamt-v2-symbol-logo.svg",
        "image_dark":  "_static/pycsamt-v2-symbol-logo.svg",
        "alt_text":    "pyCSAMT",
    },
    # Navbar
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["version-switcher", "navbar-icon-links", "theme-switcher"],
    "navbar_persistent": ["search-button"],
    "header_links_before_dropdown": 6,
    # Icon links (fill in real URLs when repo is public)
    "icon_links": [
        {
            "name": "GitHub",
            "url":  "https://github.com/earthai-tech/pycsamt",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url":  "https://pypi.org/project/pycsamt/",
            "icon": "fa-brands fa-python",
        },
    ],
    # Sidebar
    "show_toc_level": 2,
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    "primary_sidebar_end": [],
    # Code highlighting
    "pygments_light_style": "tango",
    "pygments_dark_style":  "monokai",
    # Footer
    "footer_start": ["copyright"],
    "footer_end":   ["sphinx-version", "theme-version"],
    # Announcement banner (remove before release)
    "announcement": (
        "pyCSAMT v2 docs are under active construction. "
        "APIs may change before the stable release."
    ),
    # Search
    "search_bar_text": "Search pyCSAMT docs...",
    "switcher": {
        "json_url": "_static/version_switcher.json",
        "version_match": switcher_version_match,
    },
}

html_title            = f"pyCSAMT {version}"
html_short_title      = "pyCSAMT"
html_favicon          = "_static/pycsamt-v2-symbol.ico"
html_static_path      = ["_static"]
html_css_files        = ["custom.css"]
html_show_sourcelink  = True
html_show_sphinx      = False
html_copy_source      = False

# -- LaTeX / PDF ---------------------------------------------------------------
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
}
latex_documents = [
    (
        master_doc,
        "pycsamt.tex",
        "pyCSAMT Documentation",
        author,
        "manual",
    ),
]

# -- Man pages -----------------------------------------------------------------
man_pages = [
    (master_doc, "pycsamt", "pyCSAMT Documentation", [author], 1),
]

# -- Texinfo -------------------------------------------------------------------
texinfo_documents = [
    (
        master_doc,
        "pycsamt",
        "pyCSAMT Documentation",
        author,
        "pycsamt",
        "EM geophysics toolkit for MT, AMT, CSAMT, and CSEM.",
        "Science",
    ),
]

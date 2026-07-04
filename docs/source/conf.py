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
# Repo root, so gallery example scripts can find bundled sample data
# (e.g. data/AMT/WILLY_DATA/) without guessing their own cwd.
os.environ.setdefault("PYCSAMT_DOCS_REPO_ROOT", os.path.abspath("../../"))

# sphinx-gallery executes example scripts at build time; force a headless
# backend before anything else has a chance to pick an interactive one.
import matplotlib
matplotlib.use("Agg")

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
    "sphinx_toggleprompt",  # show/hide >>> prompts on doctest blocks
    "sphinx_gallery.gen_gallery",  # execute example scripts, embed output
]

# -- sphinx-copybutton -----------------------------------------------------
# Strip interpreter prompts when copying code (>>> , ... , $ ).
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

# -- sphinx-toggleprompt ---------------------------------------------------
# Shift the toggle left so it does not overlap the copy button.
toggleprompt_offset_right = 40

# -- sphinx-gallery ---------------------------------------------------------
# Each subpackage that wants runnable, image-producing examples gets its own
# ``<subpackage>/examples`` source dir and a matching ``<subpackage>/auto_examples``
# generated dir. Add new pairs here as more subpackages grow a gallery.
sphinx_gallery_conf = {
    "examples_dirs":  ["emtools/examples"],
    "gallery_dirs":   ["emtools/auto_examples"],
    # Match "plot_*" regardless of path separator — the sphinx-gallery
    # default (r"/plot_") only matches POSIX paths and silently executes
    # zero files on Windows.
    "filename_pattern": r"[\\/]plot_",
    "download_all_examples": False,
    "plot_gallery": True,
    "remove_config_comments": True,
    "backreferences_dir": None,
    # Skip shared helper modules (_datasets.py, _synthetic.py) — they are
    # imported by plot_*.py scripts, not standalone examples.
    "ignore_pattern": r"(__init__|_datasets|_synthetic)\.py",
}

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
    # Square symbol + HTML text: crisper than the wordmark lockup at navbar
    # size, adapts to dark mode, and stays selectable/accessible.
    "logo": {
        "image_light": "_static/logo/pycsamt-v2-symbol.svg",
        "image_dark":  "_static/logo/pycsamt-v2-symbol.svg",
        "text":        "pyCSAMT",
        "alt_text":    "pyCSAMT",
    },
    # Navbar: version switcher sits next to the brand (numpy/pandas style);
    # search + icon links + theme toggle stay inline at the far right so
    # nothing wraps into a second header row.
    "navbar_start": ["navbar-logo", "version-switcher"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links", "theme-switcher"],
    "navbar_persistent": ["search-button"],
    "header_links_before_dropdown": 7,
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
    # Code highlighting — NOTE: pydata-sphinx-theme spells these
    # "pygment_*" (no "s"); the "pygments_*" spelling is ignored and
    # dark mode falls back to black-on-dark unreadable tokens.
    "pygment_light_style": "tango",
    "pygment_dark_style":  "monokai",
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
html_favicon          = "_static/logo/pycsamt-v2-symbol.ico"
html_static_path      = ["_static"]
# Static assets are organised by kind: _static/css/, _static/js/, _static/logo/
# (paths below are relative to the _static output root).
html_css_files        = ["css/custom.css", "css/pycsamt-home.css"]
html_js_files         = [
    ("js/pycsamt-home.js", {"defer": "defer"}),
    ("js/mega-menu.js",    {"defer": "defer"}),
]
# The landing page is a full-width, hand-designed layout: no primary sidebar
# (the secondary one is removed via file metadata in index.rst).
html_sidebars         = {"index": []}
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

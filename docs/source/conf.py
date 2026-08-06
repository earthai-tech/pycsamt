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

import json
import os
import re
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
# Local Sphinx extensions (gallery_hub, ...).
sys.path.insert(0, os.path.abspath("_ext"))

# -- Project information -------------------------------------------------------
project = "pyCSAMT"
author = "Kouadio K. Laurent"
copyright = "2022-2026, earthai-tech"

try:
    from pycsamt import __version__ as release

    version = ".".join(release.split(".")[:2])  # e.g. "2.0"
except Exception:
    version = release = "2.0"

switcher_version_match = f"v{version.split('.')[0]}"

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
    "myst_parser",  # parse .md files (e.g. CHANGELOG, README)
    "sphinx_copybutton",  # copy-button on code blocks
    "sphinx_design",  # grid / tab / card directives
    "sphinx_toggleprompt",  # show/hide >>> prompts on doctest blocks
    "sphinx_gallery.gen_gallery",  # execute example scripts, embed output
    # Local (docs/source/_ext)
    "gallery_hub",  # compact animated-card landing page for examples/
    "public_api_catalog",  # grouped public modules/classes/functions tables
    "code_dropdown",  # accessible collapsible long-source listings
]

# -- sphinx-copybutton -----------------------------------------------------
# Strip interpreter prompts when copying code (>>> , ... , $ ).
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

# -- sphinx-toggleprompt ---------------------------------------------------
# Show Python console prompts (>>> / ...) by default so reviewer-facing
# doctest examples read like interactive sessions. Readers can still hide
# prompts with the toggle button, and sphinx-copybutton still strips prompts
# from copied code.
toggleprompt_default_hidden = "false"
# Shift the toggle left so it does not overlap the copy button.
toggleprompt_offset_right = 40

# -- sphinx-gallery ---------------------------------------------------------
# One site-wide gallery, scikit-learn style: ``docs/examples`` is the source
# tree, each subdirectory (emtools/, site/, ...) is a section with its own
# README.txt header, and the generated gallery lands at ``examples/`` so the
# "Examples" navbar entry *is* the card-grid landing page. To add a new
# section, create ``docs/examples/<name>/README.txt`` plus ``plot_*.py``
# scripts and list the section in ``subsection_order`` below.
from sphinx_gallery.sorting import (
    ExplicitOrder,
    FileNameSortKey,
)

sphinx_gallery_conf = {
    "examples_dirs": ["../examples"],
    "gallery_dirs": ["examples"],
    # Section order on the Examples landing page; "*" catches any section
    # added later that is not listed explicitly yet.
    "subsection_order": ExplicitOrder(
        [
            "../examples/survey_diagnostics",
            "../examples/emtools",
            "../examples/forward_modelling",
            "../examples/ai_inversion",
            "../examples/map_3d",
            "../examples/interp",
            "../examples/site",
            "../examples/topo",
            "../examples/corrections",
            "../examples/inversion",
            "../examples/pipeline",
            "../examples/agents",
            "../examples/iot",
            "*",
        ]
    ),
    # Card order *within* each section: by file name, so numeric-prefixed
    # scripts (forward_modelling/plot_1_*, plot_2_* …) read in workflow order.
    "within_subsection_order": FileNameSortKey,
    # Interactive Plotly maps (map_3d section) are embedded via capture_repr:
    # each cell ends with a bare ``fig`` and sphinx-gallery captures its
    # ``_repr_html_`` as an inline interactive figure. Card thumbnails for
    # those pages come from static PNGs via ``# sphinx_gallery_thumbnail_path``.
    "capture_repr": ("_repr_html_", "__repr__"),
    # Match "plot_*" regardless of path separator — the sphinx-gallery
    # default (r"/plot_") only matches POSIX paths and silently executes
    # zero files on Windows.
    "filename_pattern": r"[\\/]plot_",
    "download_all_examples": False,
    "plot_gallery": True,
    "remove_config_comments": True,
    "backreferences_dir": None,
    # Skip any underscore-prefixed helper module (_datasets.py, _synthetic.py,
    # __init__.py, ...) — they are imported by plot_*.py scripts, not
    # standalone examples. NOTE: sections that ship a same-named helper
    # (e.g. two _datasets.py) share it via sys.modules across the whole
    # build, so such helpers must be kept byte-identical supersets.
    "ignore_pattern": r"(^|[\\/])_",
    # NOTE: do NOT enable "parallel" here. The inversion section is an
    # ordered chain (plot_4 writes a validation report that plot_5 reads,
    # plot_5 writes the run-folder audit that plot_6 reads); parallel
    # workers execute out of order and break it on fresh CI machines.
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "_templates", "Thumbs.db", ".DS_Store"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"

# -- viewcode --------------------------------------------------------------------
# Only highlight pycsamt's own modules. The default (True) follows every
# imported object and ends up syntax-highlighting third-party internals
# (pandas.core.frame, ...) — minutes of build time for pages nobody links to.
viewcode_follow_imported_members = False

# -- autodoc / autosummary -----------------------------------------------------
autosummary_generate = True  # stubs not yet written; generate on build
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"

# -- Napoleon (NumPy-style docstrings) -----------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True

# -- numpydoc ------------------------------------------------------------------
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = False
numpydoc_class_members_toctree = False
numpydoc_validation_checks = set()

# -- intersphinx ---------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "sklearn": ("https://scikit-learn.org/stable", None),
}

# -- MathJax -------------------------------------------------------------------
mathjax3_config = {
    "tex": {"tags": "ams"},
}

# -- todo ----------------------------------------------------------------------
todo_include_todos = True  # flip to False before public release

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
        "image_dark": "_static/logo/pycsamt-v2-symbol.svg",
        "text": "pyCSAMT",
        "alt_text": "pyCSAMT",
    },
    # Navbar: version switcher sits next to the brand (numpy/pandas style);
    # search + icon links + theme toggle stay inline at the far right so
    # nothing wraps into a second header row.
    "navbar_start": ["navbar-logo", "version-switcher"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links", "theme-switcher"],
    "navbar_persistent": ["search-button"],
    "header_links_before_dropdown": 6,
    # Icon links (fill in real URLs when repo is public)
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/earthai-tech/pycsamt",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "Issue tracker",
            "url": "https://github.com/earthai-tech/pycsamt/issues",
            "icon": "fa-solid fa-bug",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/pycsamt/",
            "icon": "fa-brands fa-python",
        },
        {
            "name": "Stack Overflow",
            "url": "https://stackoverflow.com/questions/tagged/pycsamt",
            "icon": "fa-brands fa-stack-overflow",
        },
    ],
    # Sidebar
    "show_toc_level": 2,
    # Per-page override lives in _configure_secondary_sidebar() below; this is
    # the fallback for any page that handler does not touch.
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    # NOTE: do NOT add "sidebar-collapse" here. That component only exists
    # in unreleased pydata-sphinx-theme builds (>0.16.1) that also inject
    # their own collapse toggle at the top of the primary sidebar; on the
    # latest PyPI release (0.16.1, see the "docs" extra pin) the component
    # template does not exist at all, and referencing it here throws
    # "'sidebar-collapse.html' not found" and hard-fails the whole build.
    "primary_sidebar_end": [],
    # Code highlighting.
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
    # Footer
    "footer_start": ["copyright"],
    # "netlify-badge" satisfies the Netlify Open Source Plan's requirement
    # of a link back to Netlify on every page (see _templates/netlify-badge.html).
    "footer_end": ["sphinx-version", "theme-version", "netlify-badge"],
    # Search
    "search_bar_text": "Search pyCSAMT docs...",
    "switcher": {
        # Absolute on purpose: the theme resolves a relative json_url
        # against the address-bar URL, so directory-style URLs like
        # /examples/ 404 and the pill sticks at "Choose version".
        "json_url": "https://pycsamt.org/_static/version_switcher.json",
        "version_match": switcher_version_match,
    },
}

# Release-note badges (sphinx-design bdg roles), usable in any page as
# |Feature| |Fix| ... — same convention as the geoprior-v3 docs.
rst_epilog = """
.. |Feature| replace:: :bdg-success:`Feature`
.. |New| replace:: :bdg-success:`New`
.. |Fix| replace:: :bdg-info:`Fix`
.. |Enhancement| replace:: :bdg-info:`Enhancement`
.. |Perf| replace:: :bdg-info:`Perf`
.. |Breaking| replace:: :bdg-danger:`Breaking`
.. |Security| replace:: :bdg-danger:`Security`
.. |API Change| replace:: :bdg-warning:`API Change`
.. |Deprecated| replace:: :bdg-warning:`Deprecated`
.. |Docs| replace:: :bdg-secondary:`Docs`
.. |Build| replace:: :bdg-primary:`Build`
.. |Tests| replace:: :bdg-primary:`Tests`

.. |Compatibility| raw:: html

   <span class="sd-badge pyc-bg-compat">Compatibility</span>
"""

# Canonical URL of the published docs (Netlify + custom domain); the theme
# emits <link rel="canonical"> tags from this for search engines.
html_baseurl = "https://pycsamt.org/"

html_title = f"pyCSAMT {version}"
html_short_title = "pyCSAMT"
html_favicon = "_static/logo/pycsamt-v2-symbol.ico"
html_static_path = ["_static"]
# Static assets are organised by kind: _static/css/, _static/js/, _static/logo/
# (paths below are relative to the _static output root).
html_css_files = [
    "css/custom.css",
    "css/pycsamt-home.css",
    "css/gallery.css",
    "css/gallery-hub.css",
    "css/code-action.css",
    "css/code-dropdown.css",
    "css/hosted-preview.css",
    "css/whats-new.css",
    "css/faq.css",
]
html_js_files = [
    ("js/pycsamt-home.js", {"defer": "defer"}),
    ("js/mega-menu.js", {"defer": "defer"}),
    ("js/gallery-hub.js", {"defer": "defer"}),
    ("js/code-action.js", {"defer": "defer"}),
    ("js/code-dropdown.js", {"defer": "defer"}),
    ("js/api-search.js", {"defer": "defer"}),
    ("js/hosted-preview.js", {"defer": "defer"}),
    ("js/whats-new.js", {"defer": "defer"}),
    ("js/faq.js", {"defer": "defer"}),
]
# The landing page is a full-width, hand-designed layout: no primary sidebar
# (the secondary one is removed via file metadata in index.rst).
html_sidebars = {"index": []}
html_show_sourcelink = True
html_show_sphinx = False
html_copy_source = False


# Pages that render full-width, with no secondary ("On this page") sidebar at
# all: the hand-designed home page, the API reference whose object table needs
# every pixel, and the generated examples-gallery grids (matched by pattern).
_NO_SECONDARY_SIDEBAR = {"index", "api/index", "faq"}

# The support card is a site-wide appeal, not page content: it belongs on the
# section landings people arrive at, not on all ~300 leaf pages.
_SIDEBAR_WITH_SUPPORT = ["page-toc", "support-box", "edit-this-page"]
_SIDEBAR_PLAIN = ["page-toc", "edit-this-page"]

# Reference pages built from a directive rather than from sections: they have
# no headings, so their page-toc renders empty and -- with edit-this-page also
# empty -- the theme would drop the sidebar entirely. Carrying the support card
# keeps the column, which is why they are listed here rather than left to the
# structural test below.
_SUPPORT_ON_HEADINGLESS = {"glossary"}


def _has_guide_outline(doctree):
    """True when the page body is one of the numbered guide outlines.

    Structural rather than name-based on purpose: ``resources`` is every bit a
    section landing but is not called ``index``, and a name test silently left
    it with an empty page-toc, no support card, and therefore -- since the
    theme drops a sidebar with nothing in it -- no right sidebar at all.
    """
    if doctree is None:
        return False
    from docutils import nodes

    return any(
        "pycsamt-guide-toc" in node.get("classes", [])
        for node in doctree.findall(nodes.compound)
    )


def _is_landing_page(pagename, doctree):
    """True for pages that should carry the support card: any page whose body
    is a guide outline, the API landing (a grid of cards instead), and the
    heading-less reference pages that would otherwise lose their sidebar."""
    return (
        pagename == "api_landing"
        or pagename in _SUPPORT_ON_HEADINGLESS
        or _has_guide_outline(doctree)
    )


def _configure_secondary_sidebar(
    app, pagename, templatename, context, doctree
):
    """Decide the right-hand sidebar per page: drop it, or pick its contents.

    Runs before the theme's own ``set_secondary_sidebar_items`` (priority 400
    vs the default 500) so the per-page item list below is the one the theme
    reads.  Doing this here rather than via a ``secondary_sidebar_items`` dict
    is deliberate: the theme warns on every page matching two wildcard
    patterns, and "all pages" plus "index pages" necessarily overlap.
    """
    if pagename in _NO_SECONDARY_SIDEBAR or (
        pagename.startswith("examples/") and pagename.endswith("index")
    ):
        # context["meta"] can be present but None (e.g. templated/special
        # pages), so a plain .get("meta", {}) is not enough — normalise it.
        meta = context.get("meta") or {}
        context["meta"] = meta
        meta["html_theme.sidebar_secondary.remove"] = ""
        return

    context["theme_secondary_sidebar_items"] = (
        _SIDEBAR_WITH_SUPPORT
        if _is_landing_page(pagename, doctree)
        else _SIDEBAR_PLAIN
    )


# "What's new" badge (see _static/js/whats-new.js): docs/source/changelog/
# (one file per major-version series, e.g. v2.rst) is the single source of
# truth for the latest released version and its date, so the badge data is
# derived from it at build time rather than hand-kept in a separate JSON
# file that could go stale. A version section only counts once it carries a
# "*Released YYYY-MM-DD.*" line -- see the Pre-release checklist in
# development/documentation_build.rst. Every series file is scanned (not
# just the current one) and the highest released version wins, so this
# logic needs no manual update when a new series file appears.
_CHANGELOG_DIR = os.path.join(os.path.dirname(__file__), "changelog")
_CHANGELOG_VERSION_RE = re.compile(
    r"^(?P<version>\d+\.\d+\.\d+)\s+(?:\|[^|]+\|\s*)+$", re.MULTILINE
)
_CHANGELOG_RELEASED_RE = re.compile(
    r"^\*Released (?P<date>\d{4}-\d{2}-\d{2})\.\*\s*$", re.MULTILINE
)


def _version_key(version):
    # _CHANGELOG_VERSION_RE only ever captures a plain "\d+.\d+.\d+" (an
    # rc-suffixed heading like "2.0.0rc1 |Feature|" never matches -- there's
    # no whitespace between the digits and "rc1"), so a plain numeric tuple
    # is all that's needed here.
    return tuple(int(part) for part in version.split("."))


def _latest_release_info():
    try:
        paths = sorted(os.listdir(_CHANGELOG_DIR))
    except OSError:
        return None

    best = None
    for name in paths:
        if not name.endswith(".rst"):
            continue
        try:
            with open(os.path.join(_CHANGELOG_DIR, name), encoding="utf-8") as fh:
                text = fh.read()
        except OSError:
            continue

        for version_match in _CHANGELOG_VERSION_RE.finditer(text):
            # The release date, if present, immediately follows the version
            # heading -- search only that window so an older section's date
            # can never be mistaken for this one's.
            window = text[version_match.end() : version_match.end() + 500]
            released_match = _CHANGELOG_RELEASED_RE.search(window)
            if not released_match:
                continue

            version = version_match.group("version")
            if best is None or _version_key(version) > _version_key(best["version"]):
                best = {
                    "version": version,
                    "date": released_match.group("date"),
                    "url": f"release_notes/v{version}.html",
                }

    return best


def _write_whats_new_json(app, exception):
    if exception is not None:
        return
    info = _latest_release_info()
    if info is None:
        return
    out_dir = os.path.join(app.outdir, "_static", "data")
    os.makedirs(out_dir, exist_ok=True)
    with open(
        os.path.join(out_dir, "whats_new.json"), "w", encoding="utf-8"
    ) as fh:
        json.dump(info, fh, indent=2)


def setup(app):
    app.connect(
        "html-page-context", _configure_secondary_sidebar, priority=400
    )
    app.connect("build-finished", _write_whats_new_json)


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

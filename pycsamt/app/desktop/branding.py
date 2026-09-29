# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
branding — single source of truth for pycsamt-desktop's identity.

App id/name/version, resource paths, and external links used to drift
independently across ``__main__.py``, ``about_dialog.py``, and (soon) the
splash screen and installer scripts — each hardcoding its own copy of the
version string or icon path. Import from here instead of re-deriving any of
this locally.

Deliberately import-light (only ``pathlib`` and an optional ``pycsamt``
import) so this module stays safe to import very early — before
``QApplication`` exists, from a frozen build's bootstrap code, or from a
packaging script that has no Qt installed at all.
"""

from __future__ import annotations

from pathlib import Path

APP_ID = "pycsamt"
APP_NAME = "pycsamt"
APP_DISPLAY_NAME = "pyCSAMT"
ORG_NAME = "earthai-tech"
ORG_DOMAIN = "pycsamt.org"
TAGLINE = "Python  MT · AMT · CSAMT · CSEM  —  Geophysical Processing Suite"

AUTHOR_NAME = "Laurent Kouadio"
AUTHOR_EMAIL = "etanoyau@gmail.com"
LICENSE_SPDX = "LGPL-3.0"

URL_DOCS = "https://pycsamt.org/"
URL_GITHUB = "https://github.com/earthai-tech/pycsamt"

_FALLBACK_VERSION = "0.0.0"


def get_version() -> str:
    """Return the installed ``pycsamt`` version, or a clearly-fake fallback.

    The fallback is deliberately not a plausible version number (unlike the
    historical ``"2.0.0"`` fallbacks this replaces) so a broken import shows
    up as obviously wrong rather than silently looking like an old release.
    """
    try:
        import pycsamt

        return getattr(pycsamt, "__version__", _FALLBACK_VERSION)
    except Exception:
        return _FALLBACK_VERSION


# ── Resource paths ──────────────────────────────────────────────────────────
RESOURCES_DIR = Path(__file__).parent / "resources"
ICONS_DIR = RESOURCES_DIR / "icons"
LOGO_ICO = ICONS_DIR / "pycsamt.logo.ico"
LOGO_SVG = ICONS_DIR / "pycsamt_logo.svg"

# Copied from docs/source/_static/applications/desktop/pycsamt-v2-splash.png
# -- kept in sync manually, not symlinked (the docs build and the packaged
# app ship from different trees). Phase 12 wires this into a QSplashScreen
# shown in __main__.py before the heavy-import block; nothing shows it yet.
SPLASH_IMAGE = RESOURCES_DIR / "pycsamt-v2-splash.png"

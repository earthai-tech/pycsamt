# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt v2 — Python toolbox for audio-magnetotellurics (AMT & CSAMT).

Documentation: https://pycsamt.readthedocs.io/en/latest/
Source:        https://github.com/WEgeophysics/pycsamt
"""

import importlib
import warnings

# ── Version ──────────────────────────────────────────────────────────────────
try:
    from importlib.metadata import version as _pkg_version
    __version__ = _pkg_version(__name__)
except Exception:
    __version__ = "2.0.0"

# ── Logging ───────────────────────────────────────────────────────────────────
from pycsamt.log._config import init_logging
init_logging()

# ── Hard-dependency check ────────────────────────────────────────────────────
_required = [
    ("numpy",       "numpy"),
    ("pandas",      "pandas"),
    ("scipy",       "scipy"),
    ("matplotlib",  "matplotlib"),
    ("tqdm",        "tqdm"),
    ("scikit-learn","sklearn"),
]
_missing = []
for _pkg_name, _import_name in _required:
    try:
        importlib.import_module(_import_name)
    except ImportError:
        _missing.append(_pkg_name)

if _missing:
    warnings.warn(
        "pycsamt: some dependencies are not installed; "
        "affected features will raise ImportError at use time:\n  "
        + ", ".join(_missing),
        ImportWarning,
        stacklevel=2,
    )
del _required, _pkg_name, _import_name, _missing

# ── Eagerly loaded subpackages ────────────────────────────────────────────────
# Only the two packages that are always needed at the module level are loaded
# here; everything else is imported on demand by the user.
from . import interp  # noqa: F401  geological interpretation / export
from . import inversion  # noqa: F401  physics-based EM inversion API
from . import tdem    # noqa: F401  time-domain EM → EDICollection

# ── Backend API ───────────────────────────────────────────────────────────────
from pycsamt.backends import (
    auto_detect,
    get_backend,
    get_backend_instance,
    list_backends,
    set_backend,
)

# ── Graceful errors for removed v1 names ─────────────────────────────────────
_REMOVED = {
    "geodrill":    "Use 'pycsamt.interp' instead.",
    "bases":       "Removed in v2 (was an internal v1 serialisation helper).",
    "_csamtpylog": "Removed in v2; use 'pycsamt.log.logger.get_logger'.",
    "_path":       "Removed in v2; use standard pathlib paths.",
}


def __getattr__(name: str):
    if name in _REMOVED:
        raise AttributeError(
            f"'pycsamt.{name}' was removed in pycsamt v2.  {_REMOVED[name]}"
        )
    raise AttributeError(f"module 'pycsamt' has no attribute {name!r}")


# ── Public API ────────────────────────────────────────────────────────────────
__all__ = [
    "__version__",
    # subpackages (importable as `import pycsamt.seg` etc.)
    "ai",
    "backends",
    "emtools",
    "forward",
    "gis",
    "interp",
    "inversion",
    "io",
    "jones",
    "models",
    "seg",
    "site",
    "tdem",
    "z",
    "zonge",
    # backend helpers
    "auto_detect",
    "get_backend",
    "get_backend_instance",
    "list_backends",
    "set_backend",
]

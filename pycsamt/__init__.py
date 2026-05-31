# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: GPL-3.0
"""
pycsamt: A Python Toolbox for audio-magnetotellurics (AMT & CSAMT)

Provides fast & robust data processing, signal recovery, and
**drilling-location optimization** for AMT and CSAMT surveys.

Documentation: https://pycsamt.readthedocs.io/en/latest/
Source:        https://github.com/WEgeophysics/pycsamt
"""

import sys
import warnings
import importlib
import types

# Version
try:
    from importlib.metadata import version
    __version__ = version(__name__)
except Exception:
    __version__ = "2.0.0"

# Configure logging early

from pycsamt.log._config import init_logging
init_logging()  # loads p.configlog.yml or falls back to basicConfig


# Dependency resolution
_required_dependencies = [
    ("numpy", None),
    ("pandas", None),
    ("scipy", None),
    ("matplotlib", None),
    ("tqdm", None),
    ("scikit-learn", "sklearn"),
    ("numba", None),
    ("mtpy", None),
    ("pyyaml", None),
]
_missing_deps = []

def _lazy_import(module, alias=None):
    """Lazily load a module to defer import time."""
    def _loader():
        return importlib.import_module(module)
    globals()[alias or module] = _loader

for pkg, alias in _required_dependencies:
    try:
        if alias:
            _lazy_import(alias, alias)
        else:
            _lazy_import(pkg)
    except ImportError as e:
        _missing_deps.append(f"{pkg}: {e}")

if _missing_deps:
    warnings.warn(
        "Missing pycsamt dependencies; some features may be unavailable:\n"
        + "\n".join(_missing_deps),
        ImportWarning
    )

# Expose installer utility for on-the-fly fixes
from .utils.generic import ensure_package # Noqa 


# interp: geological interpretation and export (replaces geodrill from v1)
from . import interp  # noqa: F401

# tdem: time-domain EM → frequency-domain (EDICollection) conversion
from . import tdem  # noqa: F401

# __getattr__ for clean error on removed v1 modules
def __getattr__(name):
    if name == "geodrill":
        raise AttributeError(
            "'geodrill' was removed in pycsamt v2. "
            "Use 'pycsamt.interp' instead."
        )
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


# Backend API — always importable, no ML framework required
from pycsamt.backends import (
    get_backend, set_backend, auto_detect,
    list_backends, get_backend_instance,
)

# Public API
__all__ = [
    "__version__",
    "is_installing",
    # backend control
    "get_backend", "set_backend", "auto_detect",
    "list_backends", "get_backend_instance",
    "interp",
    "tdem",
]





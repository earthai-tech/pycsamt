# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt v2 — Python toolbox for audio-magnetotellurics (AMT & CSAMT).

Documentation: https://pycsamt.org/
Source:        https://github.com/earthai-tech/pycsamt
"""

from __future__ import annotations

import importlib

# ── Version ───────────────────────────────────────────────────────────────────
try:
    from importlib.metadata import version as _pkg_version

    __version__ = _pkg_version(__name__)
except Exception:
    __version__ = "2.3.1"

# ── Removed v1 names ──────────────────────────────────────────────────────────
# Defined BEFORE logging (and before __getattr__) so this dict is always
# present in the module dict, even if the logging import fails mid-reload.
_REMOVED = {
    "geodrill": "Use 'pycsamt.interp' instead.",
    "bases": "Removed in v2 (was an internal v1 serialisation helper).",
    "_csamtpylog": "Removed in v2; use 'pycsamt.log.logger.get_logger'.",
    "_path": "Removed in v2; use standard pathlib paths.",
}

# ── Lazy registry ─────────────────────────────────────────────────────────────
# All pure-Python constants — no external imports needed.
_SUBPACKAGES = frozenset(
    {
        "ai",
        "agents",
        "backends",
        "emtools",
        "forward",
        "geology",
        "gis",
        "interp",
        "inversion",
        "io",
        "jones",
        "log",
        "map",
        "models",
        "pipeline",
        "seg",
        "site",
        "tdem",
        "ts",
        "z",
        "zonge",
    }
)

_LAZY_SYMBOLS: dict[str, str] = {
    # pipeline
    "Pipeline": ".pipeline",
    "Step": ".pipeline",
    "configure_pipe": ".pipeline",
    "reset_pipe": ".pipeline",
    "PYCSAMT_PIPE": ".pipeline",
    # backends
    "auto_detect": ".backends",
    "get_backend": ".backends",
    "get_backend_instance": ".backends",
    "list_backends": ".backends",
    "set_backend": ".backends",
    # workflow sessions
    "Session": ".session",
    "work_session": ".session",
    "Normalize": ".session",
    "normalize_session": ".session",
}


def __getattr__(name: str):
    if name in _SUBPACKAGES:
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod

    if name in _LAZY_SYMBOLS:
        mod = importlib.import_module(_LAZY_SYMBOLS[name], __name__)
        sym = getattr(mod, name)
        globals()[name] = sym
        return sym

    if name in _REMOVED:
        raise AttributeError(
            f"'pycsamt.{name}' was removed in pycsamt v2.  {_REMOVED[name]}"
        )
    raise AttributeError(f"module 'pycsamt' has no attribute {name!r}")


# ── Logging ───────────────────────────────────────────────────────────────────
# Placed AFTER __getattr__ and _REMOVED so any import side-effect that
# triggers __getattr__ (e.g. via Werkzeug hot-reload) is always safe.
# Wrapped in try/except: yaml not yet installed won't break the import.
try:
    from .log._config import init_logging as _init_logging

    _init_logging()
    del _init_logging
except Exception:
    pass


# ── Public API ────────────────────────────────────────────────────────────────
__all__ = [
    "__version__",
    # subpackages
    "ai",
    "agents",
    "backends",
    "emtools",
    "forward",
    "geology",
    "gis",
    "interp",
    "inversion",
    "io",
    "jones",
    "map",
    "models",
    "pipeline",
    "seg",
    "site",
    "tdem",
    "ts",
    "z",
    "zonge",
    # pipeline shortcuts
    "Pipeline",
    "Step",
    "configure_pipe",
    "reset_pipe",
    "PYCSAMT_PIPE",
    # backend helpers
    "auto_detect",
    "get_backend",
    "get_backend_instance",
    "list_backends",
    "set_backend",
    # workflow sessions
    "Session",
    "work_session",
    "Normalize",
    "normalize_session",
]

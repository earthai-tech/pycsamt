# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Backend detection utilities — zero ML-framework imports.

All probing is done via :mod:`importlib.util` so neither PyTorch nor
TensorFlow is actually imported (and not required as a package-level
dependency) when this module is loaded.
"""

from __future__ import annotations

import importlib.util

__all__ = [
    "detect_available_backends",
    "probe_backend",
    "get_backend_versions",
]


def probe_backend(name: str) -> bool:
    """
    Return ``True`` if *name* is importable in the current environment.

    No import side-effects — uses :func:`importlib.util.find_spec` only.

    Parameters
    ----------
    name : str
        ``'torch'`` or ``'tensorflow'``.

    Returns
    -------
    bool
    """
    _aliases = {
        "torch": "torch",
        "pytorch": "torch",
        "tensorflow": "tensorflow",
        "tf": "tensorflow",
        "keras": "keras",
    }
    pkg = _aliases.get(name.lower(), name)
    return importlib.util.find_spec(pkg) is not None


def detect_available_backends() -> list[str]:
    """
    Return a list of available deep-learning backends.

    Returns
    -------
    backends : list of str
        Ordered list of backend names that are currently importable.
        Possible values: ``'torch'``, ``'tensorflow'``.
        Empty list if neither framework is installed.

    Examples
    --------
    >>> from pycsamt.backends._detect import detect_available_backends
    >>> available = detect_available_backends()
    >>> print(available)  # e.g. ['torch'] or ['torch', 'tensorflow']
    """
    found: list[str] = []
    if probe_backend("torch"):
        found.append("torch")
    if probe_backend("tensorflow"):
        found.append("tensorflow")
    return found


def get_backend_versions() -> dict[str, str | None]:
    """
    Return version strings for installed backends.

    Only actually imports the framework if it was already detected;
    avoids triggering slow framework initialisation.

    Returns
    -------
    versions : dict
        ``{'torch': '2.x.x', 'tensorflow': None}`` — ``None`` means
        not installed.
    """
    result: dict[str, str | None] = {
        "torch": None,
        "tensorflow": None,
    }
    if probe_backend("torch"):
        try:
            import torch

            result["torch"] = torch.__version__
        except Exception:
            pass
    if probe_backend("tensorflow"):
        try:
            import tensorflow as tf

            result["tensorflow"] = tf.__version__
        except Exception:
            pass
    return result

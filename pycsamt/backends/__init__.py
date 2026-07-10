# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.backends
================

Backend abstraction layer for deep-learning frameworks.

pycsamt's AI/ML modules support both PyTorch and TensorFlow.  Neither
framework is a required dependency of the package — they are loaded
lazily only when the AI module is actually used.

Configuration
-------------
The active backend is resolved using this priority chain:

1. Explicit call to :func:`set_backend` in the current Python session.
2. ``PYCSAMT_AI_BACKEND`` environment variable
   (``torch`` / ``tensorflow`` / ``auto``).
3. ``~/.pycsamt/config.json`` key ``"ai_backend"``.
4. Auto-detection: first available framework in ``['torch', 'tensorflow']``.

Quick start
-----------
>>> import pycsamt
>>> pycsamt.list_backends()
{'torch': True, 'tensorflow': False}
>>> pycsamt.set_backend('torch')      # or 'tensorflow' / 'auto'
>>> pycsamt.get_backend()
'torch'

Environment variable
--------------------
Set before importing pycsamt to override the default::

    PYCSAMT_AI_BACKEND=tensorflow python my_script.py
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from ._config import BackendConfig
from ._detect import (
    detect_available_backends,
    get_backend_versions,
    probe_backend,
)

__all__ = [
    "get_backend",
    "set_backend",
    "auto_detect",
    "list_backends",
    "get_backend_instance",
    "detect_available_backends",
    "probe_backend",
    "get_backend_versions",
]

_CFG = BackendConfig()


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def get_backend() -> str:
    """
    Return the name of the currently active backend.

    Returns
    -------
    name : str
        ``'torch'``, ``'tensorflow'``, or ``'none'`` when no DL framework
        is installed.

    Examples
    --------
    >>> from pycsamt.backends import get_backend
    >>> get_backend()
    'torch'
    """
    return _CFG.backend_name


def set_backend(name: str, *, persist: bool = False) -> None:
    """
    Set the active AI backend.

    Parameters
    ----------
    name : {'torch', 'tensorflow', 'auto'}
        ``'auto'`` triggers detection and selects the first available
        framework.
    persist : bool, default False
        If ``True``, save the choice to ``~/.pycsamt/config.json`` so
        it survives across Python sessions.

    Raises
    ------
    ValueError
        If *name* is not ``'torch'``, ``'tensorflow'``, or ``'auto'``.
    ImportError
        If the requested backend is not installed.

    Examples
    --------
    >>> from pycsamt.backends import set_backend
    >>> set_backend('torch')
    >>> set_backend('tensorflow', persist=True)   # saves to config file
    """
    _CFG.set(name)
    if persist:
        _CFG.write_config_file(_CFG.backend_name)


def auto_detect() -> str:
    """
    Detect the best available backend and activate it.

    Returns
    -------
    name : str
        The name of the detected and activated backend.

    Raises
    ------
    RuntimeError
        If no compatible framework is installed.

    Examples
    --------
    >>> from pycsamt.backends import auto_detect
    >>> backend = auto_detect()
    >>> print(backend)
    torch
    """
    _CFG.reset()
    return _CFG.backend_name


def list_backends() -> dict[str, bool]:
    """
    Return availability status for all known backends.

    Returns
    -------
    availability : dict
        ``{'torch': bool, 'tensorflow': bool}``

    Examples
    --------
    >>> from pycsamt.backends import list_backends
    >>> list_backends()
    {'torch': True, 'tensorflow': False}
    """
    return {
        "torch": probe_backend("torch"),
        "tensorflow": probe_backend("tensorflow"),
    }


def get_backend_instance() -> Any:
    """
    Return the concrete :class:`~pycsamt.backends._base.NeuralBackend`
    instance for the active backend, or ``None`` when no DL framework is
    installed.

    Returns
    -------
    backend : NeuralBackend or None
        :class:`~pycsamt.backends._torch.TorchBackend`,
        :class:`~pycsamt.backends._tensorflow.TensorFlowBackend`, or
        ``None`` when no framework is available.
    """
    name = get_backend()
    if name == "torch":
        from ._torch import TorchBackend
        return TorchBackend()
    if name == "tensorflow":
        from ._tensorflow import TensorFlowBackend
        return TensorFlowBackend()
    if name == "none":
        return None
    raise RuntimeError(f"Unknown backend {name!r}")

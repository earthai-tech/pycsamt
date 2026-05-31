# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Backend-aware utility helpers shared by AI processing and inversion modules.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

__all__ = [
    "resolve_device",
    "get_weights",
    "set_weights",
    "active_backend",
]


def resolve_device(device: Optional[str] = None) -> str:
    """
    Return a device string appropriate for the active backend.

    Falls back to ``'cpu'`` / ``'/CPU:0'`` when no backend is available.
    """
    try:
        from pycsamt.backends import get_backend_instance
        return get_backend_instance().resolve_device(device)
    except RuntimeError:
        return device or "cpu"


def get_weights(model: Any) -> Dict[str, np.ndarray]:
    """Return model weights via the active backend."""
    from pycsamt.backends import get_backend_instance
    return get_backend_instance().get_weights(model)


def set_weights(model: Any, weights: Dict[str, np.ndarray]) -> None:
    """Restore model weights via the active backend."""
    from pycsamt.backends import get_backend_instance
    get_backend_instance().set_weights(model, weights)


def active_backend() -> str:
    """
    Return the name of the active backend (``'torch'`` or ``'tensorflow'``).

    Raises ``RuntimeError`` if no backend is installed.
    """
    from pycsamt.backends import get_backend
    return get_backend()

# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Backend registry for :mod:`pycsamt.inversion`."""

from __future__ import annotations

from importlib import import_module

__all__ = ["available_backends", "get_backend"]


_BACKENDS = {
    "builtin": ("pycsamt.inversion.backends.builtin", "Builtin1DBackend"),
    "occam2d": ("pycsamt.inversion.backends.occam2d", "Occam2DBackend"),
    "modem": ("pycsamt.inversion.backends.modem", "ModEMBackend"),
    "simpeg": ("pycsamt.inversion.backends.simpeg", "SimPEGBackend"),
    "pygimli": ("pycsamt.inversion.backends.pygimli", "PyGIMLiBackend"),
}


def available_backends() -> list[str]:
    """Return registered backend names."""
    return sorted(_BACKENDS)


def get_backend(name: str):
    """Return the backend class for *name* using lazy imports."""
    key = str(name).lower()
    if key not in _BACKENDS:
        raise KeyError(f"Unknown inversion backend {name!r}.")
    module_name, class_name = _BACKENDS[key]
    module = import_module(module_name)
    return getattr(module, class_name)

# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.models
==============

Inversion model subpackages.  Each subpackage wraps an external inversion
code and provides a Python interface that consumes EDI files (single source
of truth) and produces ready-to-run inversion inputs, a runner that calls
the compiled binary, and a post-processing layer for results.

Subpackages
-----------
occam2d
    2-D smooth-model MT inversion using the Occam2DMT v3.0 Fortran code
    (Constable et al. 1987 / deGroot-Hedlin & Constable 1990).
mare2dem
    2.5-D FEM MT and CSEM inversion using MARE2DEM (Key 2016). Source
    is not bundled; use ``SourceManager`` to download and build. Includes
    download/build helpers, input builders, subprocess runner, result
    loader, and plot utilities.
"""

from importlib import import_module as _imp
import types as _types
import sys as _sys

_subpackages = ["occam2d", "mare2dem"]

for _pkg in _subpackages:
    try:
        globals()[_pkg] = _imp(f"pycsamt.models.{_pkg}")
    except ImportError:
        _dummy = _types.ModuleType(f"pycsamt.models.{_pkg}")
        _sys.modules[f"pycsamt.models.{_pkg}"] = _dummy
        globals()[_pkg] = _dummy

__all__ = _subpackages

# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""2-D TDEM inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["TDEM2DInversion"]


class TDEM2DInversion(ConfiguredInversionWorkflow):
    """Convenience workflow for TDEM 2-D inversion backends."""

    method = "tdem"
    dimension = "2d"
    backend = "builtin"

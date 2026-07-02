# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""2-D MT inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["MT2DInversion"]


class MT2DInversion(ConfiguredInversionWorkflow):
    """Built-in stitched 2-D MT inversion workflow."""

    method = "mt"
    dimension = "2d"
    backend = "builtin"

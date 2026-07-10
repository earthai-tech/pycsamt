# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""2-D AMT inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["AMT2DInversion"]


class AMT2DInversion(ConfiguredInversionWorkflow):
    """Convenience workflow preconfigured for built-in stitched 2-D AMT."""

    method = "amt"
    dimension = "2d"
    backend = "builtin"

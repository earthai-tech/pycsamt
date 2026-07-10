# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""2-D EMAP inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["EMAP2DInversion"]


class EMAP2DInversion(ConfiguredInversionWorkflow):
    """Convenience workflow preconfigured for built-in stitched 2-D EMAP."""

    method = "emap"
    dimension = "2d"
    backend = "builtin"

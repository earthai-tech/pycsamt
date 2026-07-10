# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""1-D AMT inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["AMT1DInversion"]


class AMT1DInversion(ConfiguredInversionWorkflow):
    """Convenience workflow preconfigured for built-in 1-D AMT inversion."""

    method = "amt"
    dimension = "1d"
    backend = "builtin"

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""1-D MT inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["MT1DInversion"]


class MT1DInversion(ConfiguredInversionWorkflow):
    """Built-in 1-D MT inversion workflow."""

    method = "mt"
    dimension = "1d"
    backend = "builtin"

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""3-D MT inversion convenience wrappers."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["MT3DInversion", "SimPEGMT3DInversion"]


class MT3DInversion(ConfiguredInversionWorkflow):
    """ModEM 3-D MT inversion workflow."""

    method = "mt"
    dimension = "3d"
    backend = "modem"


class SimPEGMT3DInversion(ConfiguredInversionWorkflow):
    """SimPEG 3-D MT inversion workflow."""

    method = "mt"
    dimension = "3d"
    backend = "simpeg"

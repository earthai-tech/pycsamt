# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""1-D TDEM inversion convenience wrapper."""

from __future__ import annotations

from .._method import ConfiguredInversionWorkflow

__all__ = ["TDEM1DInversion"]


class TDEM1DInversion(ConfiguredInversionWorkflow):
    """Convenience workflow for TDEM 1-D inversion backends.

    The default backend is ``"builtin"`` for a lightweight layered-earth
    least-squares inversion. Use ``backend="simpeg"`` when a SimPEG-based
    implementation is available in the environment.
    """

    method = "tdem"
    dimension = "1d"
    backend = "builtin"

# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Physics-based EM inversion API for pyCSAMT.

The :mod:`pycsamt.inversion` package provides a backend-neutral interface
for EM inversion. It includes runnable built-in layered-earth inversions
for MT/AMT/CSAMT and TDEM soundings, stitched station-by-station 2-D
sections, and adapter backends for existing Occam2D and ModEM workflows.
SimPEG and pyGIMLi are registered lazily for optional numerical backends.

Examples
--------
>>> from pycsamt.inversion import InversionConfig, InversionWorkflow
>>> cfg = InversionConfig(
...     method="mt",
...     dimension="1d",
...     backend="builtin",
...     data={"freqs": [1, 10], "rho_a": [100, 120], "phase": [45, 47]},
... )
>>> result = InversionWorkflow(cfg).run()  # doctest: +SKIP
"""

from .config import InversionConfig
from .data import EMData
from .mesh import InversionMesh
from .model import ReferenceModel, StartingModel
from .regularization import Regularization
from .results import InversionResult
from .workflow import InversionWorkflow, run_inversion
from .backends import available_backends, get_backend
from . import export, plot

__all__ = [
    "EMData",
    "InversionConfig",
    "InversionMesh",
    "InversionResult",
    "InversionWorkflow",
    "ReferenceModel",
    "Regularization",
    "StartingModel",
    "available_backends",
    "export",
    "get_backend",
    "plot",
    "run_inversion",
]

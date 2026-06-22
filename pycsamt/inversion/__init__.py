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
from .mesh import (
    InversionMesh,
    build_1d_tensor_mesh,
    build_3d_tensor_mesh,
    build_fd2d_grid,
    core_rho_from_start,
    depth_widths,
)
from .model import ReferenceModel, StartingModel
from .objective import (
    ErrorModel,
    component_errors,
    component_mask,
    error_model_from_config,
)
from .regularization import (
    Regularization,
    pygimli_lambda,
    regularization_from_config,
    regularization_residual,
    regularization_weight,
)
from .results import InversionHistory, InversionResult, InversionUncertainty
from .workflow import InversionWorkflow, run_inversion
from .backends import available_backends, get_backend
from .mt import MT1DInversion, MT2DInversion, MT3DInversion, SimPEGMT3DInversion
from . import export, plot

__all__ = [
    "EMData",
    "ErrorModel",
    "InversionHistory",
    "InversionConfig",
    "InversionMesh",
    "InversionResult",
    "InversionUncertainty",
    "InversionWorkflow",
    "ReferenceModel",
    "Regularization",
    "StartingModel",
    "available_backends",
    "MT1DInversion",
    "MT2DInversion",
    "MT3DInversion",
    "SimPEGMT3DInversion",
    "build_1d_tensor_mesh",
    "build_3d_tensor_mesh",
    "build_fd2d_grid",
    "component_errors",
    "component_mask",
    "core_rho_from_start",
    "depth_widths",
    "error_model_from_config",
    "export",
    "get_backend",
    "plot",
    "pygimli_lambda",
    "regularization_from_config",
    "regularization_residual",
    "regularization_weight",
    "run_inversion",
]

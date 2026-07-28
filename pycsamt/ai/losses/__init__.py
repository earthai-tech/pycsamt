# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Loss terms for model, spatial, response, boundary, and
uncertainty constraints.

Loss implementations belong here when they express inversion semantics
rather than a specific training loop.  Backend-specific adapters may
use NumPy, PyTorch, or TensorFlow internally, but importing this
package must remain safe when optional machine-learning dependencies
are unavailable.

Status
------
All five submodules are implemented: masked, weighted data-fit
losses on canonical resistivity grids (``L_model``), gradient/
total-variation regularizers (``L_grad_x``, ``L_grad_z``, ``L_TV``),
uncertainty-normalized complex-impedance residuals (``L_response``),
masked-target boundary constraints, and heteroscedastic Gaussian
NLL / coverage-calibration diagnostics for the M9 uncertainty
milestone in the AI-inversion plan.
"""

from .boundary import BoundaryLoss, boundary_condition_loss
from .model import (
    ModelLoss,
    ModelLossResult,
    depth_weights,
    model_huber_loss,
    model_l1_loss,
    model_l2_loss,
)
from .response import (
    ResponseLoss,
    ResponseLossResult,
    response_loss_from_contracts,
    response_residual_loss,
)
from .spatial import (
    SpatialLoss,
    SpatialLossResult,
    gradient_smoothness_loss,
    total_variation_loss,
)
from .uncertainty import (
    UncertaintyLoss,
    UncertaintyLossResult,
    calibration_loss,
    gaussian_nll_loss,
)

__all__ = [
    "ModelLossResult",
    "ModelLoss",
    "model_l1_loss",
    "model_l2_loss",
    "model_huber_loss",
    "depth_weights",
    "SpatialLossResult",
    "SpatialLoss",
    "gradient_smoothness_loss",
    "total_variation_loss",
    "ResponseLossResult",
    "ResponseLoss",
    "response_residual_loss",
    "response_loss_from_contracts",
    "BoundaryLoss",
    "boundary_condition_loss",
    "UncertaintyLossResult",
    "UncertaintyLoss",
    "gaussian_nll_loss",
    "calibration_loss",
]

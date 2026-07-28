# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Scientific validation gates for AI-assisted EM inversion.

This package contains synthetic recovery metrics, complex-response
residual diagnostics, uncertainty calibration, and out-of-distribution
checks.  Plotting remains in :mod:`pycsamt.ai.plot`; this package
owns the numerical evidence behind plots.

Status
------
All four submodules are implemented: masked RMSE/MAE/R2, structural
similarity (SSIM), and depth-profile recovery diagnostics for
known-truth synthetic grids (M0 baseline metrics); complex-impedance
residual reports broken down by station/frequency/component
(Inversion/Field rows of the validation matrix); Gaussian
reliability curves with coverage, calibration penalty, and sharpness
(Uncertainty row); and Mahalanobis/k-NN out-of-distribution scoring
against a reference realization set (OOD sensitivity).
"""

from .calibration import (
    ReliabilityCurve,
    empirical_coverage,
    predictive_sharpness,
    reliability_curve,
)
from .ood import OODReport, flag_out_of_distribution, ood_score
from .recovery import (
    RecoveryReport,
    depth_profile_mae,
    depth_profile_rmse,
    recovery_report,
    structural_similarity,
)
from .residuals import (
    ResponseResidualReport,
    response_residual_report,
    response_residual_report_from_contracts,
)

__all__ = [
    "RecoveryReport",
    "recovery_report",
    "structural_similarity",
    "depth_profile_rmse",
    "depth_profile_mae",
    "ResponseResidualReport",
    "response_residual_report",
    "response_residual_report_from_contracts",
    "ReliabilityCurve",
    "reliability_curve",
    "empirical_coverage",
    "predictive_sharpness",
    "OODReport",
    "ood_score",
    "flag_out_of_distribution",
]

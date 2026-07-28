# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Domain-gap and noise simulation for AI-assisted EM inversion (M3).

This package turns a clean :class:`~pycsamt.ai.data.contracts.SurveyData`
into realistic training data by injecting heteroscedastic noise, dropout,
static shift, galvanic distortion, coordinate perturbation, and outliers
(:mod:`~pycsamt.ai.domain_gap.simulator`); by fitting plausible parameter
ranges from a real survey's own QC diagnostics, AMT, CSAMT, MT, or
otherwise (:mod:`~pycsamt.ai.domain_gap.survey_fit`); and by comparing
simulated and field feature distributions quantitatively
(:mod:`~pycsamt.ai.domain_gap.report`).
"""

from __future__ import annotations

from .audit import (
    DimensionalitySummary,
    FrequencyGridReport,
    StationExclusion,
    SurveyAuditReport,
    audit_survey,
)
from .report import (
    DistributionComparisonReport,
    FeatureComparison,
    compare_feature_distributions,
    compare_survey_distributions,
)
from .simulator import (
    SEVERITY_PRESETS,
    CorruptionConfig,
    CorruptionRecord,
    add_heteroscedastic_noise,
    apply_corruption_suite,
    apply_dropout,
    apply_error_floor,
    apply_galvanic_distortion,
    apply_static_shift,
    inject_outliers,
    perturb_coordinates,
)
from .survey_fit import (
    fit_corruption_config,
    fit_distortion_priors_from_sites,
    survey_data_from_sites,
)

__all__ = [
    "StationExclusion",
    "FrequencyGridReport",
    "DimensionalitySummary",
    "SurveyAuditReport",
    "audit_survey",
    "SEVERITY_PRESETS",
    "CorruptionConfig",
    "CorruptionRecord",
    "add_heteroscedastic_noise",
    "apply_corruption_suite",
    "apply_dropout",
    "apply_error_floor",
    "apply_galvanic_distortion",
    "apply_static_shift",
    "inject_outliers",
    "perturb_coordinates",
    "fit_corruption_config",
    "fit_distortion_priors_from_sites",
    "survey_data_from_sites",
    "DistributionComparisonReport",
    "FeatureComparison",
    "compare_feature_distributions",
    "compare_survey_distributions",
]

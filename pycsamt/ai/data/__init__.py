# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Canonical data contracts for reproducible AI-assisted EM inversion.

This package owns survey tensors, validity masks, coordinate metadata,
normalisation state, realization-level dataset splits, and dataset manifests.
It complements :mod:`pycsamt.ai.training.dataset`, which owns framework-facing
training dataset wrappers.

The first implementation milestone will add dependency-light canonical data
objects here.  Importing this package intentionally has no optional machine-
learning or forward-solver side effects.
"""

from .contracts import (
    ImpedanceConvention,
    SurveyCoverage,
    SurveyData,
    merge_surveys,
)
from .manifest import (
    ArtifactRecord,
    DatasetManifest,
    canonical_hash,
    sha256_file,
)
from .normalization import ComplexZScore, NormalizedSurvey
from .splits import RealizationSplit, realization_folds, split_realizations

__all__ = [
    "SurveyData",
    "SurveyCoverage",
    "ImpedanceConvention",
    "merge_surveys",
    "ComplexZScore",
    "NormalizedSurvey",
    "RealizationSplit",
    "split_realizations",
    "realization_folds",
    "DatasetManifest",
    "ArtifactRecord",
    "canonical_hash",
    "sha256_file",
]

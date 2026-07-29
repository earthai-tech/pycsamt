# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reproducible experiment configuration and artifact provenance.

Experiment records will capture configuration hashes, seeds, realization
splits, frequency and depth grids, normalisation state, software versions,
metrics, and artifact checksums.  Runtime outputs and large artifacts are not
package resources and must remain outside source control.
"""

from .config import (
    AcceptanceCriterion,
    DatasetReference,
    ExperimentConfig,
    GateEvaluation,
    SeedPlan,
)

__all__ = [
    "SeedPlan",
    "DatasetReference",
    "AcceptanceCriterion",
    "GateEvaluation",
    "ExperimentConfig",
]

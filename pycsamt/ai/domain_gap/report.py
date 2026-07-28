# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Quantitative comparison between simulated and field feature distributions.

The M3 acceptance gate requires that "simulated and field feature
distributions are compared quantitatively" rather than judged by eye. This
module computes per-feature summary statistics and a two-sample
Kolmogorov-Smirnov test between a simulated (corrupted synthetic) survey
and a real field survey, on features derived from the shared
:class:`~pycsamt.ai.data.contracts.SurveyData` contract so the comparison
never depends on how either survey was produced.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..data.contracts import SurveyData

__all__ = [
    "FeatureComparison",
    "DistributionComparisonReport",
    "compare_feature_distributions",
    "compare_survey_distributions",
]

_DEFAULT_FEATURES = (
    "log_impedance_magnitude",
    "phase_deg",
    "error_to_magnitude_ratio",
)


def _feature_values(survey: SurveyData, feature: str) -> np.ndarray:
    valid = survey.valid
    z = survey.impedance[valid]
    if feature == "log_impedance_magnitude":
        magnitude = np.abs(z)
        return np.log10(magnitude[magnitude > 0])
    if feature == "phase_deg":
        return np.degrees(np.angle(z))
    if feature == "error_to_magnitude_ratio":
        if survey.impedance_error is None:
            return np.empty(0)
        error = survey.impedance_error[valid]
        magnitude = np.abs(z)
        ratio = error[magnitude > 0] / magnitude[magnitude > 0]
        return ratio[np.isfinite(ratio)]
    raise ValueError(
        f"unknown feature {feature!r}; choose from {_DEFAULT_FEATURES} "
        "or supply values directly."
    )


@dataclass(frozen=True)
class FeatureComparison:
    """Quantitative comparison of one feature between two distributions.

    Parameters
    ----------
    feature : str
        Name of the compared feature.
    simulated_stats, field_stats : mapping
        ``count``, ``mean``, ``std``, ``median`` of each sample.
    ks_statistic, ks_pvalue : float
        Two-sample Kolmogorov-Smirnov statistic and p-value; ``NaN`` when
        either sample is empty.
    mean_difference : float
        ``simulated mean - field mean``.
    std_ratio : float
        ``simulated std / field std``; ``NaN`` when the field std is zero.
    """

    feature: str
    simulated_stats: Mapping[str, float]
    field_stats: Mapping[str, float]
    ks_statistic: float
    ks_pvalue: float
    mean_difference: float
    std_ratio: float

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            All fields, with statistics mappings converted to plain dicts.

        Examples
        --------
        >>> import numpy as np
        >>> comparison = compare_feature_distributions(
        ...     np.array([1.0, 2.0, 3.0]),
        ...     np.array([1.0, 2.0, 3.0]),
        ...     feature="custom",
        ... )
        >>> comparison.to_dict()["feature"]
        'custom'
        """
        return {
            "feature": self.feature,
            "simulated_stats": dict(self.simulated_stats),
            "field_stats": dict(self.field_stats),
            "ks_statistic": self.ks_statistic,
            "ks_pvalue": self.ks_pvalue,
            "mean_difference": self.mean_difference,
            "std_ratio": self.std_ratio,
        }


@dataclass(frozen=True)
class DistributionComparisonReport:
    """Collection of :class:`FeatureComparison` results across features.

    Parameters
    ----------
    comparisons : mapping
        Feature name to :class:`FeatureComparison`.
    """

    comparisons: Mapping[str, FeatureComparison]

    def worst_feature(self) -> str:
        """Return the feature name with the largest KS statistic.

        Returns
        -------
        str
            Feature name whose simulated/field distributions differ most,
            ignoring features with a ``NaN`` statistic (empty samples).

        Raises
        ------
        ValueError
            If every feature has a ``NaN`` KS statistic.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.ai.data.contracts import SurveyData
        >>> z = np.full((2, 4, 1), 100 + 50j)
        >>> survey = SurveyData(
        ...     z, np.linspace(100, 1, 4), ["A", "B"], ["xy"], np.zeros((2, 2))
        ... )
        >>> report = compare_survey_distributions(survey, survey)
        >>> report.worst_feature() in report.comparisons
        True
        """
        finite = {
            name: comparison.ks_statistic
            for name, comparison in self.comparisons.items()
            if np.isfinite(comparison.ks_statistic)
        }
        if not finite:
            raise ValueError("no feature produced a finite KS statistic.")
        return max(finite, key=finite.get)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            Mapping of feature name to its comparison dict.

        Examples
        --------
        >>> import numpy as np
        >>> from pycsamt.ai.data.contracts import SurveyData
        >>> z = np.full((2, 4, 1), 100 + 50j)
        >>> survey = SurveyData(
        ...     z, np.linspace(100, 1, 4), ["A", "B"], ["xy"], np.zeros((2, 2))
        ... )
        >>> report = compare_survey_distributions(survey, survey)
        >>> sorted(report.to_dict())
        ['error_to_magnitude_ratio', 'log_impedance_magnitude', 'phase_deg']
        """
        return {
            name: comparison.to_dict()
            for name, comparison in self.comparisons.items()
        }


def _summary(values: np.ndarray) -> dict[str, float]:
    if values.size == 0:
        return {
            "count": 0,
            "mean": float("nan"),
            "std": float("nan"),
            "median": float("nan"),
        }
    return {
        "count": int(values.size),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
        "median": float(np.median(values)),
    }


def compare_feature_distributions(
    simulated: np.ndarray | Sequence[float],
    field: np.ndarray | Sequence[float],
    *,
    feature: str,
) -> FeatureComparison:
    """Compare two 1-D samples of the same feature quantitatively.

    Parameters
    ----------
    simulated, field : array-like
        Feature values already extracted from each survey (see
        :func:`compare_survey_distributions` to extract them from
        :class:`~pycsamt.ai.data.contracts.SurveyData` directly).
    feature : str
        Label recorded on the returned :class:`FeatureComparison`.

    Returns
    -------
    FeatureComparison
        Summary statistics, mean difference, std ratio, and two-sample KS
        test between the samples.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> comparison = compare_feature_distributions(
    ...     rng.normal(0, 1, 200), rng.normal(0, 1, 200), feature="demo"
    ... )
    >>> comparison.ks_pvalue > 0.01
    True
    """
    simulated = np.asarray(simulated, dtype=float)
    field = np.asarray(field, dtype=float)
    simulated = simulated[np.isfinite(simulated)]
    field = field[np.isfinite(field)]

    ks_statistic = float("nan")
    ks_pvalue = float("nan")
    if simulated.size > 0 and field.size > 0:
        from scipy import stats

        result = stats.ks_2samp(simulated, field)
        ks_statistic = float(result.statistic)
        ks_pvalue = float(result.pvalue)

    sim_stats = _summary(simulated)
    field_stats = _summary(field)
    mean_difference = sim_stats["mean"] - field_stats["mean"]
    field_std = field_stats["std"]
    std_ratio = (
        sim_stats["std"] / field_std
        if np.isfinite(field_std) and field_std != 0.0
        else float("nan")
    )
    return FeatureComparison(
        feature=feature,
        simulated_stats=sim_stats,
        field_stats=field_stats,
        ks_statistic=ks_statistic,
        ks_pvalue=ks_pvalue,
        mean_difference=mean_difference,
        std_ratio=std_ratio,
    )


def compare_survey_distributions(
    simulated: SurveyData,
    field: SurveyData,
    *,
    features: Sequence[str] = _DEFAULT_FEATURES,
) -> DistributionComparisonReport:
    """Compare simulated and field surveys across several canonical features.

    Parameters
    ----------
    simulated : SurveyData
        Simulated (e.g. corrupted synthetic) survey.
    field : SurveyData
        Real field survey, ideally sharing the simulated survey's frequency
        band and component set for a meaningful comparison.
    features : sequence of str, default features
        Any of ``"log_impedance_magnitude"``, ``"phase_deg"``, or
        ``"error_to_magnitude_ratio"``. The error-ratio feature is skipped
        (empty sample) for a survey without declared errors.

    Returns
    -------
    DistributionComparisonReport
        One :class:`FeatureComparison` per requested feature.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((2, 4, 1), 100 + 50j)
    >>> survey = SurveyData(
    ...     z, np.linspace(100, 1, 4), ["A", "B"], ["xy"], np.zeros((2, 2))
    ... )
    >>> report = compare_survey_distributions(survey, survey)
    >>> report.comparisons["phase_deg"].mean_difference
    0.0
    """
    comparisons = {}
    for feature in features:
        sim_values = _feature_values(simulated, feature)
        field_values = _feature_values(field, feature)
        comparisons[feature] = compare_feature_distributions(
            sim_values, field_values, feature=feature
        )
    return DistributionComparisonReport(comparisons=comparisons)

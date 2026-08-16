# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Empirical field-calibrated corruption for synthetic EM surveys."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from ..data.contracts import SurveyData
from ..data.manifest import canonical_hash

__all__ = [
    "EmpiricalCorruptionResult",
    "apply_empirical_corruption",
]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    """Return a read-only array copy."""
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _seed(value: int) -> int:
    """Validate and normalize an unsigned random seed."""
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError("seed must be an integer")
    result = int(value)
    if result < 0 or result >= 2**64:
        raise ValueError("seed must lie in [0, 2**64)")
    return result


def _profiles(
    values: Any,
    field_frequencies_hz: np.ndarray,
    name: str,
) -> np.ndarray:
    """Validate station-by-frequency empirical profiles."""
    array = np.asarray(values, dtype=float)
    expected_frequency = field_frequencies_hz.size
    if (
        array.ndim != 2
        or array.shape[0] < 1
        or array.shape[1] != expected_frequency
        or not np.all(np.isfinite(array))
    ):
        raise ValueError(
            f"{name} must be finite with shape (n_station, "
            f"{expected_frequency})"
        )
    return array


def _bounded_profiles(
    values: Any,
    field_frequencies_hz: np.ndarray,
    name: str,
) -> np.ndarray:
    """Validate empirical reliability profiles in ``[0, 1]``."""
    array = _profiles(values, field_frequencies_hz, name)
    if np.any((array < 0) | (array > 1)):
        raise ValueError(f"{name} must lie in [0, 1]")
    return array


def _interpolate_profiles(
    profiles: np.ndarray,
    source_frequency_hz: np.ndarray,
    target_frequency_hz: np.ndarray,
) -> np.ndarray:
    """Interpolate profiles linearly in log-frequency with edge clamping."""
    order = np.argsort(source_frequency_hz)
    source = np.log10(source_frequency_hz[order])
    target = np.log10(target_frequency_hz)
    return np.stack(
        [np.interp(target, source, row[order]) for row in profiles]
    )


def _quantile_values(
    specification: Mapping[str, Sequence[float]],
    component: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return validated quantile levels and values for one component."""
    if component not in specification:
        raise ValueError(f"relative-error quantiles missing {component!r}")
    value = specification[component]
    if isinstance(value, Mapping):
        levels = np.asarray(value["levels"], dtype=float)
        quantiles = np.asarray(value["values"], dtype=float)
    else:
        raise TypeError("relative-error quantiles must be mappings")
    if (
        levels.ndim != 1
        or quantiles.shape != levels.shape
        or levels.size < 2
        or not np.all(np.isfinite(levels))
        or not np.all(np.isfinite(quantiles))
        or not np.all(np.diff(levels) > 0)
        or not np.all(np.diff(quantiles) >= 0)
        or levels[0] < 0
        or levels[-1] > 1
        or np.any(quantiles < 0)
    ):
        raise ValueError("relative-error quantile curves are invalid")
    return levels, quantiles


@dataclass(frozen=True)
class EmpiricalCorruptionResult:
    """One corrupted survey and every sampled latent array.

    Parameters
    ----------
    survey : SurveyData
        Corrupted survey with declared empirical impedance errors.
    seed : int
        Parent random seed.
    field_station_indices : ndarray of int
        Empirical field station selected for each synthetic station.
    static_log10_resistivity_factor : ndarray
        Sampled ``log10(rho_observed / rho_smooth)``.
    relative_error_fraction : ndarray
        Sampled impedance-error fractions by component.
    noise_realization : ndarray of complex
        Additive complex noise in impedance units.
    missing_mask : ndarray of bool
        True where an observation was removed.
    measurement_reliability, dimensionality_reliability, observation_reliability : ndarray
        Separate and combined reliability factors.
    """

    survey: SurveyData
    seed: int
    field_station_indices: np.ndarray
    static_log10_resistivity_factor: np.ndarray
    relative_error_fraction: np.ndarray
    noise_realization: np.ndarray
    missing_mask: np.ndarray
    measurement_reliability: np.ndarray
    dimensionality_reliability: np.ndarray
    observation_reliability: np.ndarray

    def __post_init__(self) -> None:
        if not isinstance(self.survey, SurveyData):
            raise TypeError("survey must be a SurveyData")
        object.__setattr__(self, "seed", _seed(self.seed))
        shapes = {
            "relative_error_fraction": self.survey.shape,
            "noise_realization": self.survey.shape,
            "missing_mask": self.survey.shape,
            "measurement_reliability": self.survey.shape,
            "dimensionality_reliability": self.survey.shape,
            "observation_reliability": self.survey.shape,
        }
        for name, shape in shapes.items():
            value = np.asarray(getattr(self, name))
            if value.shape != shape:
                raise ValueError(f"{name} must be shaped {shape}")
            object.__setattr__(self, name, _readonly(value))
        static_shape = self.survey.shape[:2]
        static = np.asarray(self.static_log10_resistivity_factor, dtype=float)
        if static.shape != static_shape or not np.all(np.isfinite(static)):
            raise ValueError(
                "static_log10_resistivity_factor must be finite and shaped "
                f"{static_shape}"
            )
        indices = np.asarray(self.field_station_indices, dtype=int)
        if indices.shape != (self.survey.n_stations,) or np.any(indices < 0):
            raise ValueError("field_station_indices have an invalid shape")
        object.__setattr__(
            self,
            "static_log10_resistivity_factor",
            _readonly(static),
        )
        object.__setattr__(
            self, "field_station_indices", _readonly(indices, dtype=int)
        )

    @property
    def record_hash(self) -> str:
        """Return a deterministic hash of seed and sampled summaries.

        Returns
        -------
        str
            SHA-256 digest of JSON-compatible provenance.
        """
        return canonical_hash(self.to_dict())

    def to_dict(self) -> dict[str, Any]:
        """Return compact JSON-compatible corruption provenance.

        Returns
        -------
        dict
            Seed, selected field stations, and sampled ranges.
        """
        return {
            "schema_version": 1,
            "seed": self.seed,
            "field_station_indices": self.field_station_indices.tolist(),
            "static_log10_resistivity_factor_range": [
                float(np.min(self.static_log10_resistivity_factor)),
                float(np.max(self.static_log10_resistivity_factor)),
            ],
            "relative_error_fraction_range": [
                float(np.min(self.relative_error_fraction)),
                float(np.max(self.relative_error_fraction)),
            ],
            "missing_fraction": float(np.mean(self.missing_mask)),
            "measurement_reliability_mean": float(
                np.mean(self.measurement_reliability)
            ),
            "dimensionality_reliability_mean": float(
                np.mean(self.dimensionality_reliability)
            ),
            "observation_reliability_mean": float(
                np.mean(self.observation_reliability)
            ),
        }


def apply_empirical_corruption(
    survey: SurveyData,
    *,
    field_frequencies_hz: Any,
    static_log10_resistivity_profiles: Any,
    measurement_reliability_profiles: Any,
    dimensionality_reliability_profiles: Any,
    observation_reliability_profiles: Any,
    relative_error_quantiles: Mapping[str, Sequence[float]],
    seed: int,
    missing_rate_by_component: Mapping[str, float] | None = None,
) -> EmpiricalCorruptionResult:
    """Apply jointly sampled field profiles and empirical error quantiles.

    Parameters
    ----------
    survey : SurveyData
        Clean synthetic survey.
    field_frequencies_hz : array-like
        Positive unique frequencies of the empirical field profiles.
    static_log10_resistivity_profiles : array-like
        Field ``log10(rho_observed / rho_smooth)`` profiles shaped
        ``(n_field_station, n_field_frequency)``.
    measurement_reliability_profiles, dimensionality_reliability_profiles, observation_reliability_profiles : array-like
        Aligned empirical reliability profiles in ``[0, 1]``.
    relative_error_quantiles : mapping
        Component names mapped to ``{"levels": ..., "values": ...}``
        monotone empirical quantile curves.
    seed : int
        Parent random seed.
    missing_rate_by_component : mapping or None, optional
        Empirical independent missing probabilities. Missing observations are
        set to complex NaN and marked invalid.

    Returns
    -------
    EmpiricalCorruptionResult
        Corrupted survey and complete sampled provenance.

    Notes
    -----
    Static shift is supplied in apparent-resistivity space and therefore
    applied to impedance as ``10**(0.5 * log10_rho_factor)``.

    Examples
    --------
    >>> clean = SurveyData(
    ...     np.ones((2, 2, 1), complex),
    ...     [10.0, 1.0],
    ...     ["A", "B"],
    ...     ["zxy"],
    ...     [[0, 0], [1, 0]],
    ... )
    >>> result = apply_empirical_corruption(
    ...     clean,
    ...     field_frequencies_hz=[10.0, 1.0],
    ...     static_log10_resistivity_profiles=[[0.0, 0.0]],
    ...     measurement_reliability_profiles=[[0.8, 0.7]],
    ...     dimensionality_reliability_profiles=[[1.0, 0.5]],
    ...     observation_reliability_profiles=[[0.8, 0.35]],
    ...     relative_error_quantiles={
    ...         "zxy": {"levels": [0, 1], "values": [0.01, 0.05]}
    ...     },
    ...     seed=0,
    ... )
    >>> result.survey.shape
    (2, 2, 1)
    """
    if not isinstance(survey, SurveyData):
        raise TypeError("survey must be a SurveyData")
    seed = _seed(seed)
    field_frequency = np.asarray(field_frequencies_hz, dtype=float)
    if (
        field_frequency.ndim != 1
        or field_frequency.size < 2
        or not np.all(np.isfinite(field_frequency))
        or np.any(field_frequency <= 0)
        or np.unique(field_frequency).size != field_frequency.size
    ):
        raise ValueError(
            "field_frequencies_hz must be positive, finite, and unique"
        )
    static = _profiles(
        static_log10_resistivity_profiles,
        field_frequency,
        "static_log10_resistivity_profiles",
    )
    measurement = _bounded_profiles(
        measurement_reliability_profiles,
        field_frequency,
        "measurement_reliability_profiles",
    )
    dimensionality = _bounded_profiles(
        dimensionality_reliability_profiles,
        field_frequency,
        "dimensionality_reliability_profiles",
    )
    combined = _bounded_profiles(
        observation_reliability_profiles,
        field_frequency,
        "observation_reliability_profiles",
    )
    if not (
        static.shape
        == measurement.shape
        == dimensionality.shape
        == combined.shape
    ):
        raise ValueError("all empirical profile arrays must have one shape")

    rng = np.random.default_rng(seed)
    indices = rng.integers(0, static.shape[0], size=survey.n_stations)
    target_frequency = survey.frequencies_hz
    static_target = _interpolate_profiles(
        static[indices], field_frequency, target_frequency
    )
    measurement_target = _interpolate_profiles(
        measurement[indices], field_frequency, target_frequency
    )
    dimensionality_target = _interpolate_profiles(
        dimensionality[indices], field_frequency, target_frequency
    )
    combined_target = _interpolate_profiles(
        combined[indices], field_frequency, target_frequency
    )

    impedance_factor = np.power(10.0, 0.5 * static_target)
    shifted = survey.impedance * impedance_factor[:, :, None]
    relative_error = np.empty(survey.shape, dtype=float)
    for component_index, component in enumerate(survey.components):
        levels, values = _quantile_values(
            relative_error_quantiles, component
        )
        probability = rng.uniform(
            levels[0], levels[-1], size=survey.shape[:2]
        )
        relative_error[:, :, component_index] = np.interp(
            probability, levels, values
        )

    declared_error = relative_error * np.abs(shifted)
    standardized_noise = (
        rng.normal(size=survey.shape) + 1j * rng.normal(size=survey.shape)
    ) / np.sqrt(2.0)
    noise = declared_error * standardized_noise
    corrupted = np.where(survey.valid, shifted + noise, survey.impedance)

    missing_rates = {} if missing_rate_by_component is None else dict(
        missing_rate_by_component
    )
    missing = np.zeros(survey.shape, dtype=bool)
    for component_index, component in enumerate(survey.components):
        rate = float(missing_rates.get(component, 0.0))
        if not np.isfinite(rate) or not 0 <= rate <= 1:
            raise ValueError("missing rates must be finite and lie in [0, 1]")
        if rate > 0:
            missing[:, :, component_index] = (
                rng.random(survey.shape[:2]) < rate
            )
    missing |= ~survey.valid
    corrupted = np.where(missing, complex(np.nan, np.nan), corrupted)
    declared_error = np.where(missing, np.nan, declared_error)
    corrupted_survey = replace(
        survey,
        impedance=corrupted,
        impedance_error=declared_error,
    )
    shape = survey.shape
    return EmpiricalCorruptionResult(
        survey=corrupted_survey,
        seed=seed,
        field_station_indices=indices,
        static_log10_resistivity_factor=static_target,
        relative_error_fraction=relative_error,
        noise_realization=noise,
        missing_mask=missing,
        measurement_reliability=np.broadcast_to(
            measurement_target[:, :, None], shape
        ),
        dimensionality_reliability=np.broadcast_to(
            dimensionality_target[:, :, None], shape
        ),
        observation_reliability=np.broadcast_to(
            combined_target[:, :, None], shape
        ),
    )

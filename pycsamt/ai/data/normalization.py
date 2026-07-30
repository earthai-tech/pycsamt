# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Mask-aware, train-fitted normalization for complex MT/AMT surveys.

Normalization state is a scientific artifact.  It is fitted only from
explicitly supplied surveys, records the frequency/component axes and complex
impedance convention, and is then reused unchanged for validation, test, and
field data.  Invalid observations remain identifiable through explicit masks.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from .contracts import ImpedanceConvention, SurveyData
from .manifest import canonical_hash

__all__ = ["NormalizedSurvey", "ComplexZScore"]


def _readonly(value: Any, dtype: Any | None = None) -> np.ndarray:
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _surveys(value: SurveyData | Sequence[SurveyData]) -> list[SurveyData]:
    items = [value] if isinstance(value, SurveyData) else list(value)
    if not items:
        raise ValueError("at least one training survey is required.")
    if any(not isinstance(item, SurveyData) for item in items):
        raise TypeError("every training item must be a SurveyData instance.")
    return items


@dataclass(frozen=True)
class NormalizedSurvey:
    """Normalized real/imaginary feature channels and their validity state.

    Parameters
    ----------
    values : ndarray, shape (n_station, n_frequency, n_component, 2)
        Normalized channels ordered as real then imaginary.
    valid : ndarray of bool
        Mask with the same shape as ``values``. Both channels of an impedance
        observation normally share one validity state.
    frequencies_hz : ndarray, shape (n_frequency,)
        Frequency axis used by the fitted normalizer.
    station_names, components : sequence of str
        Axis labels corresponding to ``values``.
    errors : ndarray or None, optional
        Normalized absolute standard errors with the same shape as ``values``.
    state_hash : str or None, optional
        Digest of the normalization state that produced the features.

    Examples
    --------
    ``NormalizedSurvey`` objects are normally created with
    :meth:`ComplexZScore.transform_survey`:

    >>> values = np.zeros((1, 2, 1, 2))
    >>> result = NormalizedSurvey(
    ...     values, np.ones_like(values, bool), [10, 1], ["S"], ["xy"]
    ... )
    >>> result.shape
    (1, 2, 1, 2)
    >>> result.n_valid_observations
    2
    """

    values: np.ndarray
    valid: np.ndarray
    frequencies_hz: np.ndarray
    station_names: tuple[str, ...]
    components: tuple[str, ...]
    errors: np.ndarray | None = None
    state_hash: str | None = None

    def __post_init__(self) -> None:
        values = np.asarray(self.values, dtype=float)
        valid = np.asarray(self.valid, dtype=bool)
        if values.ndim != 4 or values.shape[-1] != 2:
            raise ValueError(
                "values must have shape (station, frequency, component, 2)."
            )
        if valid.shape != values.shape:
            raise ValueError("valid must have the same shape as values.")
        if not np.all(np.isfinite(values)):
            raise ValueError(
                "normalized values must be finite; use a finite fill value."
            )
        n_station, n_frequency, n_component, _ = values.shape
        frequencies = np.asarray(self.frequencies_hz, dtype=float)
        if frequencies.shape != (n_frequency,) or not np.all(
            np.isfinite(frequencies)
        ):
            raise ValueError(
                "frequencies_hz does not match the normalized frequency axis."
            )
        station_names = tuple(str(name) for name in self.station_names)
        components = tuple(str(name) for name in self.components)
        if (
            len(station_names) != n_station
            or len(set(station_names)) != n_station
        ):
            raise ValueError(
                "station_names must uniquely label the station axis."
            )
        if (
            len(components) != n_component
            or len(set(components)) != n_component
        ):
            raise ValueError(
                "components must uniquely label the component axis."
            )
        errors = None
        if self.errors is not None:
            errors = np.asarray(self.errors, dtype=float)
            if errors.shape != values.shape:
                raise ValueError("errors must have the same shape as values.")
            if np.any(valid & (~np.isfinite(errors) | (errors <= 0))):
                raise ValueError(
                    "valid normalized errors must be finite and positive."
                )
        state_hash = None if self.state_hash is None else str(self.state_hash)
        if state_hash is not None and (
            len(state_hash) != 64
            or any(
                character not in "0123456789abcdef" for character in state_hash
            )
        ):
            raise ValueError("state_hash must be a lowercase SHA-256 digest.")
        object.__setattr__(self, "values", _readonly(values))
        object.__setattr__(self, "valid", _readonly(valid, bool))
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "station_names", station_names)
        object.__setattr__(self, "components", components)
        object.__setattr__(
            self, "errors", None if errors is None else _readonly(errors)
        )
        object.__setattr__(self, "state_hash", state_hash)

    @property
    def shape(self) -> tuple[int, int, int, int]:
        """Return the normalized feature shape.

        Returns
        -------
        tuple of int
            ``(n_station, n_frequency, n_component, 2)``.

        Examples
        --------
        >>> x = np.zeros((2, 3, 1, 2))
        >>> n = NormalizedSurvey(
        ...     x, np.ones_like(x, bool), [100, 10, 1], ["A", "B"], ["xy"]
        ... )
        >>> n.shape
        (2, 3, 1, 2)
        """
        return self.values.shape

    @property
    def n_valid_observations(self) -> int:
        """Return the number of valid complex observations.

        Returns
        -------
        int
            Count on the impedance grid, not the doubled channel count.

        Examples
        --------
        >>> x = np.zeros((1, 1, 1, 2))
        >>> n = NormalizedSurvey(x, np.ones_like(x, bool), [1], ["S"], ["xy"])
        >>> n.n_valid_observations
        1
        """
        return int(np.count_nonzero(np.all(self.valid, axis=-1)))

    def flatten(self) -> tuple[np.ndarray, np.ndarray]:
        """Flatten frequency, component, and channel axes for dense models.

        Returns
        -------
        values : ndarray, shape (n_station, n_feature)
            Read-only station-major feature matrix.
        valid : ndarray of bool, shape (n_station, n_feature)
            Read-only feature mask in identical order.

        Examples
        --------
        >>> x = np.zeros((2, 3, 2, 2))
        >>> n = NormalizedSurvey(
        ...     x,
        ...     np.ones_like(x, bool),
        ...     [100, 10, 1],
        ...     ["A", "B"],
        ...     ["xy", "yx"],
        ... )
        >>> n.flatten()[0].shape
        (2, 12)
        """
        return (
            _readonly(self.values.reshape(self.shape[0], -1)),
            _readonly(self.valid.reshape(self.shape[0], -1), bool),
        )


@dataclass(frozen=True)
class ComplexZScore:
    """Immutable per-frequency/component complex z-score state.

    Real and imaginary impedance parts are standardized independently. The
    stored statistic shape is ``(frequency, component, channel)`` where the
    final channels are real and imaginary.

    Parameters
    ----------
    mean, scale : ndarray
        Per-feature location and positive scale arrays.
    frequencies_hz : ndarray
        Exact fitted frequency grid and order.
    components : sequence of str
        Exact fitted component order.
    eps : float, default=1e-8
        Minimum allowed scale.
    count : ndarray or None, optional
        Number of valid training observations supporting each channel.
    weight_sum : ndarray or None, optional
        Sum of fitting weights supporting each channel.
    weighting : {"uniform", "inverse_variance"}, default="uniform"
        Statistic weighting policy.
    ddof : int, default=0
        Delta degrees of freedom used for uniform variance.
    convention : ImpedanceConvention or None, optional
        Complex convention bound to the fitted state. ``None`` is accepted only
        for legacy schema-1 states and skips convention compatibility checks.
    training_survey_count, training_station_count : int, optional
        Audit counts describing the data used during fitting.

    Examples
    --------
    Fit on training data and reuse the same state for later data:

    >>> z = np.array([[[1 + 2j]], [[3 + 4j]]])
    >>> training = SurveyData(z, [1], ["A", "B"], ["xy"], [[0, 0], [1, 0]])
    >>> normalizer = ComplexZScore.fit(training)
    >>> features, mask = normalizer.transform(training)
    >>> np.allclose(features[:, 0, 0, 0], [-1, 1])
    True
    >>> mask.all()
    True
    """

    mean: np.ndarray
    scale: np.ndarray
    frequencies_hz: np.ndarray
    components: tuple[str, ...]
    eps: float = 1e-8
    count: np.ndarray | None = None
    weight_sum: np.ndarray | None = None
    weighting: str = "uniform"
    ddof: int = 0
    convention: ImpedanceConvention | None = None
    training_survey_count: int | None = None
    training_station_count: int | None = None

    def __post_init__(self) -> None:
        mean = np.asarray(self.mean, dtype=float)
        scale = np.asarray(self.scale, dtype=float)
        frequencies = np.asarray(self.frequencies_hz, dtype=float)
        components = tuple(str(name).strip() for name in self.components)
        expected = (len(frequencies), len(components), 2)
        if mean.shape != expected or scale.shape != expected:
            raise ValueError(f"mean and scale must have shape {expected}.")
        if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(scale)):
            raise ValueError("normalization statistics must be finite.")
        eps = float(self.eps)
        if np.any(scale <= 0) or not np.isfinite(eps) or eps <= 0:
            raise ValueError("scale and eps must be finite and positive.")
        if (
            frequencies.ndim != 1
            or not np.all(np.isfinite(frequencies))
            or np.any(frequencies <= 0)
        ):
            raise ValueError(
                "frequencies_hz must be a finite positive 1-D array."
            )
        differences = np.diff(frequencies)
        if differences.size and not (
            np.all(differences > 0) or np.all(differences < 0)
        ):
            raise ValueError(
                "frequencies_hz must be strictly monotonic and unique."
            )
        if (
            not components
            or any(not name for name in components)
            or len(set(components)) != len(components)
        ):
            raise ValueError("components must be non-empty and unique.")
        if self.weighting not in {"uniform", "inverse_variance"}:
            raise ValueError(
                "weighting must be 'uniform' or 'inverse_variance'."
            )
        if (
            not isinstance(self.ddof, int)
            or isinstance(self.ddof, bool)
            or self.ddof < 0
        ):
            raise ValueError("ddof must be a non-negative integer.")
        if self.weighting != "uniform" and self.ddof != 0:
            raise ValueError(
                "ddof must be zero for inverse-variance weighting."
            )
        count = (
            np.ones(expected, dtype=np.int64)
            if self.count is None
            else np.asarray(self.count)
        )
        if (
            count.shape != expected
            or not np.issubdtype(count.dtype, np.integer)
            or np.any(count <= self.ddof)
        ):
            raise ValueError(
                f"count must be an integer array shaped {expected} with values greater than ddof."
            )
        weight_sum = (
            count.astype(float)
            if self.weight_sum is None
            else np.asarray(self.weight_sum, dtype=float)
        )
        if (
            weight_sum.shape != expected
            or not np.all(np.isfinite(weight_sum))
            or np.any(weight_sum <= 0)
        ):
            raise ValueError(
                f"weight_sum must be finite, positive, and shaped {expected}."
            )
        convention = self.convention
        if isinstance(convention, Mapping):
            convention = ImpedanceConvention.from_dict(convention)
        if convention is not None and not isinstance(
            convention, ImpedanceConvention
        ):
            raise TypeError(
                "convention must be an ImpedanceConvention or None."
            )
        for name, value in (
            ("training_survey_count", self.training_survey_count),
            ("training_station_count", self.training_station_count),
        ):
            if value is not None and (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value <= 0
            ):
                raise ValueError(f"{name} must be a positive integer or None.")
        if (
            self.training_station_count is not None
            and int(np.max(count)) > self.training_station_count
        ):
            raise ValueError(
                "training_station_count cannot be smaller than feature counts."
            )
        for name, value, dtype in (
            ("mean", mean, float),
            ("scale", scale, float),
            ("frequencies_hz", frequencies, float),
            ("count", count, np.int64),
            ("weight_sum", weight_sum, float),
        ):
            object.__setattr__(self, name, _readonly(value, dtype))
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "eps", eps)
        object.__setattr__(self, "convention", convention)

    @classmethod
    def fit(
        cls,
        surveys: SurveyData | Sequence[SurveyData],
        *,
        eps: float = 1e-8,
        weighting: str = "uniform",
        ddof: int = 0,
    ) -> ComplexZScore:
        """Fit statistics from explicitly supplied training surveys.

        Parameters
        ----------
        surveys : SurveyData or sequence of SurveyData
            Training surveys sharing frequency grid/order, component order,
            and impedance convention. Stations are pooled across surveys.
        eps : float, default=1e-8
            Positive scale floor for constant features.
        weighting : {"uniform", "inverse_variance"}, default="uniform"
            Use equal weights or inverse squared impedance errors. The latter
            requires error arrays on every training survey.
        ddof : int, default=0
            Delta degrees of freedom for uniform variance. It must be smaller
            than every feature's valid observation count and must be zero for
            inverse-variance weighting.

        Returns
        -------
        ComplexZScore
            Immutable fitted state for reuse on validation/test/field data.

        Raises
        ------
        ValueError
            If axes or conventions differ, a feature has insufficient valid
            observations, or requested weights are unavailable.
        TypeError
            If an input is not :class:`SurveyData`.

        Examples
        --------
        >>> z = np.array([[[1 + 1j]], [[2 + 3j]], [[5 + 7j]]])
        >>> survey = SurveyData(
        ...     z, [1], ["a", "b", "c"], ["xy"], [[0, 0], [1, 0], [2, 0]]
        ... )
        >>> state = ComplexZScore.fit(survey, ddof=1)
        >>> state.training_station_count
        3
        >>> state.count[0, 0, 0]
        3
        """
        items = _surveys(surveys)
        if weighting not in {"uniform", "inverse_variance"}:
            raise ValueError(
                "weighting must be 'uniform' or 'inverse_variance'."
            )
        if not isinstance(ddof, int) or isinstance(ddof, bool) or ddof < 0:
            raise ValueError("ddof must be a non-negative integer.")
        if weighting != "uniform" and ddof != 0:
            raise ValueError(
                "ddof must be zero for inverse-variance weighting."
            )
        reference = items[0]
        values = []
        masks = []
        weights = []
        for survey in items:
            reference.assert_compatible(survey, require_crs=False)
            values.append(
                np.stack(
                    [survey.impedance.real, survey.impedance.imag], axis=-1
                )
            )
            mask = np.repeat(survey.valid[..., None], 2, axis=-1)
            masks.append(mask)
            if weighting == "inverse_variance":
                if survey.impedance_error is None:
                    raise ValueError(
                        "inverse-variance weighting requires impedance_error on every survey."
                    )
                error = np.repeat(
                    survey.impedance_error[..., None], 2, axis=-1
                )
                weights.append(np.where(mask, 1.0 / np.square(error), 0.0))
            else:
                weights.append(mask.astype(float))
        x = np.concatenate(values, axis=0)
        mask = np.concatenate(masks, axis=0)
        weight = np.concatenate(weights, axis=0)
        count = mask.sum(axis=0)
        if np.any(count <= ddof):
            raise ValueError(
                "every feature needs more valid training observations than ddof."
            )
        weight_sum = weight.sum(axis=0)
        if np.any(~np.isfinite(weight_sum)) or np.any(weight_sum <= 0):
            raise ValueError(
                "training weights must have a positive finite sum for every feature."
            )
        safe_x = np.where(mask, x, 0.0)
        mean = (weight * safe_x).sum(axis=0) / weight_sum
        squared = np.where(mask, np.square(x - mean), 0.0)
        if weighting == "uniform":
            variance = (weight * squared).sum(axis=0) / (count - ddof)
        else:
            variance = (weight * squared).sum(axis=0) / weight_sum
        scale = np.maximum(np.sqrt(variance), float(eps))
        return cls(
            mean=mean,
            scale=scale,
            frequencies_hz=reference.frequencies_hz,
            components=reference.components,
            eps=eps,
            count=count,
            weight_sum=weight_sum,
            weighting=weighting,
            ddof=ddof,
            convention=reference.convention,
            training_survey_count=len(items),
            training_station_count=sum(item.n_stations for item in items),
        )

    @property
    def state_hash(self) -> str:
        """Return a deterministic digest of the complete fitted state.

        Returns
        -------
        str
            Lowercase SHA-256 digest suitable for artifact provenance.

        Examples
        --------
        >>> state = ComplexZScore(
        ...     np.zeros((1, 1, 2)), np.ones((1, 1, 2)), [1], ["xy"]
        ... )
        >>> len(state.state_hash)
        64
        """
        return canonical_hash(self.to_dict())

    @property
    def feature_names(self) -> tuple[str, ...]:
        """Return flattened feature names in canonical array order.

        Returns
        -------
        tuple of str
            Names ordered by frequency, component, then real/imaginary channel.

        Examples
        --------
        >>> state = ComplexZScore(
        ...     np.zeros((1, 1, 2)), np.ones((1, 1, 2)), [10], ["xy"]
        ... )
        >>> state.feature_names
        ('10Hz:xy:real', '10Hz:xy:imag')
        """
        return tuple(
            f"{frequency:g}Hz:{component}:{channel}"
            for frequency in self.frequencies_hz
            for component in self.components
            for channel in ("real", "imag")
        )

    def validate_survey(self, survey: SurveyData) -> None:
        """Validate survey axes and complex convention against fitted state.

        Parameters
        ----------
        survey : SurveyData
            Candidate survey for transformation.

        Returns
        -------
        None
            Successful return means the survey is transform-compatible.

        Raises
        ------
        TypeError
            If ``survey`` is not :class:`SurveyData`.
        ValueError
            If frequency values/order, component order, or a recorded complex
            convention differs.

        Examples
        --------
        >>> survey = SurveyData(
        ...     np.ones((1, 1, 1), complex), [1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> state = ComplexZScore.fit(survey)
        >>> state.validate_survey(survey) is None
        True
        """
        if not isinstance(survey, SurveyData):
            raise TypeError("survey must be a SurveyData instance.")
        if not np.array_equal(survey.frequencies_hz, self.frequencies_hz):
            raise ValueError(
                "survey frequency grid/order differs from fitted state."
            )
        if survey.components != self.components:
            raise ValueError(
                "survey component order differs from fitted state."
            )
        if (
            self.convention is not None
            and survey.convention != self.convention
        ):
            raise ValueError(
                "survey impedance convention differs from fitted state."
            )

    def transform(
        self,
        survey: SurveyData,
        *,
        fill_value: float = 0.0,
        clip: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Normalize impedance and return features with an explicit mask.

        Parameters
        ----------
        survey : SurveyData
            Compatible survey to transform.
        fill_value : float, default=0.0
            Finite value placed in invalid feature channels.
        clip : float or None, optional
            Symmetric positive z-score limit. Clipping can improve numerical
            robustness but makes exact inversion impossible for clipped values.

        Returns
        -------
        features : ndarray
            Shape ``(station, frequency, component, 2)`` with real then
            imaginary channels.
        valid : ndarray of bool
            Same shape, retaining the observation mask independently of fill.

        Examples
        --------
        >>> z = np.array([[[1 + 2j]], [[3 + 4j]]])
        >>> survey = SurveyData(z, [1], ["A", "B"], ["xy"], [[0, 0], [1, 0]])
        >>> features, valid = ComplexZScore.fit(survey).transform(survey)
        >>> features.shape, valid.all()
        ((2, 1, 1, 2), True)
        """
        self.validate_survey(survey)
        fill = float(fill_value)
        if not np.isfinite(fill):
            raise ValueError("fill_value must be finite.")
        if clip is not None and (not np.isfinite(clip) or clip <= 0):
            raise ValueError("clip must be finite and positive or None.")
        x = np.stack([survey.impedance.real, survey.impedance.imag], axis=-1)
        transformed = (x - self.mean) / self.scale
        if clip is not None:
            transformed = np.clip(transformed, -float(clip), float(clip))
        mask = np.repeat(survey.valid[..., None], 2, axis=-1)
        return np.where(mask, transformed, fill), mask

    def transform_errors(
        self, survey: SurveyData, *, fill_value: float = 1.0
    ) -> tuple[np.ndarray, np.ndarray]:
        """Propagate absolute impedance errors into normalized channel units.

        Parameters
        ----------
        survey : SurveyData
            Compatible survey containing ``impedance_error``.
        fill_value : float, default=1.0
            Positive finite error assigned to invalid channels.

        Returns
        -------
        errors : ndarray
            Normalized errors shaped ``(station, frequency, component, 2)``.
            The same scalar complex-impedance error is divided by the separate
            real and imaginary channel scales.
        valid : ndarray of bool
            Expanded observation mask.

        Raises
        ------
        ValueError
            If the survey has no impedance errors or fill is invalid.

        Examples
        --------
        >>> z = np.array([[[1 + 1j]], [[3 + 3j]]])
        >>> s = SurveyData(
        ...     z,
        ...     [1],
        ...     ["A", "B"],
        ...     ["xy"],
        ...     [[0, 0], [1, 0]],
        ...     impedance_error=np.ones_like(z.real),
        ... )
        >>> errors, valid = ComplexZScore.fit(s).transform_errors(s)
        >>> errors.shape, valid.all()
        ((2, 1, 1, 2), True)
        """
        self.validate_survey(survey)
        if survey.impedance_error is None:
            raise ValueError("survey does not contain impedance_error.")
        fill = float(fill_value)
        if not np.isfinite(fill) or fill <= 0:
            raise ValueError("fill_value must be finite and positive.")
        error = (
            np.repeat(survey.impedance_error[..., None], 2, axis=-1)
            / self.scale
        )
        mask = np.repeat(survey.valid[..., None], 2, axis=-1)
        return np.where(mask, error, fill), mask

    def transform_survey(
        self,
        survey: SurveyData,
        *,
        fill_value: float = 0.0,
        error_fill_value: float = 1.0,
        clip: float | None = None,
    ) -> NormalizedSurvey:
        """Create a labeled immutable normalized-survey container.

        Parameters
        ----------
        survey : SurveyData
            Compatible survey.
        fill_value, error_fill_value : float, optional
            Finite values for invalid features and invalid normalized errors.
        clip : float or None, optional
            Symmetric feature clipping threshold.

        Returns
        -------
        NormalizedSurvey
            Features, mask, optional propagated errors, axes, and state digest.

        Examples
        --------
        >>> z = np.array([[[1 + 2j]], [[3 + 4j]]])
        >>> s = SurveyData(z, [1], ["A", "B"], ["xy"], [[0, 0], [1, 0]])
        >>> result = ComplexZScore.fit(s).transform_survey(s)
        >>> result.station_names
        ('A', 'B')
        """
        features, mask = self.transform(
            survey, fill_value=fill_value, clip=clip
        )
        errors = None
        if survey.impedance_error is not None:
            errors, _ = self.transform_errors(
                survey, fill_value=error_fill_value
            )
        return NormalizedSurvey(
            values=features,
            valid=mask,
            frequencies_hz=survey.frequencies_hz,
            station_names=survey.station_names,
            components=survey.components,
            errors=errors,
            state_hash=self.state_hash,
        )

    def inverse_transform(
        self,
        values: np.ndarray,
        *,
        valid: np.ndarray | None = None,
        invalid_value: complex = complex(np.nan, np.nan),
    ) -> np.ndarray:
        """Convert normalized channels back to complex SI impedance.

        Parameters
        ----------
        values : ndarray
            Shape ``(station, frequency, component, 2)``.
        valid : ndarray of bool, optional
            Mask with the same shape. An impedance is retained only if both
            real and imaginary channels are valid.
        invalid_value : complex, default=nan+nanj
            Value assigned where either channel is invalid.

        Returns
        -------
        ndarray of complex
            Shape ``(station, frequency, component)`` in V/A.

        Examples
        --------
        >>> z = np.array([[[1 + 2j]], [[3 + 4j]]])
        >>> s = SurveyData(z, [1], ["A", "B"], ["xy"], [[0, 0], [1, 0]])
        >>> state = ComplexZScore.fit(s)
        >>> features, mask = state.transform(s)
        >>> np.allclose(state.inverse_transform(features, valid=mask), z)
        True
        """
        x = np.asarray(values, dtype=float)
        if x.ndim != 4 or x.shape[1:] != self.mean.shape:
            raise ValueError(
                f"values must have shape (station, {self.mean.shape[0]}, {self.mean.shape[1]}, 2)."
            )
        if not np.all(np.isfinite(x)):
            raise ValueError("normalized values must be finite.")
        raw = x * self.scale + self.mean
        complex_values = raw[..., 0] + 1j * raw[..., 1]
        if valid is not None:
            mask = np.asarray(valid, dtype=bool)
            if mask.shape != x.shape:
                raise ValueError("valid must have the same shape as values.")
            try:
                invalid = complex(invalid_value)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "invalid_value must be complex-compatible."
                ) from exc
            complex_values = np.where(
                np.all(mask, axis=-1), complex_values, invalid
            )
        return complex_values

    def to_dict(self) -> dict[str, Any]:
        """Return a complete JSON-serializable schema-2 state.

        Returns
        -------
        dict
            Fitted statistics, axes, counts, fitting policy, convention, and
            training audit counts.

        Examples
        --------
        >>> state = ComplexZScore(
        ...     np.zeros((1, 1, 2)), np.ones((1, 1, 2)), [1], ["xy"]
        ... )
        >>> state.to_dict()["schema_version"]
        2
        """
        return {
            "schema_version": 2,
            "kind": "complex_cartesian_zscore",
            "mean": self.mean.tolist(),
            "scale": self.scale.tolist(),
            "frequencies_hz": self.frequencies_hz.tolist(),
            "components": list(self.components),
            "eps": self.eps,
            "count": self.count.tolist(),
            "weight_sum": self.weight_sum.tolist(),
            "weighting": self.weighting,
            "ddof": self.ddof,
            "convention": None
            if self.convention is None
            else self.convention.to_dict(),
            "training_survey_count": self.training_survey_count,
            "training_station_count": self.training_station_count,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ComplexZScore:
        """Restore a schema-1 or schema-2 fitted state.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`, or the earlier
            schema-1 Cartesian z-score representation.

        Returns
        -------
        ComplexZScore
            Validated immutable runtime state.

        Raises
        ------
        ValueError
            If the schema discriminator is unsupported or statistics violate
            the normalization contract.

        Examples
        --------
        >>> original = ComplexZScore(
        ...     np.zeros((1, 1, 2)), np.ones((1, 1, 2)), [1], ["xy"]
        ... )
        >>> restored = ComplexZScore.from_dict(original.to_dict())
        >>> restored.state_hash == original.state_hash
        True
        """
        version = data.get("schema_version")
        if (
            version not in {1, 2}
            or data.get("kind") != "complex_cartesian_zscore"
        ):
            raise ValueError("unsupported ComplexZScore state.")
        if version == 1:
            mean = np.asarray(data["mean"], dtype=float)
            count = np.ones(mean.shape, dtype=np.int64)
            return cls(
                mean=mean,
                scale=np.asarray(data["scale"], dtype=float),
                frequencies_hz=np.asarray(data["frequencies_hz"], dtype=float),
                components=tuple(data["components"]),
                eps=float(data["eps"]),
                count=count,
                weight_sum=count.astype(float),
                convention=None,
            )
        convention_data = data.get("convention")
        return cls(
            mean=np.asarray(data["mean"], dtype=float),
            scale=np.asarray(data["scale"], dtype=float),
            frequencies_hz=np.asarray(data["frequencies_hz"], dtype=float),
            components=tuple(data["components"]),
            eps=float(data["eps"]),
            count=np.asarray(data["count"], dtype=np.int64),
            weight_sum=np.asarray(data["weight_sum"], dtype=float),
            weighting=data["weighting"],
            ddof=int(data["ddof"]),
            convention=None
            if convention_data is None
            else ImpedanceConvention.from_dict(convention_data),
            training_survey_count=data.get("training_survey_count"),
            training_station_count=data.get("training_station_count"),
        )

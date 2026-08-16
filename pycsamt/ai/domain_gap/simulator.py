# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reproducible corruption of canonical surveys for domain-gap training.

Every function in this module maps a validated
:class:`~pycsamt.ai.data.contracts.SurveyData` to a new, independently
validated :class:`~pycsamt.ai.data.contracts.SurveyData`. Corruptions are
composed from a single integer seed so a training run can be reproduced
exactly from its :class:`CorruptionConfig` and seed alone.

Two families of corruption are distinguished:

Systematic (bias-like)
    :func:`apply_static_shift` and :func:`apply_galvanic_distortion` model
    near-surface, frequency-independent effects that multiply the true
    impedance by a real matrix. They change the signal, not its declared
    uncertainty, beyond the linear scaling of that uncertainty.

Random (noise-like)
    :func:`add_heteroscedastic_noise`, :func:`apply_error_floor`,
    :func:`apply_dropout`, and :func:`inject_outliers` model acquisition
    noise, missing data, and undetected bad readings. These update
    ``impedance_error`` and ``valid`` so that masks and errors stay
    consistent with the injected corruption, per the M3 acceptance gate.

:func:`apply_corruption_suite` composes all of the above in a fixed,
documented order using seeds spawned from a single parent seed, and returns
a :class:`CorruptionRecord` describing exactly what was sampled.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from ..data.contracts import SurveyData
from ..data.manifest import canonical_hash

__all__ = [
    "CorruptionConfig",
    "CorruptionRecord",
    "SEVERITY_PRESETS",
    "add_heteroscedastic_noise",
    "apply_error_floor",
    "apply_static_shift",
    "apply_galvanic_distortion",
    "apply_dropout",
    "inject_outliers",
    "perturb_coordinates",
    "apply_corruption_suite",
]

_COMPONENT_INDEX: dict[str, tuple[int, int]] = {
    "xx": (0, 0),
    "xy": (0, 1),
    "yx": (1, 0),
    "yy": (1, 1),
}


def _seed(value: int) -> int:
    if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
        raise TypeError("seed must be an integer.")
    result = int(value)
    if result < 0 or result >= 2**64:
        raise ValueError("seed must be in [0, 2**64).")
    return result


def _rate(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or not (0.0 <= result <= 1.0):
        raise ValueError(f"{name} must be finite and in [0, 1].")
    return result


def _non_negative(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return result


def _range(
    value: tuple[float, float], name: str, *, allow_negative: bool = False
) -> tuple[float, float]:
    lo, hi = float(value[0]), float(value[1])
    if not (np.isfinite(lo) and np.isfinite(hi)):
        raise ValueError(f"{name} bounds must be finite.")
    if not allow_negative and (lo < 0.0 or hi < 0.0):
        raise ValueError(f"{name} bounds must be non-negative.")
    if lo > hi:
        raise ValueError(f"{name} lower bound cannot exceed the upper bound.")
    return (lo, hi)


def _child_rng(rng: np.random.Generator) -> np.random.Generator:
    return np.random.default_rng(rng.integers(0, 2**63 - 1))


def _present_components(survey: SurveyData) -> dict[str, int]:
    return {
        name: survey.components.index(name)
        for name in _COMPONENT_INDEX
        if name in survey.components
    }


def _z_matrix(survey: SurveyData) -> np.ndarray:
    """Return a dense (station, frequency, 2, 2) view; absent parts are 0."""
    present = _present_components(survey)
    out = np.zeros(survey.shape[:2] + (2, 2), dtype=complex)
    for name, col in present.items():
        i, j = _COMPONENT_INDEX[name]
        out[:, :, i, j] = survey.impedance[:, :, col]
    return out, present


@dataclass(frozen=True)
class CorruptionConfig:
    """Parameter ranges for one corruption pass over a :class:`SurveyData`.

    All defaults are zero/no-op so ``CorruptionConfig()`` is the clean
    synthetic control set required by the M3 gate.

    Parameters
    ----------
    noise_level_range : (float, float), default=(0.0, 0.0)
        Bounds on the relative heteroscedastic noise standard deviation
        sampled independently per station/frequency observation.
    error_floor_fraction : float, default=0.0
        Minimum declared ``impedance_error`` as a fraction of ``|Z|``.
    static_shift_log10_sigma : float, default=0.0
        Std. dev. of the log\\ :sub:`10` per-station static-shift factor
        applied identically across all frequencies.
    distortion_gain_log10_sigma, distortion_twist_deg_sigma, distortion_shear_sigma, distortion_anisotropy_sigma : float, default=0.0
        Std. dev. of the per-station Groom-Bailey-style gain, twist, shear,
        and anisotropy parameters of the injected galvanic distortion.
    station_dropout_rate, frequency_dropout_rate, random_dropout_rate : float, default=0.0
        Probability that an entire station, an entire frequency (across all
        stations), or an individual observation is marked missing.
    outlier_rate : float, default=0.0
        Fraction of remaining valid observations perturbed by a large,
        undetected multiplicative shift.
    outlier_log10_shift_range : (float, float), default=(0.5, 1.5)
        Bounds on the magnitude (in log\\ :sub:`10` decades) of injected
        outliers; the sign is randomized.
    coordinate_sigma_m, elevation_sigma_m : float, default=0.0
        Std. dev. of Gaussian perturbation applied to station horizontal
        coordinates and elevation, respectively.

    Examples
    --------
    >>> config = CorruptionConfig(noise_level_range=(0.01, 0.05))
    >>> config.config_hash() == CorruptionConfig(
    ...     noise_level_range=(0.01, 0.05)
    ... ).config_hash()
    True
    """

    noise_level_range: tuple[float, float] = (0.0, 0.0)
    error_floor_fraction: float = 0.0
    static_shift_log10_sigma: float = 0.0
    distortion_gain_log10_sigma: float = 0.0
    distortion_twist_deg_sigma: float = 0.0
    distortion_shear_sigma: float = 0.0
    distortion_anisotropy_sigma: float = 0.0
    station_dropout_rate: float = 0.0
    frequency_dropout_rate: float = 0.0
    random_dropout_rate: float = 0.0
    outlier_rate: float = 0.0
    outlier_log10_shift_range: tuple[float, float] = (0.5, 1.5)
    coordinate_sigma_m: float = 0.0
    elevation_sigma_m: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "noise_level_range",
            _range(self.noise_level_range, "noise_level_range"),
        )
        object.__setattr__(
            self,
            "error_floor_fraction",
            _non_negative(self.error_floor_fraction, "error_floor_fraction"),
        )
        for name in (
            "static_shift_log10_sigma",
            "distortion_gain_log10_sigma",
            "distortion_twist_deg_sigma",
            "distortion_shear_sigma",
            "distortion_anisotropy_sigma",
            "coordinate_sigma_m",
            "elevation_sigma_m",
        ):
            object.__setattr__(
                self, name, _non_negative(getattr(self, name), name)
            )
        for name in (
            "station_dropout_rate",
            "frequency_dropout_rate",
            "random_dropout_rate",
            "outlier_rate",
        ):
            object.__setattr__(self, name, _rate(getattr(self, name), name))
        object.__setattr__(
            self,
            "outlier_log10_shift_range",
            _range(
                self.outlier_log10_shift_range, "outlier_log10_shift_range"
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable, order-stable representation.

        Returns
        -------
        dict
            Field values with a schema discriminator.

        Examples
        --------
        >>> CorruptionConfig().to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "noise_level_range": list(self.noise_level_range),
            "error_floor_fraction": self.error_floor_fraction,
            "static_shift_log10_sigma": self.static_shift_log10_sigma,
            "distortion_gain_log10_sigma": self.distortion_gain_log10_sigma,
            "distortion_twist_deg_sigma": self.distortion_twist_deg_sigma,
            "distortion_shear_sigma": self.distortion_shear_sigma,
            "distortion_anisotropy_sigma": self.distortion_anisotropy_sigma,
            "station_dropout_rate": self.station_dropout_rate,
            "frequency_dropout_rate": self.frequency_dropout_rate,
            "random_dropout_rate": self.random_dropout_rate,
            "outlier_rate": self.outlier_rate,
            "outlier_log10_shift_range": list(self.outlier_log10_shift_range),
            "coordinate_sigma_m": self.coordinate_sigma_m,
            "elevation_sigma_m": self.elevation_sigma_m,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CorruptionConfig:
        """Restore and validate a serialized configuration.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`.

        Returns
        -------
        CorruptionConfig
            Validated immutable configuration.

        Examples
        --------
        >>> state = CorruptionConfig(error_floor_fraction=0.02).to_dict()
        >>> CorruptionConfig.from_dict(state).error_floor_fraction
        0.02
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported CorruptionConfig schema version.")
        payload = dict(data)
        payload.pop("schema_version")
        payload["noise_level_range"] = tuple(payload["noise_level_range"])
        payload["outlier_log10_shift_range"] = tuple(
            payload["outlier_log10_shift_range"]
        )
        return cls(**payload)

    def config_hash(self) -> str:
        """Return the SHA-256 digest of this configuration's canonical JSON.

        Returns
        -------
        str
            Lowercase 64-character hexadecimal digest.

        Examples
        --------
        >>> len(CorruptionConfig().config_hash())
        64
        """
        return canonical_hash(self.to_dict())


@dataclass(frozen=True)
class CorruptionRecord:
    """Provenance of one applied :func:`apply_corruption_suite` call.

    Parameters
    ----------
    config : CorruptionConfig
        Configuration the sampled parameters were drawn from.
    seed : int
        Parent seed used to spawn every corruption step's generator.
    severity : str or None
        Name of the severity preset used, when applicable.
    sampled : mapping
        Concrete, JSON-serializable per-step summary statistics (e.g. the
        number of dropped stations, or the mean sampled distortion gain).
    """

    config: CorruptionConfig
    seed: int
    sampled: Mapping[str, Any] = field(default_factory=dict)
    severity: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable provenance record.

        Returns
        -------
        dict
            Schema-versioned config, seed, severity, and sampled summary.

        Examples
        --------
        >>> record = CorruptionRecord(CorruptionConfig(), seed=0)
        >>> record.to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "seed": self.seed,
            "severity": self.severity,
            "config_hash": self.config.config_hash(),
            "config": self.config.to_dict(),
            "sampled": dict(self.sampled),
        }


def add_heteroscedastic_noise(
    survey: SurveyData,
    *,
    level_range: tuple[float, float] = (0.02, 0.05),
    rng: np.random.Generator,
) -> SurveyData:
    """Add complex, per-observation heteroscedastic Gaussian noise.

    Parameters
    ----------
    survey : SurveyData
        Clean or already-corrupted survey.
    level_range : (float, float), default=(0.02, 0.05)
        Bounds on the relative noise standard deviation sampled
        independently for every ``(station, frequency)`` pair and shared
        across components at that pair, mimicking correlated instrument
        noise.
    rng : numpy.random.Generator
        Source of randomness; callers control reproducibility.

    Returns
    -------
    SurveyData
        New survey with perturbed impedance and an ``impedance_error`` that
        combines any pre-existing error with the injected noise in
        quadrature.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((1, 2, 1), 100 + 0j)
    >>> survey = SurveyData(z, [10.0, 1.0], ["S"], ["xy"], [[0.0, 0.0]])
    >>> noisy = add_heteroscedastic_noise(
    ...     survey, level_range=(0.05, 0.05), rng=np.random.default_rng(0)
    ... )
    >>> noisy.impedance_error is not None
    True
    """
    lo, hi = _range(level_range, "level_range")
    if lo == 0.0 and hi == 0.0:
        return survey
    n_station, n_frequency, n_component = survey.shape
    sigma_sf = rng.uniform(lo, hi, size=(n_station, n_frequency))
    sigma = np.broadcast_to(sigma_sf[:, :, None], survey.shape)
    magnitude = np.abs(survey.impedance)

    noise = (
        rng.normal(size=survey.shape) + 1j * rng.normal(size=survey.shape)
    ) / np.sqrt(2.0)
    injected_error = sigma * magnitude
    z_noisy = np.where(
        survey.valid,
        survey.impedance + injected_error * noise,
        survey.impedance,
    )

    if survey.impedance_error is None:
        error = np.where(survey.valid, injected_error, np.nan)
    else:
        error = np.where(
            survey.valid,
            np.sqrt(survey.impedance_error**2 + injected_error**2),
            survey.impedance_error,
        )
    return replace(survey, impedance=z_noisy, impedance_error=error)


def apply_error_floor(
    survey: SurveyData, *, floor_fraction: float
) -> SurveyData:
    """Clamp the declared error to a minimum fraction of ``|Z|``.

    Parameters
    ----------
    survey : SurveyData
        Survey whose error floor should be enforced.
    floor_fraction : float
        Minimum ``impedance_error`` as a fraction of ``|Z|``. Zero is a
        no-op.

    Returns
    -------
    SurveyData
        New survey with an error array at least as large as
        ``floor_fraction * |Z|`` on valid observations.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((1, 1, 1), 100 + 0j)
    >>> survey = SurveyData(
    ...     z,
    ...     [1.0],
    ...     ["S"],
    ...     ["xy"],
    ...     [[0.0, 0.0]],
    ...     impedance_error=np.ones((1, 1, 1)),
    ... )
    >>> floored = apply_error_floor(survey, floor_fraction=0.5)
    >>> floored.impedance_error[0, 0, 0]
    50.0
    """
    floor_fraction = _non_negative(floor_fraction, "floor_fraction")
    if floor_fraction == 0.0:
        return survey
    floor = floor_fraction * np.abs(survey.impedance)
    if survey.impedance_error is None:
        error = np.where(survey.valid, floor, np.nan)
    else:
        error = np.maximum(survey.impedance_error, floor)
    return replace(survey, impedance_error=error)


def apply_static_shift(
    survey: SurveyData,
    *,
    log10_sigma: float,
    rng: np.random.Generator,
    return_info: bool = False,
) -> SurveyData:
    """Multiply every present component by a per-station real factor.

    Static shift is modelled as a frequency-independent real scalar
    ``c_s = 10 ** N(0, log10_sigma)`` per station, applied identically to
    every impedance component and to the declared error, consistent with
    the linear scaling of a real multiplicative distortion.

    Parameters
    ----------
    survey : SurveyData
        Survey to distort.
    log10_sigma : float
        Std. dev. of the log\\ :sub:`10` static-shift factor. Zero is a
        no-op.
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    SurveyData
        New survey with static shift applied.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((2, 1, 1), 100 + 0j)
    >>> survey = SurveyData(z, [1.0], ["A", "B"], ["xy"], [[0, 0], [1, 0]])
    >>> shifted = apply_static_shift(
    ...     survey, log10_sigma=0.1, rng=np.random.default_rng(0)
    ... )
    >>> shifted.shape == survey.shape
    True
    """
    log10_sigma = _non_negative(log10_sigma, "log10_sigma")
    if log10_sigma == 0.0:
        factor = np.ones(survey.n_stations)
        out = survey
    else:
        factor = 10.0 ** rng.normal(0.0, log10_sigma, size=survey.n_stations)
        z = survey.impedance * factor[:, None, None]
        error = (
            None
            if survey.impedance_error is None
            else survey.impedance_error * factor[:, None, None]
        )
        out = replace(survey, impedance=z, impedance_error=error)
    if return_info:
        return out, {"static_shift_factor": factor}
    return out


def _distortion_matrix(
    gain: np.ndarray,
    twist_deg: np.ndarray,
    shear: np.ndarray,
    anisotropy: np.ndarray,
) -> np.ndarray:
    """Build per-station real 2x2 distortion matrices.

    ``D = gain * R(twist) @ [[1 + anisotropy, shear],
    [shear, 1 - anisotropy]]``, a simplified Groom & Bailey (1989)-style
    parameterisation used here only to inject distortion, not to
    reproduce :func:`pycsamt.emtools.gb.groom_bailey_table`'s
    decomposition exactly.
    """
    twist_rad = np.deg2rad(twist_deg)
    cos_t, sin_t = np.cos(twist_rad), np.sin(twist_rad)
    rotation = np.stack(
        [
            np.stack([cos_t, -sin_t], axis=-1),
            np.stack([sin_t, cos_t], axis=-1),
        ],
        axis=-2,
    )
    shear = np.clip(shear, -0.9, 0.9)
    anisotropy = np.clip(anisotropy, -0.9, 0.9)
    shape_matrix = np.stack(
        [
            np.stack([1.0 + anisotropy, shear], axis=-1),
            np.stack([shear, 1.0 - anisotropy], axis=-1),
        ],
        axis=-2,
    )
    return gain[:, None, None] * (rotation @ shape_matrix)


def apply_galvanic_distortion(
    survey: SurveyData,
    *,
    gain_log10_sigma: float = 0.0,
    twist_deg_sigma: float = 0.0,
    shear_sigma: float = 0.0,
    anisotropy_sigma: float = 0.0,
    rng: np.random.Generator,
    return_info: bool = False,
) -> SurveyData:
    """Inject a per-station real Groom & Bailey-style distortion matrix.

    Parameters
    ----------
    survey : SurveyData
        Survey to distort. Must expose at least one impedance component; any
        of ``xx, xy, yx, yy`` absent from :attr:`SurveyData.components` is
        treated as zero for the purpose of building the dense 2x2 impedance
        used internally, which is an approximation for surveys that only
        store off-diagonal components.
    gain_log10_sigma, twist_deg_sigma, shear_sigma, anisotropy_sigma : float
        Std. dev. of the per-station gain (log\\ :sub:`10`), twist (degrees),
        shear, and anisotropy parameters. Zero for every parameter is a
        no-op.
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    SurveyData
        New survey with distorted impedance; declared error, if any, is
        scaled by the sampled gain as a first-order approximation.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.ones((1, 1, 2), dtype=complex) * (1 + 1j)
    >>> survey = SurveyData(z, [1.0], ["S"], ["xy", "yx"], [[0, 0]])
    >>> distorted = apply_galvanic_distortion(
    ...     survey, twist_deg_sigma=10.0, rng=np.random.default_rng(0)
    ... )
    >>> distorted.shape == survey.shape
    True
    """
    n_station = survey.n_stations
    if (
        gain_log10_sigma == 0.0
        and twist_deg_sigma == 0.0
        and shear_sigma == 0.0
        and anisotropy_sigma == 0.0
    ):
        out = survey
        sampled = {
            "gain": np.ones(n_station),
            "twist_deg": np.zeros(n_station),
            "shear": np.zeros(n_station),
            "anisotropy": np.zeros(n_station),
        }
        return (out, sampled) if return_info else out

    gain = 10.0 ** rng.normal(0.0, gain_log10_sigma, size=n_station)
    twist_deg = rng.normal(0.0, twist_deg_sigma, size=n_station)
    shear = rng.normal(0.0, shear_sigma, size=n_station)
    anisotropy = rng.normal(0.0, anisotropy_sigma, size=n_station)
    D = _distortion_matrix(gain, twist_deg, shear, anisotropy)

    z_dense, present = _z_matrix(survey)
    distorted = np.einsum("sab,sfbc->sfac", D, z_dense)

    impedance = np.array(survey.impedance)
    for name, col in present.items():
        i, j = _COMPONENT_INDEX[name]
        impedance[:, :, col] = distorted[:, :, i, j]

    error = None
    if survey.impedance_error is not None:
        error = survey.impedance_error * gain[:, None, None]
    sampled = {
        "gain": gain,
        "twist_deg": twist_deg,
        "shear": shear,
        "anisotropy": anisotropy,
    }
    out = replace(survey, impedance=impedance, impedance_error=error)
    return (out, sampled) if return_info else out


def apply_dropout(
    survey: SurveyData,
    *,
    station_rate: float = 0.0,
    frequency_rate: float = 0.0,
    random_rate: float = 0.0,
    rng: np.random.Generator,
    return_info: bool = False,
) -> SurveyData:
    """Mark stations, frequencies, or individual observations as missing.

    Dropped observations are invalidated by setting the impedance to
    ``NaN``; :class:`~pycsamt.ai.data.contracts.SurveyData` construction
    then recomputes :attr:`~pycsamt.ai.data.contracts.SurveyData.valid`
    from finiteness, so masks stay authoritative automatically.

    Parameters
    ----------
    survey : SurveyData
        Survey to thin out.
    station_rate, frequency_rate, random_rate : float, default=0.0
        Independent probabilities that a whole station (all frequencies and
        components), a whole frequency (all stations and components), or an
        individual observation is dropped. Effects are combined (a station
        or frequency dropout wins over a random one at the same cell).
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    SurveyData
        New survey with additional invalid observations.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.ones((4, 4, 1), dtype=complex)
    >>> survey = SurveyData(
    ...     z,
    ...     [4.0, 3.0, 2.0, 1.0],
    ...     ["A", "B", "C", "D"],
    ...     ["xy"],
    ...     np.zeros((4, 2)),
    ... )
    >>> thinned = apply_dropout(
    ...     survey, station_rate=1.0, rng=np.random.default_rng(0)
    ... )
    >>> thinned.n_valid
    0
    """
    drop = np.zeros(survey.shape, dtype=bool)
    n_station, n_frequency, _ = survey.shape

    dropped_stations: list[int] = []
    if station_rate > 0.0:
        mask = rng.random(n_station) < station_rate
        drop |= mask[:, None, None]
        dropped_stations = np.flatnonzero(mask).tolist()

    dropped_frequencies: list[int] = []
    if frequency_rate > 0.0:
        mask = rng.random(n_frequency) < frequency_rate
        drop |= mask[None, :, None]
        dropped_frequencies = np.flatnonzero(mask).tolist()

    if random_rate > 0.0:
        drop |= rng.random(survey.shape) < random_rate

    impedance = np.where(drop, complex(np.nan, np.nan), survey.impedance)
    out = replace(survey, impedance=impedance)
    sampled = {
        "dropped_station_count": len(dropped_stations),
        "dropped_frequency_count": len(dropped_frequencies),
        "dropped_fraction": float(np.mean(drop)),
    }
    return (out, sampled) if return_info else out


def inject_outliers(
    survey: SurveyData,
    *,
    rate: float = 0.0,
    log10_shift_range: tuple[float, float] = (0.5, 1.5),
    rng: np.random.Generator,
    return_info: bool = False,
) -> SurveyData:
    """Perturb a random fraction of valid observations by a large factor.

    Outliers remain marked ``valid`` and keep their existing declared
    error, simulating a bad reading that quality control failed to flag —
    the case a robust inverter must tolerate.

    Parameters
    ----------
    survey : SurveyData
        Survey to perturb.
    rate : float, default=0.0
        Fraction of currently valid observations perturbed. Zero is a
        no-op.
    log10_shift_range : (float, float), default=(0.5, 1.5)
        Bounds on the outlier magnitude in log\\ :sub:`10` decades; the
        sign is randomized per outlier.
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    SurveyData
        New survey with a subset of valid impedance values shifted by
        ``10 ** (+/- shift)``.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((1, 10, 1), 100 + 0j)
    >>> survey = SurveyData(
    ...     z, np.arange(10.0, 0.0, -1.0), ["S"], ["xy"], [[0, 0]]
    ... )
    >>> corrupted = inject_outliers(
    ...     survey, rate=0.5, rng=np.random.default_rng(0)
    ... )
    >>> corrupted.shape == survey.shape
    True
    """
    rate = _rate(rate, "rate")
    if rate == 0.0:
        return (survey, {"n_outliers": 0}) if return_info else survey
    lo, hi = _range(log10_shift_range, "log10_shift_range")
    valid_idx = np.flatnonzero(survey.valid)
    n_outliers = int(round(rate * valid_idx.size))
    if n_outliers == 0:
        return (survey, {"n_outliers": 0}) if return_info else survey
    chosen = rng.choice(valid_idx, size=n_outliers, replace=False)
    shift = rng.uniform(lo, hi, size=n_outliers)
    sign = rng.choice([-1.0, 1.0], size=n_outliers)
    factor = 10.0 ** (sign * shift)

    impedance = np.array(survey.impedance)
    flat = impedance.reshape(-1)
    flat[chosen] = flat[chosen] * factor
    impedance = flat.reshape(survey.shape)
    out = replace(survey, impedance=impedance)
    return (out, {"n_outliers": n_outliers}) if return_info else out


def perturb_coordinates(
    survey: SurveyData,
    *,
    coordinate_sigma_m: float = 0.0,
    elevation_sigma_m: float = 0.0,
    rng: np.random.Generator,
) -> SurveyData:
    """Add Gaussian noise to station coordinates and elevation.

    Parameters
    ----------
    survey : SurveyData
        Survey whose station geometry should be perturbed.
    coordinate_sigma_m, elevation_sigma_m : float, default=0.0
        Std. dev. of Gaussian noise added to the horizontal (x, y)
        coordinates and to elevation, respectively. Elevation entries that
        are ``NaN`` (unknown) stay ``NaN``.
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    SurveyData
        New survey with perturbed :attr:`SurveyData.coordinates_m`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.ones((1, 1, 1), dtype=complex)
    >>> survey = SurveyData(z, [1.0], ["S"], ["xy"], [[0.0, 0.0, 100.0]])
    >>> moved = perturb_coordinates(
    ...     survey, coordinate_sigma_m=5.0, rng=np.random.default_rng(0)
    ... )
    >>> moved.coordinates_m.shape
    (1, 3)
    """
    coordinates = np.array(survey.coordinates_m)
    if coordinate_sigma_m > 0.0:
        coordinates[:, :2] += rng.normal(
            0.0, coordinate_sigma_m, size=(survey.n_stations, 2)
        )
    if elevation_sigma_m > 0.0:
        finite = np.isfinite(coordinates[:, 2])
        coordinates[finite, 2] += rng.normal(
            0.0, elevation_sigma_m, size=int(finite.sum())
        )
    return replace(survey, coordinates_m=coordinates)


SEVERITY_PRESETS: Mapping[str, CorruptionConfig] = {
    "clean": CorruptionConfig(),
    "in_distribution": CorruptionConfig(
        noise_level_range=(0.01, 0.03),
        error_floor_fraction=0.02,
        static_shift_log10_sigma=0.05,
        distortion_gain_log10_sigma=0.02,
        distortion_twist_deg_sigma=3.0,
        distortion_shear_sigma=0.05,
        random_dropout_rate=0.02,
        outlier_rate=0.0,
        coordinate_sigma_m=1.0,
    ),
    "severe": CorruptionConfig(
        noise_level_range=(0.05, 0.15),
        error_floor_fraction=0.05,
        static_shift_log10_sigma=0.15,
        distortion_gain_log10_sigma=0.05,
        distortion_twist_deg_sigma=10.0,
        distortion_shear_sigma=0.15,
        station_dropout_rate=0.05,
        frequency_dropout_rate=0.05,
        random_dropout_rate=0.05,
        outlier_rate=0.02,
        coordinate_sigma_m=5.0,
        elevation_sigma_m=2.0,
    ),
    "held_out_corruption": CorruptionConfig(
        noise_level_range=(0.10, 0.25),
        error_floor_fraction=0.05,
        static_shift_log10_sigma=0.30,
        distortion_gain_log10_sigma=0.10,
        distortion_twist_deg_sigma=25.0,
        distortion_shear_sigma=0.40,
        distortion_anisotropy_sigma=0.30,
        station_dropout_rate=0.15,
        frequency_dropout_rate=0.15,
        random_dropout_rate=0.10,
        outlier_rate=0.08,
        outlier_log10_shift_range=(1.0, 2.5),
        coordinate_sigma_m=15.0,
        elevation_sigma_m=5.0,
    ),
}
"""Named corruption suites for M3's clean/in-distribution/severe/OOD split."""


def apply_corruption_suite(
    survey: SurveyData,
    config: CorruptionConfig | None = None,
    *,
    severity: str | None = None,
    seed: int,
) -> tuple[SurveyData, CorruptionRecord]:
    """Apply the full, ordered M3 corruption pipeline from a single seed.

    Steps run in this fixed order: static shift, galvanic distortion,
    heteroscedastic noise, error floor, dropout, outliers, coordinate
    perturbation. Systematic distortions are applied to the clean signal
    before random noise and missingness, matching how these effects
    compose physically.

    Parameters
    ----------
    survey : SurveyData
        Clean survey to corrupt.
    config : CorruptionConfig, optional
        Explicit configuration. Mutually exclusive with ``severity``.
    severity : str, optional
        Name of an entry in :data:`SEVERITY_PRESETS` to use as ``config``.
    seed : int
        Parent seed. Each step draws from an independently spawned child
        generator so adding a new step does not change earlier steps'
        draws.

    Returns
    -------
    survey : SurveyData
        Corrupted survey.
    record : CorruptionRecord
        Provenance of the applied configuration and sampled parameters.

    Raises
    ------
    ValueError
        If both or neither of ``config``/``severity`` are given, or
        ``severity`` is unknown.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((3, 5, 2), 100 + 50j)
    >>> survey = SurveyData(
    ...     z,
    ...     np.linspace(100, 1, 5),
    ...     ["A", "B", "C"],
    ...     ["xy", "yx"],
    ...     np.zeros((3, 2)),
    ... )
    >>> corrupted, record = apply_corruption_suite(
    ...     survey, severity="in_distribution", seed=0
    ... )
    >>> record.severity
    'in_distribution'
    """
    if (config is None) == (severity is None):
        raise ValueError("exactly one of config or severity must be given.")
    if severity is not None:
        if severity not in SEVERITY_PRESETS:
            raise ValueError(
                f"unknown severity {severity!r}; choose from "
                f"{sorted(SEVERITY_PRESETS)}."
            )
        config = SEVERITY_PRESETS[severity]

    parent = np.random.default_rng(_seed(seed))
    sampled: dict[str, Any] = {}

    out = survey
    out, static_info = apply_static_shift(
        out,
        log10_sigma=config.static_shift_log10_sigma,
        rng=_child_rng(parent),
        return_info=True,
    )
    sampled["static_shift_factor_mean"] = float(
        np.mean(static_info["static_shift_factor"])
    )

    out, distortion_sampled = apply_galvanic_distortion(
        out,
        gain_log10_sigma=config.distortion_gain_log10_sigma,
        twist_deg_sigma=config.distortion_twist_deg_sigma,
        shear_sigma=config.distortion_shear_sigma,
        anisotropy_sigma=config.distortion_anisotropy_sigma,
        rng=_child_rng(parent),
        return_info=True,
    )
    sampled["distortion_twist_deg_mean"] = float(
        np.mean(distortion_sampled["twist_deg"])
    )

    if config.noise_level_range != (0.0, 0.0):
        out = add_heteroscedastic_noise(
            out, level_range=config.noise_level_range, rng=_child_rng(parent)
        )

    if config.error_floor_fraction > 0.0:
        out = apply_error_floor(
            out, floor_fraction=config.error_floor_fraction
        )

    out, dropout_sampled = apply_dropout(
        out,
        station_rate=config.station_dropout_rate,
        frequency_rate=config.frequency_dropout_rate,
        random_rate=config.random_dropout_rate,
        rng=_child_rng(parent),
        return_info=True,
    )
    sampled.update(dropout_sampled)

    out, outlier_sampled = inject_outliers(
        out,
        rate=config.outlier_rate,
        log10_shift_range=config.outlier_log10_shift_range,
        rng=_child_rng(parent),
        return_info=True,
    )
    sampled.update(outlier_sampled)

    if config.coordinate_sigma_m > 0.0 or config.elevation_sigma_m > 0.0:
        out = perturb_coordinates(
            out,
            coordinate_sigma_m=config.coordinate_sigma_m,
            elevation_sigma_m=config.elevation_sigma_m,
            rng=_child_rng(parent),
        )

    record = CorruptionRecord(
        config=config, seed=_seed(seed), sampled=sampled, severity=severity
    )
    return out, record

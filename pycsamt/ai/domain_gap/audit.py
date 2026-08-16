# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Survey-level audit reports for the M1 data-contract gate.

The M1 milestone requires that "an audit report accounts for every
included/excluded datum" before a survey is trusted as a training or
validation source. :func:`audit_survey` is that report: it runs a raw
``sites`` input (anything :func:`~pycsamt.emtools._core.ensure_sites`
accepts) through the same canonical bridge used everywhere else in
this package (:func:`~pycsamt.ai.domain_gap.survey_fit.survey_data_from_sites`)
*without* letting a single bad station abort the whole audit the way
that bridge's own fail-fast frequency-grid check does. Every station
dropped for missing data, and every station whose frequency grid does
not match the rest, is named and given a reason; dimensionality,
strike, static-shift, and distortion indicators come from the
existing, heavier pandas-based diagnostics in
:mod:`pycsamt.emtools.dimensionality`, :mod:`pycsamt.emtools.gb`, and
:mod:`pycsamt.emtools.ss` rather than being re-implemented here.

This module depends on those pandas-based diagnostics and therefore
lives in :mod:`pycsamt.ai.domain_gap`, not in the NumPy-only
:mod:`pycsamt.ai.data`, matching the layering
:mod:`pycsamt.ai.data.contracts` documents for itself.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from ...emtools._core import ensure_sites
from ...emtools.dimensionality import pre2d_inversion_assessment
from ..data.contracts import SurveyCoverage
from .survey_fit import (
    fit_distortion_priors_from_sites,
    survey_data_from_sites,
)

__all__ = [
    "StationExclusion",
    "FrequencyGridReport",
    "DimensionalitySummary",
    "SurveyAuditReport",
    "audit_survey",
]

_DEG_THRESHOLD_M = 1e-5  # ~1 m of latitude/longitude spread


def _station_coordinates_m(
    stations: list[Any], *, station_spacing: float
) -> tuple[np.ndarray, np.ndarray]:
    """Project station (lat, lon) to local metres; return (xy, elevation_m)."""
    n = len(stations)
    lats = np.full(n, np.nan)
    lons = np.full(n, np.nan)
    elevation = np.full(n, np.nan)
    for row, site in enumerate(stations):
        coords = getattr(site, "coords", None)
        if coords is None:
            continue
        try:
            if len(coords) >= 2:
                lats[row] = float(coords[0])
                lons[row] = float(coords[1])
            if len(coords) >= 3:
                elevation[row] = float(coords[2])
        except (TypeError, ValueError):
            continue
    ok = np.isfinite(lats) & np.isfinite(lons)
    lat_range = float(lats[ok].max() - lats[ok].min()) if ok.any() else 0.0
    lon_range = float(lons[ok].max() - lons[ok].min()) if ok.any() else 0.0
    if ok.sum() >= 2 and (
        lat_range > _DEG_THRESHOLD_M or lon_range > _DEG_THRESHOLD_M
    ):
        lat_ref = float(np.nanmean(lats[ok]))
        lon_ref = float(np.nanmean(lons[ok]))
        cos_lat = float(np.cos(np.radians(lat_ref)))
        x = np.where(ok, (lons - lon_ref) * cos_lat * 111320.0, np.nan)
        y = np.where(ok, (lats - lat_ref) * 111320.0, np.nan)
        xy = np.column_stack([x, y])
    else:
        side = int(np.ceil(np.sqrt(max(n, 1))))
        xs = np.array([(i % side) * station_spacing for i in range(n)])
        ys = np.array([(i // side) * station_spacing for i in range(n)])
        xy = np.column_stack([xs, ys])
    return xy, elevation


@dataclass(frozen=True)
class StationExclusion:
    """One station dropped before the canonical survey bridge.

    Parameters
    ----------
    station : str
        Station identifier as reported by the raw site collection.
    reason : str
        Human-readable, non-empty exclusion reason.

    Examples
    --------
    >>> StationExclusion("18-099Z", "missing freq or z array").reason
    'missing freq or z array'
    """

    station: str
    reason: str

    def __post_init__(self) -> None:
        station = str(self.station).strip()
        reason = str(self.reason).strip()
        if not station:
            raise ValueError("station cannot be empty.")
        if not reason:
            raise ValueError("reason cannot be empty.")
        object.__setattr__(self, "station", station)
        object.__setattr__(self, "reason", reason)

    def to_dict(self) -> dict[str, str]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            ``station`` and ``reason`` fields.

        Examples
        --------
        >>> StationExclusion("A", "missing freq or z array").to_dict()
        {'station': 'A', 'reason': 'missing freq or z array'}
        """
        return {"station": self.station, "reason": self.reason}


@dataclass(frozen=True)
class FrequencyGridReport:
    """Whether every included station shares one frequency grid.

    Parameters
    ----------
    matched : bool
        ``True`` only when every included station's frequency grid equals
        the reference station's, within tolerance.
    reference_station : str or None
        Station whose grid every other station was compared against.
    n_frequencies_by_station : mapping
        Station name to its own frequency count, for every included
        station, regardless of whether it matched.
    mismatched_stations : tuple of str
        Included stations whose frequency grid differs from the
        reference. Empty when ``matched`` is ``True``.

    Examples
    --------
    >>> report = FrequencyGridReport(True, "A", {"A": 10, "B": 10}, ())
    >>> report.matched
    True
    """

    matched: bool
    reference_station: str | None
    n_frequencies_by_station: Mapping[str, int]
    mismatched_stations: tuple[str, ...]

    def __post_init__(self) -> None:
        counts = {
            str(name): int(value)
            for name, value in self.n_frequencies_by_station.items()
        }
        mismatched = tuple(str(name) for name in self.mismatched_stations)
        if bool(self.matched) == bool(mismatched):
            raise ValueError(
                "matched must be true exactly when mismatched_stations is empty."
            )
        object.__setattr__(self, "matched", bool(self.matched))
        object.__setattr__(
            self,
            "reference_station",
            None
            if self.reference_station is None
            else str(self.reference_station),
        )
        object.__setattr__(self, "n_frequencies_by_station", counts)
        object.__setattr__(self, "mismatched_stations", mismatched)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            All fields as plain JSON-compatible values.

        Examples
        --------
        >>> FrequencyGridReport(True, "A", {"A": 1}, ()).to_dict()["matched"]
        True
        """
        return {
            "matched": self.matched,
            "reference_station": self.reference_station,
            "n_frequencies_by_station": dict(self.n_frequencies_by_station),
            "mismatched_stations": list(self.mismatched_stations),
        }


@dataclass(frozen=True)
class DimensionalitySummary:
    """Survey-wide aggregate of :func:`~pycsamt.emtools.dimensionality.pre2d_inversion_assessment`.

    Parameters
    ----------
    n_samples : int
        Total station-period samples the per-station table was built from.
    frac_1d, frac_2d, frac_3d : float
        Sample-weighted fraction classified 1-D, 2-D, and 3-D. They sum to
        one when ``n_samples`` is positive.
    strike_consensus_deg, strike_consensus_iqr_deg : float or None
        Median across stations of each station's consensus strike angle
        and its interquartile spread; ``None`` when no station produced a
        finite value.
    stations_recommending_3d_review : tuple of str
        Stations :func:`~pycsamt.emtools.dimensionality.pre2d_inversion_assessment`
        flagged ``"review_3d_effects_before_2d"``.

    Examples
    --------
    >>> summary = DimensionalitySummary(10, 0.7, 0.2, 0.1, 12.0, 5.0, ())
    >>> round(summary.frac_1d + summary.frac_2d + summary.frac_3d, 6)
    1.0
    """

    n_samples: int
    frac_1d: float
    frac_2d: float
    frac_3d: float
    strike_consensus_deg: float | None
    strike_consensus_iqr_deg: float | None
    stations_recommending_3d_review: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.n_samples < 0:
            raise ValueError("n_samples cannot be negative.")
        object.__setattr__(
            self,
            "stations_recommending_3d_review",
            tuple(str(name) for name in self.stations_recommending_3d_review),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            All fields as plain JSON-compatible values.

        Examples
        --------
        >>> DimensionalitySummary(1, 1.0, 0.0, 0.0, None, None).to_dict()[
        ...     "n_samples"
        ... ]
        1
        """
        return {
            "n_samples": self.n_samples,
            "frac_1d": self.frac_1d,
            "frac_2d": self.frac_2d,
            "frac_3d": self.frac_3d,
            "strike_consensus_deg": self.strike_consensus_deg,
            "strike_consensus_iqr_deg": self.strike_consensus_iqr_deg,
            "stations_recommending_3d_review": list(
                self.stations_recommending_3d_review
            ),
        }


@dataclass(frozen=True)
class SurveyAuditReport:
    """Complete M1 accounting of one survey's data quality.

    Parameters
    ----------
    n_stations_input : int
        Stations found by :func:`~pycsamt.emtools._core.ensure_sites`
        before any exclusion.
    excluded_stations : tuple of StationExclusion
        Every station dropped before the canonical bridge, with a reason.
    frequency_grid : FrequencyGridReport
        Whether every included station shares one frequency grid.
    coverage : SurveyCoverage or None
        Impedance coverage from :meth:`~pycsamt.ai.data.contracts.SurveyData.coverage`,
        or ``None`` when :attr:`frequency_grid` did not match (the
        canonical bridge cannot run without a shared grid).
    frequency_range_hz : (float, float) or None
        Minimum and maximum frequency, when available.
    error_ratio_p05, error_ratio_p50, error_ratio_p95 : float or None
        5th/50th/95th percentile of declared ``impedance_error / |Z|``
        over valid observations, when available.
    station_spacing_m : mapping or None
        ``min``, ``median``, and ``max`` consecutive station spacing in
        the survey's stored order.
    elevation_coverage : float
        Fraction of included stations with a finite elevation.
    crs_declared : bool
        Whether a coordinate reference system was supplied explicitly.
        Station coordinates for MT/AMT sites are ordinarily projected
        locally from latitude/longitude, which by itself declares no
        formal :term:`CRS`.
    dimensionality : DimensionalitySummary
        Aggregate dimensionality and strike indicators.
    static_shift_log10_sigma, distortion_gain_log10_sigma, distortion_twist_deg_sigma, distortion_shear_sigma, distortion_anisotropy_sigma : float
        Empirical spreads from :func:`~pycsamt.ai.domain_gap.survey_fit.fit_distortion_priors_from_sites`.
    generated_utc : str
        Timezone-aware ISO-8601 timestamp of when the audit ran.

    Examples
    --------
    Reports are normally produced by :func:`audit_survey`, not built
    directly:

    >>> report = SurveyAuditReport(
    ...     n_stations_input=1,
    ...     excluded_stations=(),
    ...     frequency_grid=FrequencyGridReport(True, "A", {"A": 1}, ()),
    ...     coverage=None,
    ...     frequency_range_hz=None,
    ...     error_ratio_p05=None,
    ...     error_ratio_p50=None,
    ...     error_ratio_p95=None,
    ...     station_spacing_m=None,
    ...     elevation_coverage=0.0,
    ...     crs_declared=False,
    ...     dimensionality=DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
    ...     static_shift_log10_sigma=0.0,
    ...     distortion_gain_log10_sigma=0.0,
    ...     distortion_twist_deg_sigma=0.0,
    ...     distortion_shear_sigma=0.0,
    ...     distortion_anisotropy_sigma=0.0,
    ...     generated_utc="2026-01-01T00:00:00Z",
    ... )
    >>> report.n_stations_included
    1
    """

    n_stations_input: int
    excluded_stations: tuple[StationExclusion, ...]
    frequency_grid: FrequencyGridReport
    coverage: SurveyCoverage | None
    frequency_range_hz: tuple[float, float] | None
    error_ratio_p05: float | None
    error_ratio_p50: float | None
    error_ratio_p95: float | None
    station_spacing_m: Mapping[str, float] | None
    elevation_coverage: float
    crs_declared: bool
    dimensionality: DimensionalitySummary
    static_shift_log10_sigma: float
    distortion_gain_log10_sigma: float
    distortion_twist_deg_sigma: float
    distortion_shear_sigma: float
    distortion_anisotropy_sigma: float
    generated_utc: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.n_stations_input < 0:
            raise ValueError("n_stations_input cannot be negative.")
        excluded = tuple(self.excluded_stations)
        if any(not isinstance(item, StationExclusion) for item in excluded):
            raise TypeError(
                "excluded_stations entries must be StationExclusion."
            )
        if not isinstance(self.frequency_grid, FrequencyGridReport):
            raise TypeError("frequency_grid must be a FrequencyGridReport.")
        if not isinstance(self.dimensionality, DimensionalitySummary):
            raise TypeError("dimensionality must be a DimensionalitySummary.")
        n_included = self.n_stations_input - len(excluded)
        if n_included < 0:
            raise ValueError(
                "excluded_stations cannot exceed n_stations_input."
            )
        coverage = self.coverage
        if coverage is not None and not isinstance(coverage, SurveyCoverage):
            raise TypeError("coverage must be a SurveyCoverage or None.")
        elevation_coverage = float(self.elevation_coverage)
        if not (0.0 <= elevation_coverage <= 1.0):
            raise ValueError("elevation_coverage must be in [0, 1].")
        object.__setattr__(self, "excluded_stations", excluded)
        object.__setattr__(self, "elevation_coverage", elevation_coverage)
        object.__setattr__(self, "generated_utc", str(self.generated_utc))
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def n_stations_included(self) -> int:
        """Return the number of stations that reached the canonical bridge.

        Returns
        -------
        int
            ``n_stations_input`` minus the number of exclusions.

        Examples
        --------
        Constructed directly for illustration; see :func:`audit_survey`
        for the normal way to obtain a report.

        >>> report = SurveyAuditReport(
        ...     2,
        ...     (StationExclusion("B", "missing freq or z array"),),
        ...     FrequencyGridReport(True, "A", {"A": 1}, ()),
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     0.0,
        ...     False,
        ...     DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     "2026-01-01T00:00:00Z",
        ... )
        >>> report.n_stations_included
        1
        """
        return self.n_stations_input - len(self.excluded_stations)

    def to_dict(self) -> dict[str, Any]:
        """Return a complete JSON-serializable representation.

        Returns
        -------
        dict
            Every field, with nested records converted to plain dicts.

        Examples
        --------
        >>> report = SurveyAuditReport(
        ...     1,
        ...     (),
        ...     FrequencyGridReport(True, "A", {"A": 1}, ()),
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     0.0,
        ...     False,
        ...     DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     "2026-01-01T00:00:00Z",
        ... )
        >>> report.to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "n_stations_input": self.n_stations_input,
            "n_stations_included": self.n_stations_included,
            "excluded_stations": [
                item.to_dict() for item in self.excluded_stations
            ],
            "frequency_grid": self.frequency_grid.to_dict(),
            "coverage": None
            if self.coverage is None
            else {
                "overall": self.coverage.overall,
                "by_station": self.coverage.by_station.tolist(),
                "by_frequency": self.coverage.by_frequency.tolist(),
                "by_component": self.coverage.by_component.tolist(),
            },
            "frequency_range_hz": None
            if self.frequency_range_hz is None
            else list(self.frequency_range_hz),
            "error_ratio_p05": self.error_ratio_p05,
            "error_ratio_p50": self.error_ratio_p50,
            "error_ratio_p95": self.error_ratio_p95,
            "station_spacing_m": None
            if self.station_spacing_m is None
            else dict(self.station_spacing_m),
            "elevation_coverage": self.elevation_coverage,
            "crs_declared": self.crs_declared,
            "dimensionality": self.dimensionality.to_dict(),
            "static_shift_log10_sigma": self.static_shift_log10_sigma,
            "distortion_gain_log10_sigma": self.distortion_gain_log10_sigma,
            "distortion_twist_deg_sigma": self.distortion_twist_deg_sigma,
            "distortion_shear_sigma": self.distortion_shear_sigma,
            "distortion_anisotropy_sigma": self.distortion_anisotropy_sigma,
            "generated_utc": self.generated_utc,
            "metadata": dict(self.metadata),
        }

    def write_json(self, path: str | Path, *, overwrite: bool = True) -> Path:
        """Write a deterministic, human-readable audit report file.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination JSON file.
        overwrite : bool, default=True
            Permit replacement of an existing file.

        Returns
        -------
        pathlib.Path
            Destination path.

        Raises
        ------
        FileExistsError
            If the destination exists and ``overwrite`` is false.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> report = SurveyAuditReport(
        ...     1,
        ...     (),
        ...     FrequencyGridReport(True, "A", {"A": 1}, ()),
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     0.0,
        ...     False,
        ...     DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     "2026-01-01T00:00:00Z",
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = report.write_json(Path(directory) / "audit.json")
        ...     loaded = SurveyAuditReport.read_json(path)
        >>> loaded.n_stations_input
        1
        """
        target = Path(path)
        if target.exists() and not overwrite:
            raise FileExistsError(f"audit report already exists: {target}")
        target.write_text(
            json.dumps(
                self.to_dict(), indent=2, sort_keys=True, ensure_ascii=False
            )
            + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> SurveyAuditReport:
        """Restore a report from its JSON representation.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`.

        Returns
        -------
        SurveyAuditReport
            Validated immutable report.

        Raises
        ------
        ValueError
            If the schema version is unsupported.

        Examples
        --------
        >>> report = SurveyAuditReport(
        ...     1,
        ...     (),
        ...     FrequencyGridReport(True, "A", {"A": 1}, ()),
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     0.0,
        ...     False,
        ...     DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     "2026-01-01T00:00:00Z",
        ... )
        >>> SurveyAuditReport.from_dict(report.to_dict()) == report
        True
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported SurveyAuditReport schema version.")
        coverage_data = data.get("coverage")
        coverage = (
            None
            if coverage_data is None
            else SurveyCoverage(
                overall=coverage_data["overall"],
                by_station=coverage_data["by_station"],
                by_frequency=coverage_data["by_frequency"],
                by_component=coverage_data["by_component"],
            )
        )
        frequency_range = data.get("frequency_range_hz")
        return cls(
            n_stations_input=data["n_stations_input"],
            excluded_stations=tuple(
                StationExclusion(item["station"], item["reason"])
                for item in data.get("excluded_stations", [])
            ),
            frequency_grid=FrequencyGridReport(**data["frequency_grid"]),
            coverage=coverage,
            frequency_range_hz=None
            if frequency_range is None
            else tuple(frequency_range),
            error_ratio_p05=data.get("error_ratio_p05"),
            error_ratio_p50=data.get("error_ratio_p50"),
            error_ratio_p95=data.get("error_ratio_p95"),
            station_spacing_m=data.get("station_spacing_m"),
            elevation_coverage=data["elevation_coverage"],
            crs_declared=data["crs_declared"],
            dimensionality=DimensionalitySummary(
                **{
                    **data["dimensionality"],
                    "stations_recommending_3d_review": tuple(
                        data["dimensionality"][
                            "stations_recommending_3d_review"
                        ]
                    ),
                }
            ),
            static_shift_log10_sigma=data["static_shift_log10_sigma"],
            distortion_gain_log10_sigma=data["distortion_gain_log10_sigma"],
            distortion_twist_deg_sigma=data["distortion_twist_deg_sigma"],
            distortion_shear_sigma=data["distortion_shear_sigma"],
            distortion_anisotropy_sigma=data["distortion_anisotropy_sigma"],
            generated_utc=data["generated_utc"],
            metadata=data.get("metadata", {}),
        )

    @classmethod
    def read_json(cls, path: str | Path) -> SurveyAuditReport:
        """Read and validate a UTF-8 JSON audit report.

        Parameters
        ----------
        path : str or pathlib.Path
            Existing report file written by :meth:`write_json`.

        Returns
        -------
        SurveyAuditReport
            Validated immutable report.

        Examples
        --------
        See :meth:`write_json` for a complete round trip.
        """
        return cls.from_dict(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

    def summary(self) -> str:
        """Return a compact, human-readable multi-line report.

        Returns
        -------
        str
            Plain-text accounting of inclusion, coverage, geometry,
            dimensionality, and distortion indicators.

        Examples
        --------
        >>> report = SurveyAuditReport(
        ...     1,
        ...     (),
        ...     FrequencyGridReport(True, "A", {"A": 1}, ()),
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     None,
        ...     0.0,
        ...     False,
        ...     DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None),
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     0.0,
        ...     "2026-01-01T00:00:00Z",
        ... )
        >>> print(report.summary())
        Survey audit (generated 2026-01-01T00:00:00Z)
          Stations: 1 input, 1 included, 0 excluded
          Frequency grid: matched
          CRS declared: False
          Elevation coverage: 0.0%
          Dimensionality: n=0, 1D=0.0%, 2D=0.0%, 3D=0.0%
          Static shift log10 sigma: 0.0000
          Distortion sigma: gain(log10)=0.0000, twist_deg=0.00, shear=0.0000, anisotropy=0.0000
        """
        lines = [f"Survey audit (generated {self.generated_utc})"]
        lines.append(
            f"  Stations: {self.n_stations_input} input, "
            f"{self.n_stations_included} included, "
            f"{len(self.excluded_stations)} excluded"
        )
        for item in self.excluded_stations:
            lines.append(f"    - {item.station}: {item.reason}")
        grid = self.frequency_grid
        lines.append(
            "  Frequency grid: "
            + (
                "matched"
                if grid.matched
                else f"MISMATCHED ({len(grid.mismatched_stations)} stations)"
            )
        )
        if self.coverage is not None:
            lines.append(
                f"  Impedance coverage: {self.coverage.overall * 100:.1f}%"
            )
        if self.frequency_range_hz is not None:
            lo, hi = self.frequency_range_hz
            lines.append(f"  Frequency range: {lo:g}-{hi:g} Hz")
        if self.error_ratio_p50 is not None:
            lines.append(
                "  Declared error / |Z|: p05="
                f"{self.error_ratio_p05:.4f}, p50={self.error_ratio_p50:.4f}, "
                f"p95={self.error_ratio_p95:.4f}"
            )
        if self.station_spacing_m is not None:
            lines.append(
                "  Station spacing (m): min="
                f"{self.station_spacing_m['min']:.1f}, "
                f"median={self.station_spacing_m['median']:.1f}, "
                f"max={self.station_spacing_m['max']:.1f}"
            )
        lines.append(f"  CRS declared: {self.crs_declared}")
        lines.append(
            f"  Elevation coverage: {self.elevation_coverage * 100:.1f}%"
        )
        dim = self.dimensionality
        lines.append(
            f"  Dimensionality: n={dim.n_samples}, "
            f"1D={dim.frac_1d * 100:.1f}%, 2D={dim.frac_2d * 100:.1f}%, "
            f"3D={dim.frac_3d * 100:.1f}%"
        )
        if dim.strike_consensus_deg is not None:
            lines.append(
                f"  Strike (consensus): {dim.strike_consensus_deg:.1f} deg "
                f"(IQR {dim.strike_consensus_iqr_deg:.1f} deg)"
            )
        if dim.stations_recommending_3d_review:
            lines.append(
                "  Stations flagged for 3-D review: "
                + ", ".join(dim.stations_recommending_3d_review)
            )
        lines.append(
            f"  Static shift log10 sigma: {self.static_shift_log10_sigma:.4f}"
        )
        lines.append(
            "  Distortion sigma: gain(log10)="
            f"{self.distortion_gain_log10_sigma:.4f}, "
            f"twist_deg={self.distortion_twist_deg_sigma:.2f}, "
            f"shear={self.distortion_shear_sigma:.4f}, "
            f"anisotropy={self.distortion_anisotropy_sigma:.4f}"
        )
        return "\n".join(lines)


def audit_survey(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    verbose: int = 0,
    freq_rtol: float = 1e-6,
    band: tuple[float, float] | None = None,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
    station_spacing_fallback: float = 500.0,
    metadata: Mapping[str, Any] | None = None,
) -> SurveyAuditReport:
    """Audit a raw survey and account for every included/excluded station.

    Unlike :func:`~pycsamt.ai.domain_gap.survey_fit.survey_data_from_sites`,
    which deliberately raises on a mismatched frequency grid so training
    can never proceed silently on inconsistent data, this function never
    raises for that reason: it is meant to be run *before* deciding
    whether a survey is fit for that stricter bridge, and reports a
    mismatch as a finding rather than an exception.

    Parameters
    ----------
    sites : Any
        Anything accepted by :func:`pycsamt.emtools._core.ensure_sites`.
    recursive, on_dup, verbose
        Forwarded to ``ensure_sites`` and the underlying diagnostics.
    freq_rtol : float, default=1e-6
        Relative tolerance used when comparing each station's frequency
        grid against the reference station's.
    band : (float, float), optional
        Period band in seconds forwarded to
        :func:`~pycsamt.emtools.dimensionality.pre2d_inversion_assessment`.
    skew_th, ellipt_th : float, optional
        Phase-tensor skew and ellipticity thresholds forwarded to the
        same dimensionality assessment.
    station_spacing_fallback : float, default=500.0
        Uniform-grid spacing in metres used only when no station reports
        usable latitude/longitude.
    metadata : mapping, optional
        Extra provenance recorded on the returned report.

    Returns
    -------
    SurveyAuditReport
        Complete accounting of inclusion, coverage, geometry,
        dimensionality, and distortion indicators.

    Raises
    ------
    ValueError
        If no station in ``sites`` has usable impedance data at all.

    Examples
    --------
    >>> report = audit_survey(
    ...     "data/AMT/WILLY_data/L18PLT", recursive=True, verbose=0
    ... )  # doctest: +SKIP
    >>> report.frequency_grid.matched  # doctest: +SKIP
    True
    """
    collection = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup, verbose=verbose
    )
    stations = list(collection)
    n_input = len(stations)

    excluded: list[StationExclusion] = []
    kept: list[Any] = []
    for site in stations:
        freq = getattr(site, "freq", None)
        z = getattr(site, "z", None)
        name = str(
            getattr(site, "name", f"station-{len(kept) + len(excluded)}")
        )
        if freq is None or z is None:
            excluded.append(StationExclusion(name, "missing freq or z array"))
            continue
        kept.append(site)
    if not kept:
        raise ValueError("no station in sites has usable impedance data.")

    n_freq_by_station: dict[str, int] = {}
    reference = np.asarray(kept[0].freq, dtype=float)
    reference_name = str(getattr(kept[0], "name", "station-0"))
    mismatched: list[str] = []
    for site in kept:
        name = str(getattr(site, "name", "?"))
        freq = np.asarray(site.freq, dtype=float)
        n_freq_by_station[name] = int(freq.size)
        if freq.shape != reference.shape or not np.allclose(
            freq, reference, rtol=freq_rtol, atol=0.0
        ):
            mismatched.append(name)
    grid_report = FrequencyGridReport(
        matched=not mismatched,
        reference_station=reference_name,
        n_frequencies_by_station=n_freq_by_station,
        mismatched_stations=tuple(mismatched),
    )

    coverage: SurveyCoverage | None = None
    frequency_range: tuple[float, float] | None = None
    error_p05 = error_p50 = error_p95 = None
    if grid_report.matched:
        survey = survey_data_from_sites(
            collection, recursive=False, on_dup=on_dup, verbose=verbose
        )
        coverage = survey.coverage()
        frequency_range = (
            float(np.min(survey.frequencies_hz)),
            float(np.max(survey.frequencies_hz)),
        )
        if survey.impedance_error is not None:
            valid = survey.valid
            ratio = np.abs(survey.impedance_error[valid]) / np.maximum(
                np.abs(survey.impedance[valid]), 1e-24
            )
            ratio = ratio[np.isfinite(ratio)]
            if ratio.size:
                error_p05, error_p50, error_p95 = (
                    float(np.percentile(ratio, 5)),
                    float(np.percentile(ratio, 50)),
                    float(np.percentile(ratio, 95)),
                )

    xy, elevation = _station_coordinates_m(
        kept, station_spacing=station_spacing_fallback
    )
    if len(kept) >= 2:
        steps = np.linalg.norm(np.diff(xy, axis=0), axis=1)
        steps = steps[np.isfinite(steps)]
    else:
        steps = np.empty(0)
    spacing = (
        None
        if steps.size == 0
        else {
            "min": float(np.min(steps)),
            "median": float(np.median(steps)),
            "max": float(np.max(steps)),
        }
    )
    elevation_coverage = (
        float(np.mean(np.isfinite(elevation))) if len(kept) else 0.0
    )

    dim_table = pre2d_inversion_assessment(
        collection,
        band=band,
        skew_th=skew_th,
        ellipt_th=ellipt_th,
        recursive=False,
        on_dup=on_dup,
        verbose=verbose,
        api=False,
    )
    if len(dim_table):
        n_samples = dim_table["n_samples"].to_numpy(dtype=float)
        total = float(n_samples.sum())
        if total > 0:
            frac_1d = float(
                (dim_table["frac_1d"].to_numpy() * n_samples).sum() / total
            )
            frac_2d = float(
                (dim_table["frac_2d"].to_numpy() * n_samples).sum() / total
            )
            frac_3d = float(
                (dim_table["frac_3d"].to_numpy() * n_samples).sum() / total
            )
        else:
            frac_1d = frac_2d = frac_3d = 0.0
        cons = dim_table["strike_consensus_deg"].to_numpy(dtype=float)
        cons_iqr = dim_table["strike_consensus_iqr_deg"].to_numpy(dtype=float)
        strike_deg = (
            float(np.nanmedian(cons)) if np.isfinite(cons).any() else None
        )
        strike_iqr = (
            float(np.nanmedian(cons_iqr))
            if np.isfinite(cons_iqr).any()
            else None
        )
        flagged = tuple(
            str(name)
            for name, recommendation in zip(
                dim_table["station"], dim_table["recommendation"]
            )
            if recommendation == "review_3d_effects_before_2d"
        )
        dimensionality = DimensionalitySummary(
            n_samples=int(total),
            frac_1d=frac_1d,
            frac_2d=frac_2d,
            frac_3d=frac_3d,
            strike_consensus_deg=strike_deg,
            strike_consensus_iqr_deg=strike_iqr,
            stations_recommending_3d_review=flagged,
        )
    else:
        dimensionality = DimensionalitySummary(0, 0.0, 0.0, 0.0, None, None)

    priors = fit_distortion_priors_from_sites(
        collection, recursive=False, on_dup=on_dup, verbose=verbose
    )

    generated_utc = (
        datetime.now(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )

    return SurveyAuditReport(
        n_stations_input=n_input,
        excluded_stations=tuple(excluded),
        frequency_grid=grid_report,
        coverage=coverage,
        frequency_range_hz=frequency_range,
        error_ratio_p05=error_p05,
        error_ratio_p50=error_p50,
        error_ratio_p95=error_p95,
        station_spacing_m=spacing,
        elevation_coverage=elevation_coverage,
        crs_declared=False,
        dimensionality=dimensionality,
        static_shift_log10_sigma=priors["static_shift_log10_sigma"],
        distortion_gain_log10_sigma=priors["distortion_gain_log10_sigma"],
        distortion_twist_deg_sigma=priors["distortion_twist_deg_sigma"],
        distortion_shear_sigma=priors["distortion_shear_sigma"],
        distortion_anisotropy_sigma=priors["distortion_anisotropy_sigma"],
        generated_utc=generated_utc,
        metadata=metadata or {},
    )

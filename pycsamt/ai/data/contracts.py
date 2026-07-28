# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Canonical, dependency-light data contracts for MT/AMT surveys.

The classes in this module define the boundary between survey ingestion and
all downstream AI or physics code.  They deliberately depend only on NumPy:
EDI readers, interpolation policies, neural-network frameworks, and Maxwell
solver backends belong elsewhere.

The canonical impedance axis order is ``(station, frequency, component)``.
Keeping this order explicit prevents a common class of scientifically quiet
errors in which stations, tensor components, or frequency order are swapped
without changing an array's rank.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

__all__ = [
    "ImpedanceConvention",
    "SurveyCoverage",
    "SurveyData",
    "merge_surveys",
]


def _readonly(array: Any, dtype: Any | None = None) -> np.ndarray:
    out = np.array(array, dtype=dtype, copy=True)
    out.setflags(write=False)
    return out


def _names(values: Sequence[str], expected: int, label: str) -> tuple[str, ...]:
    result = tuple(str(value).strip() for value in values)
    if len(result) != expected:
        raise ValueError(f"{label} must contain {expected} entries; got {len(result)}.")
    if any(not value for value in result):
        raise ValueError(f"{label} cannot contain empty names.")
    if len(set(result)) != len(result):
        raise ValueError(f"{label} must be unique.")
    return result


def _json_metadata(value: Mapping[str, Any]) -> Mapping[str, Any]:
    data = dict(value)
    try:
        encoded = json.dumps(data, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "metadata must contain finite JSON-serializable values."
        ) from exc
    return MappingProxyType(json.loads(encoded))


@dataclass(frozen=True)
class ImpedanceConvention:
    """Describe the sign, units, and rotation convention of impedance data.

    Parameters
    ----------
    time_dependence : {"exp(+iwt)", "exp(-iwt)"}, default="exp(+iwt)"
        Fourier time convention used for the stored complex impedance. Solver
        predictions must use the same convention before residuals are formed.
    units : {"V/A"}, default="V/A"
        Physical unit of electric field divided by magnetic field. The current
        canonical contract accepts SI impedance only.
    rotation_deg : float, default=0.0
        Clockwise rotation already applied to the horizontal tensor axes, in
        degrees. Values are normalized to the interval ``[0, 360)``.
    coordinate_orientation : str, default="x_north_y_east"
        Human-readable definition of the horizontal tensor axes.

    Examples
    --------
    Record a tensor rotated clockwise into a geological strike frame:

    >>> convention = ImpedanceConvention(rotation_deg=32.0)
    >>> convention.rotation_deg
    32.0
    >>> convention.to_dict()["time_dependence"]
    'exp(+iwt)'

    Notes
    -----
    This object records a convention; it does not rotate or conjugate data.
    Such transformations must be explicit preprocessing operations that create
    a new :class:`SurveyData` object with updated provenance.
    """

    time_dependence: str = "exp(+iwt)"
    units: str = "V/A"
    rotation_deg: float = 0.0
    coordinate_orientation: str = "x_north_y_east"

    def __post_init__(self) -> None:
        if self.time_dependence not in {"exp(+iwt)", "exp(-iwt)"}:
            raise ValueError("time_dependence must be 'exp(+iwt)' or 'exp(-iwt)'.")
        if self.units != "V/A":
            raise ValueError("units must be 'V/A' for canonical impedance.")
        rotation = float(self.rotation_deg)
        if not np.isfinite(rotation):
            raise ValueError("rotation_deg must be finite.")
        orientation = str(self.coordinate_orientation).strip()
        if not orientation:
            raise ValueError("coordinate_orientation cannot be empty.")
        object.__setattr__(self, "rotation_deg", rotation % 360.0)
        object.__setattr__(self, "coordinate_orientation", orientation)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            Convention fields with a schema discriminator.

        Examples
        --------
        >>> state = ImpedanceConvention().to_dict()
        >>> state["schema_version"], state["units"]
        (1, 'V/A')
        """
        return {
            "schema_version": 1,
            "time_dependence": self.time_dependence,
            "units": self.units,
            "rotation_deg": self.rotation_deg,
            "coordinate_orientation": self.coordinate_orientation,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ImpedanceConvention:
        """Restore and validate a serialized convention.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`.

        Returns
        -------
        ImpedanceConvention
            Validated immutable convention.

        Raises
        ------
        ValueError
            If the schema version is unsupported or a field is invalid.

        Examples
        --------
        >>> state = ImpedanceConvention(rotation_deg=15).to_dict()
        >>> ImpedanceConvention.from_dict(state).rotation_deg
        15.0
        """
        if data.get("schema_version") != 1:
            raise ValueError("unsupported ImpedanceConvention schema version.")
        return cls(
            time_dependence=data["time_dependence"],
            units=data["units"],
            rotation_deg=data["rotation_deg"],
            coordinate_orientation=data["coordinate_orientation"],
        )


@dataclass(frozen=True)
class SurveyCoverage:
    """Summarize the usable fraction of a survey along each data axis.

    Parameters
    ----------
    overall : float
        Fraction of valid impedance observations across the complete cube.
    by_station : ndarray, shape (n_station,)
        Valid fraction for each station.
    by_frequency : ndarray, shape (n_frequency,)
        Valid fraction for each frequency.
    by_component : ndarray, shape (n_component,)
        Valid fraction for each impedance component.
    tipper_overall : float or None, optional
        Fraction of valid tipper observations, or ``None`` when absent.

    Examples
    --------
    Coverage is normally obtained from :meth:`SurveyData.coverage`:

    >>> coverage = SurveyCoverage(1.0, [1.0], [1.0], [1.0])
    >>> coverage.complete
    True
    """

    overall: float
    by_station: np.ndarray
    by_frequency: np.ndarray
    by_component: np.ndarray
    tipper_overall: float | None = None

    def __post_init__(self) -> None:
        overall = float(self.overall)
        tipper = None if self.tipper_overall is None else float(self.tipper_overall)
        values = [overall] + ([] if tipper is None else [tipper])
        arrays = {}
        for name in ("by_station", "by_frequency", "by_component"):
            array = np.asarray(getattr(self, name), dtype=float)
            if array.ndim != 1 or array.size == 0:
                raise ValueError(f"{name} must be a non-empty 1-D array.")
            values.extend(array.tolist())
            arrays[name] = _readonly(array)
        if not np.all(np.isfinite(values)) or np.any(
            (np.asarray(values) < 0) | (np.asarray(values) > 1)
        ):
            raise ValueError("coverage fractions must be finite and in [0, 1].")
        object.__setattr__(self, "overall", overall)
        object.__setattr__(self, "tipper_overall", tipper)
        for name, array in arrays.items():
            object.__setattr__(self, name, array)

    @property
    def complete(self) -> bool:
        """Whether every impedance observation is usable.

        Returns
        -------
        bool
            ``True`` only when impedance coverage is exactly one. Tipper
            coverage is not included because tipper is optional.

        Examples
        --------
        >>> SurveyCoverage(0.5, [0.5], [0.5], [0.5]).complete
        False
        """
        return self.overall == 1.0


@dataclass(frozen=True)
class SurveyData:
    """Validated MT/AMT observations on a common survey grid.

    Parameters
    ----------
    impedance : array-like of complex, shape (n_station, n_frequency, n_component)
        Complex impedance in V/A. Invalid entries may be NaN but must be false
        in ``valid`` after construction.
    frequencies_hz : array-like, shape (n_frequency,)
        Positive, finite, unique, strictly monotonic frequencies.
    station_names, components : sequence of str
        Names corresponding exactly to the station and component axes.
    coordinates_m : array-like, shape (n_station, 2 or 3)
        Projected x/y coordinates and optional elevation in metres. A missing
        third column is represented by NaN elevation.
    impedance_error : array-like, optional
        Positive absolute standard errors with the same shape as impedance.
        Entries without usable errors are invalidated.
    valid : array-like of bool, optional
        Explicit observation mask. It is combined with finite-value and error
        checks; invalid data are never silently imputed.
    tipper, tipper_error, tipper_valid : array-like, optional
        Optional complex magnetic transfer functions shaped
        ``(n_station, n_frequency, 2)`` for Tx and Ty.
    crs : str, optional
        Coordinate reference system identifier. Projected coordinates should
        normally provide an EPSG or WKT identifier.
    metadata : mapping, optional
        Finite JSON-serializable provenance only.
    convention : ImpedanceConvention, optional
        Explicit complex sign, SI unit, tensor-axis, and rotation convention.

    Attributes
    ----------
    impedance : ndarray
        Read-only complex impedance cube.
    valid : ndarray of bool
        Read-only authoritative mask for usable impedance observations.

    Examples
    --------
    Construct a two-station survey with descending frequency order. Two-
    column coordinates are accepted and expanded with unknown elevations:

    >>> z = np.ones((2, 3, 2), dtype=complex) * (1 + 2j)
    >>> survey = SurveyData(
    ...     impedance=z,
    ...     frequencies_hz=[100.0, 10.0, 1.0],
    ...     station_names=["S01", "S02"],
    ...     components=["xy", "yx"],
    ...     coordinates_m=[[0.0, 0.0], [100.0, 0.0]],
    ...     crs="EPSG:32630",
    ... )
    >>> survey.shape
    (2, 3, 2)
    >>> survey.frequency_order
    'descending'

    Notes
    -----
    Construction copies all numerical inputs and marks them read-only. The
    dataclass is therefore safe to share between training, validation, and
    reporting code without accidental in-place mutation.
    """

    impedance: np.ndarray
    frequencies_hz: np.ndarray
    station_names: tuple[str, ...]
    components: tuple[str, ...]
    coordinates_m: np.ndarray
    impedance_error: np.ndarray | None = None
    valid: np.ndarray | None = None
    tipper: np.ndarray | None = None
    tipper_error: np.ndarray | None = None
    tipper_valid: np.ndarray | None = None
    crs: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)
    convention: ImpedanceConvention = field(default_factory=ImpedanceConvention)

    def __post_init__(self) -> None:
        z = np.asarray(self.impedance)
        if z.ndim != 3:
            raise ValueError(
                "impedance must have shape (station, frequency, component)."
            )
        if not np.issubdtype(z.dtype, np.complexfloating):
            z = z.astype(np.complex128)
        n_station, n_frequency, n_component = z.shape
        if min(z.shape) < 1:
            raise ValueError("impedance axes cannot be empty.")

        frequencies = np.asarray(self.frequencies_hz, dtype=float)
        if frequencies.shape != (n_frequency,):
            raise ValueError(f"frequencies_hz must have shape ({n_frequency},).")
        if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
            raise ValueError("frequencies_hz must be finite and positive.")
        differences = np.diff(frequencies)
        if differences.size and not (
            np.all(differences > 0) or np.all(differences < 0)
        ):
            raise ValueError("frequencies_hz must be strictly monotonic and unique.")

        station_names = _names(self.station_names, n_station, "station_names")
        components = _names(self.components, n_component, "components")
        coordinates = np.asarray(self.coordinates_m, dtype=float)
        if coordinates.shape == (n_station, 2):
            coordinates = np.column_stack([coordinates, np.full(n_station, np.nan)])
        if coordinates.shape != (n_station, 3):
            raise ValueError(f"coordinates_m must have shape ({n_station}, 2 or 3).")
        if not np.all(np.isfinite(coordinates[:, :2])):
            raise ValueError("coordinate x/y values must be finite.")

        finite_z = np.isfinite(z.real) & np.isfinite(z.imag)
        valid = (
            np.ones(z.shape, dtype=bool)
            if self.valid is None
            else np.asarray(self.valid, dtype=bool)
        )
        if valid.shape != z.shape:
            raise ValueError("valid must have the same shape as impedance.")
        valid = valid & finite_z

        error = None
        if self.impedance_error is not None:
            error = np.asarray(self.impedance_error, dtype=float)
            if error.shape != z.shape:
                raise ValueError(
                    "impedance_error must have the same shape as impedance."
                )
            valid = valid & np.isfinite(error) & (error > 0)

        tipper, tipper_error, tipper_valid = self._validate_tipper(
            n_station, n_frequency
        )
        crs = None if self.crs is None else str(self.crs).strip() or None
        convention = self.convention
        if isinstance(convention, Mapping):
            convention = ImpedanceConvention.from_dict(convention)
        if not isinstance(convention, ImpedanceConvention):
            raise TypeError("convention must be an ImpedanceConvention.")

        object.__setattr__(self, "impedance", _readonly(z))
        object.__setattr__(self, "frequencies_hz", _readonly(frequencies))
        object.__setattr__(self, "station_names", station_names)
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "coordinates_m", _readonly(coordinates))
        object.__setattr__(
            self,
            "impedance_error",
            None if error is None else _readonly(error),
        )
        object.__setattr__(self, "valid", _readonly(valid, bool))
        object.__setattr__(self, "tipper", tipper)
        object.__setattr__(self, "tipper_error", tipper_error)
        object.__setattr__(self, "tipper_valid", tipper_valid)
        object.__setattr__(self, "crs", crs)
        object.__setattr__(self, "metadata", _json_metadata(self.metadata))
        object.__setattr__(self, "convention", convention)

    def _validate_tipper(self, n_station: int, n_frequency: int):
        if self.tipper is None:
            if self.tipper_error is not None or self.tipper_valid is not None:
                raise ValueError("tipper_error/tipper_valid require tipper data.")
            return None, None, None
        t = np.asarray(self.tipper)
        expected = (n_station, n_frequency, 2)
        if t.shape != expected:
            raise ValueError(f"tipper must have shape {expected}.")
        if not np.issubdtype(t.dtype, np.complexfloating):
            t = t.astype(np.complex128)
        mask = np.isfinite(t.real) & np.isfinite(t.imag)
        if self.tipper_valid is not None:
            supplied = np.asarray(self.tipper_valid, dtype=bool)
            if supplied.shape != expected:
                raise ValueError(f"tipper_valid must have shape {expected}.")
            mask &= supplied
        error = None
        if self.tipper_error is not None:
            error = np.asarray(self.tipper_error, dtype=float)
            if error.shape != expected:
                raise ValueError(f"tipper_error must have shape {expected}.")
            mask &= np.isfinite(error) & (error > 0)
        return (
            _readonly(t),
            None if error is None else _readonly(error),
            _readonly(mask, bool),
        )

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return the canonical impedance shape.

        Returns
        -------
        tuple of int
            ``(n_station, n_frequency, n_component)``.

        Examples
        --------
        >>> z = np.ones((1, 2, 1), dtype=complex)
        >>> s = SurveyData(z, [10, 1], ["S"], ["xy"], [[0, 0]])
        >>> s.shape
        (1, 2, 1)
        """
        return self.impedance.shape

    @property
    def n_valid(self) -> int:
        """Return the number of usable impedance observations.

        Returns
        -------
        int
            Count of ``True`` entries in :attr:`valid`.

        Examples
        --------
        >>> z = np.array([[[1 + 1j], [complex(np.nan, np.nan)]]])
        >>> s = SurveyData(z, [10, 1], ["S"], ["xy"], [[0, 0]])
        >>> s.n_valid
        1
        """
        return int(np.count_nonzero(self.valid))

    @property
    def n_stations(self) -> int:
        """Return the number of stations.

        Returns
        -------
        int
            Length of the station axis.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((2, 1, 1), complex),
        ...     [1],
        ...     ["a", "b"],
        ...     ["xy"],
        ...     [[0, 0], [1, 0]],
        ... )
        >>> s.n_stations
        2
        """
        return self.shape[0]

    @property
    def n_frequencies(self) -> int:
        """Return the number of frequencies.

        Returns
        -------
        int
            Length of the frequency axis.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 2, 1), complex), [10, 1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> s.n_frequencies
        2
        """
        return self.shape[1]

    @property
    def n_components(self) -> int:
        """Return the number of impedance components.

        Returns
        -------
        int
            Length of the component axis.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 1, 2), complex), [1], ["S"], ["xy", "yx"], [[0, 0]]
        ... )
        >>> s.n_components
        2
        """
        return self.shape[2]

    @property
    def frequency_order(self) -> str:
        """Return the monotonic direction of the frequency axis.

        Returns
        -------
        {"ascending", "descending"}
            Direction in which frequencies are stored.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 2, 1), complex), [1, 10], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> s.frequency_order
        'ascending'
        """
        if self.n_frequencies == 1:
            return "ascending"
        return (
            "ascending"
            if self.frequencies_hz[1] > self.frequencies_hz[0]
            else "descending"
        )

    @property
    def has_tipper(self) -> bool:
        """Whether optional Tx/Ty transfer functions are present.

        Returns
        -------
        bool
            ``True`` when :attr:`tipper` is populated.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 1, 1), complex), [1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> s.has_tipper
        False
        """
        return self.tipper is not None

    def station_index(self, name: str) -> int:
        """Return the integer position of a named station.

        Parameters
        ----------
        name : str
            Exact, case-sensitive station identifier.

        Returns
        -------
        int
            Position along the station axis.

        Raises
        ------
        KeyError
            If ``name`` does not occur in :attr:`station_names`.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((2, 1, 1), complex),
        ...     [1],
        ...     ["A", "B"],
        ...     ["xy"],
        ...     [[0, 0], [1, 0]],
        ... )
        >>> s.station_index("B")
        1
        """
        try:
            return self.station_names.index(name)
        except ValueError as exc:
            raise KeyError(f"unknown station {name!r}.") from exc

    def component_index(self, name: str) -> int:
        """Return the integer position of an impedance component.

        Parameters
        ----------
        name : str
            Exact, case-sensitive component name such as ``"xy"``.

        Returns
        -------
        int
            Position along the component axis.

        Raises
        ------
        KeyError
            If ``name`` is not stored.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 1, 2), complex), [1], ["S"], ["xy", "yx"], [[0, 0]]
        ... )
        >>> s.component_index("yx")
        1
        """
        try:
            return self.components.index(name)
        except ValueError as exc:
            raise KeyError(f"unknown component {name!r}.") from exc

    def coverage(self) -> SurveyCoverage:
        """Calculate valid-data coverage along every impedance axis.

        Returns
        -------
        SurveyCoverage
            Overall, station, frequency, and component fractions. Optional
            tipper coverage is included when tipper data exist.

        Examples
        --------
        >>> z = np.ones((2, 2, 1), dtype=complex)
        >>> mask = np.array([[[True], [False]], [[True], [True]]])
        >>> s = SurveyData(
        ...     z, [10, 1], ["A", "B"], ["xy"], [[0, 0], [1, 0]], valid=mask
        ... )
        >>> s.coverage().overall
        0.75
        >>> s.coverage().by_station.tolist()
        [0.5, 1.0]
        """
        return SurveyCoverage(
            overall=float(np.mean(self.valid)),
            by_station=np.mean(self.valid, axis=(1, 2)),
            by_frequency=np.mean(self.valid, axis=(0, 2)),
            by_component=np.mean(self.valid, axis=(0, 1)),
            tipper_overall=None
            if self.tipper_valid is None
            else float(np.mean(self.tipper_valid)),
        )

    def component_data(
        self, name: str
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        """Return values, errors, and validity mask for one component.

        Parameters
        ----------
        name : str
            Exact component name.

        Returns
        -------
        values : ndarray, shape (n_station, n_frequency)
            Read-only complex impedance view.
        errors : ndarray or None
            Read-only absolute standard-error view, when available.
        valid : ndarray of bool
            Read-only validity-mask view.

        Raises
        ------
        KeyError
            If the named component is unavailable.

        Examples
        --------
        >>> z = np.ones((1, 2, 2), dtype=complex)
        >>> s = SurveyData(z, [10, 1], ["S"], ["xy", "yx"], [[0, 0]])
        >>> values, errors, valid = s.component_data("xy")
        >>> values.shape, errors, valid.all()
        ((1, 2), None, True)
        """
        index = self.component_index(name)
        errors = (
            None if self.impedance_error is None else self.impedance_error[:, :, index]
        )
        return self.impedance[:, :, index], errors, self.valid[:, :, index]

    def select(
        self,
        *,
        stations: Sequence[int] | None = None,
        frequencies: Sequence[int] | None = None,
        components: Sequence[int] | None = None,
    ) -> SurveyData:
        """Select survey axes by integer position.

        Parameters
        ----------
        stations, frequencies, components : sequence of int, optional
            Positions to retain. Omitted axes are retained completely. The
            requested order is preserved, but duplicates are rejected by the
            canonical unique-name/frequency validation.

        Returns
        -------
        SurveyData
            New immutable survey containing the requested subset.

        Raises
        ------
        ValueError
            If an index collection is not one-dimensional or creates a
            non-monotonic frequency axis.
        IndexError
            If an index lies outside its axis.

        Examples
        --------
        >>> z = np.ones((2, 3, 2), dtype=complex)
        >>> s = SurveyData(
        ...     z, [100, 10, 1], ["A", "B"], ["xy", "yx"], [[0, 0], [1, 0]]
        ... )
        >>> subset = s.select(stations=[1], frequencies=[1, 2], components=[0])
        >>> subset.station_names, subset.frequencies_hz.tolist()
        (('B',), [10.0, 1.0])
        """
        si = (
            np.arange(self.shape[0])
            if stations is None
            else np.asarray(stations, dtype=int)
        )
        fi = (
            np.arange(self.shape[1])
            if frequencies is None
            else np.asarray(frequencies, dtype=int)
        )
        ci = (
            np.arange(self.shape[2])
            if components is None
            else np.asarray(components, dtype=int)
        )
        if si.ndim != 1 or fi.ndim != 1 or ci.ndim != 1:
            raise ValueError("selection indices must be one-dimensional.")
        index = np.ix_(si, fi, ci)
        tip_index = np.ix_(si, fi, np.arange(2))
        return replace(
            self,
            impedance=self.impedance[index],
            frequencies_hz=self.frequencies_hz[fi],
            station_names=tuple(self.station_names[i] for i in si),
            components=tuple(self.components[i] for i in ci),
            coordinates_m=self.coordinates_m[si],
            impedance_error=None
            if self.impedance_error is None
            else self.impedance_error[index],
            valid=self.valid[index],
            tipper=None if self.tipper is None else self.tipper[tip_index],
            tipper_error=None
            if self.tipper_error is None
            else self.tipper_error[tip_index],
            tipper_valid=None
            if self.tipper_valid is None
            else self.tipper_valid[tip_index],
        )

    def select_names(
        self,
        *,
        stations: Sequence[str] | None = None,
        components: Sequence[str] | None = None,
        frequency_min_hz: float | None = None,
        frequency_max_hz: float | None = None,
    ) -> SurveyData:
        """Select stations/components by name and frequencies by interval.

        Parameters
        ----------
        stations, components : sequence of str, optional
            Exact names in the desired output order.
        frequency_min_hz, frequency_max_hz : float, optional
            Inclusive physical bounds. Their numerical order is independent of
            the stored frequency direction.

        Returns
        -------
        SurveyData
            New validated subset.

        Raises
        ------
        KeyError
            If a station or component name is unknown.
        ValueError
            If bounds are invalid or select no frequencies.

        Examples
        --------
        >>> z = np.ones((2, 3, 2), dtype=complex)
        >>> s = SurveyData(
        ...     z, [100, 10, 1], ["A", "B"], ["xy", "yx"], [[0, 0], [1, 0]]
        ... )
        >>> subset = s.select_names(
        ...     stations=["B"], components=["yx"], frequency_min_hz=5
        ... )
        >>> subset.shape
        (1, 2, 1)
        """
        si = (
            None
            if stations is None
            else [self.station_index(name) for name in stations]
        )
        ci = (
            None
            if components is None
            else [self.component_index(name) for name in components]
        )
        low = -np.inf if frequency_min_hz is None else float(frequency_min_hz)
        high = np.inf if frequency_max_hz is None else float(frequency_max_hz)
        if np.isnan(low) or np.isnan(high) or low <= 0 and np.isfinite(low):
            raise ValueError("frequency bounds must be positive when supplied.")
        if np.isfinite(high) and high <= 0:
            raise ValueError("frequency bounds must be positive when supplied.")
        if low > high:
            raise ValueError("frequency_min_hz cannot exceed frequency_max_hz.")
        fi = np.flatnonzero(
            (self.frequencies_hz >= low) & (self.frequencies_hz <= high)
        )
        if fi.size == 0:
            raise ValueError("frequency bounds select no observations.")
        return self.select(stations=si, frequencies=fi, components=ci)

    def assert_compatible(
        self,
        other: SurveyData,
        *,
        require_crs: bool = True,
        require_convention: bool = True,
    ) -> None:
        """Assert that two surveys share model-facing axes and conventions.

        Parameters
        ----------
        other : SurveyData
            Survey to compare.
        require_crs : bool, default=True
            Require identical CRS identifiers.
        require_convention : bool, default=True
            Require identical impedance conventions.

        Returns
        -------
        None
            Successful return means the surveys can share a fitted normalizer
            or be concatenated along their station axes.

        Raises
        ------
        TypeError
            If ``other`` is not :class:`SurveyData`.
        ValueError
            If frequency values/order, components/order, CRS, or convention
            differ under the requested policy.

        Examples
        --------
        >>> a = SurveyData(
        ...     np.ones((1, 1, 1), complex),
        ...     [1],
        ...     ["A"],
        ...     ["xy"],
        ...     [[0, 0]],
        ...     crs="EPSG:32630",
        ... )
        >>> b = SurveyData(
        ...     np.ones((1, 1, 1), complex),
        ...     [1],
        ...     ["B"],
        ...     ["xy"],
        ...     [[1, 0]],
        ...     crs="EPSG:32630",
        ... )
        >>> a.assert_compatible(b) is None
        True
        """
        if not isinstance(other, SurveyData):
            raise TypeError("other must be a SurveyData instance.")
        if not np.array_equal(self.frequencies_hz, other.frequencies_hz):
            raise ValueError("survey frequency grids or ordering differ.")
        if self.components != other.components:
            raise ValueError("survey component names or ordering differ.")
        if require_crs and self.crs != other.crs:
            raise ValueError("survey coordinate reference systems differ.")
        if require_convention and self.convention != other.convention:
            raise ValueError("survey impedance conventions differ.")

    def with_metadata(
        self, metadata: Mapping[str, Any], *, merge: bool = True
    ) -> SurveyData:
        """Return a copy with validated provenance metadata.

        Parameters
        ----------
        metadata : mapping
            Finite JSON-serializable values to add or use as replacement.
        merge : bool, default=True
            Merge with existing metadata when true. New keys replace existing
            keys. When false, replace the complete metadata mapping.

        Returns
        -------
        SurveyData
            New survey sharing no mutable numerical or metadata state.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 1, 1), complex), [1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> tagged = s.with_metadata({"line": "L18", "rotation_deg": 0.0})
        >>> tagged.metadata["line"]
        'L18'
        """
        updated = dict(self.metadata) if merge else {}
        updated.update(dict(metadata))
        return replace(self, metadata=updated)

    def summary(self) -> dict[str, Any]:
        """Return a compact JSON-serializable survey summary.

        Returns
        -------
        dict
            Axis sizes, bounds, coverage, CRS, component names, convention,
            and optional-data flags. Raw observations are not included.

        Examples
        --------
        >>> s = SurveyData(
        ...     np.ones((1, 2, 1), complex), [10, 1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> s.summary()["frequency_range_hz"]
        [1.0, 10.0]
        """
        coverage = self.coverage()
        return {
            "shape": list(self.shape),
            "n_valid": self.n_valid,
            "coverage": coverage.overall,
            "tipper_coverage": coverage.tipper_overall,
            "frequency_range_hz": [
                float(np.min(self.frequencies_hz)),
                float(np.max(self.frequencies_hz)),
            ],
            "frequency_order": self.frequency_order,
            "components": list(self.components),
            "crs": self.crs,
            "has_errors": self.impedance_error is not None,
            "has_tipper": self.has_tipper,
            "convention": self.convention.to_dict(),
        }

    def to_npz(self, path: str | Path) -> Path:
        """Write a lossless, pickle-free compressed archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination ``.npz`` path. NumPy appends ``.npz`` when the supplied
            filename has no such suffix.

        Returns
        -------
        pathlib.Path
            Requested destination path.

        Notes
        -----
        The archive contains only numerical arrays and JSON/Unicode scalars.
        :meth:`from_npz` therefore loads it with ``allow_pickle=False``.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> s = SurveyData(
        ...     np.ones((1, 1, 1), complex), [1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = s.to_npz(Path(directory) / "survey.npz")
        ...     restored = SurveyData.from_npz(path)
        >>> restored.station_names
        ('S',)
        """
        target = Path(path)
        payload = {
            "schema_version": np.array("2"),
            "impedance": self.impedance,
            "frequencies_hz": self.frequencies_hz,
            "station_names": np.asarray(self.station_names),
            "components": np.asarray(self.components),
            "coordinates_m": self.coordinates_m,
            "valid": self.valid,
            "crs": np.array(self.crs or ""),
            "metadata_json": np.array(json.dumps(dict(self.metadata), sort_keys=True)),
            "convention_json": np.array(
                json.dumps(self.convention.to_dict(), sort_keys=True)
            ),
        }
        for name in (
            "impedance_error",
            "tipper",
            "tipper_error",
            "tipper_valid",
        ):
            value = getattr(self, name)
            if value is not None:
                payload[name] = value
        np.savez_compressed(target, **payload)
        return target

    @classmethod
    def from_npz(cls, path: str | Path) -> SurveyData:
        """Load and validate a survey archive without enabling pickle.

        Parameters
        ----------
        path : str or pathlib.Path
            Archive previously written by :meth:`to_npz`.

        Returns
        -------
        SurveyData
            Newly validated immutable survey.

        Raises
        ------
        ValueError
            If the schema is unsupported or restored values violate the data
            contract.
        OSError
            If the path cannot be read as a NumPy archive.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> original = SurveyData(
        ...     np.ones((1, 1, 1), complex), [1], ["S"], ["xy"], [[0, 0]]
        ... )
        >>> with TemporaryDirectory() as directory:
        ...     path = original.to_npz(Path(directory) / "survey.npz")
        ...     loaded = SurveyData.from_npz(path)
        >>> np.array_equal(loaded.impedance, original.impedance)
        True

        Notes
        -----
        Schema version 1 remains readable and receives the default
        :class:`ImpedanceConvention`. New writes use schema version 2.
        """
        with np.load(Path(path), allow_pickle=False) as data:
            version = str(data["schema_version"].item())
            if version not in {"1", "2"}:
                raise ValueError(f"Unsupported SurveyData schema version {version!r}.")

            def optional(name):
                return data[name].copy() if name in data.files else None

            convention = (
                ImpedanceConvention()
                if version == "1"
                else ImpedanceConvention.from_dict(
                    json.loads(str(data["convention_json"].item()))
                )
            )
            return cls(
                impedance=data["impedance"],
                frequencies_hz=data["frequencies_hz"],
                station_names=tuple(data["station_names"].tolist()),
                components=tuple(data["components"].tolist()),
                coordinates_m=data["coordinates_m"],
                impedance_error=optional("impedance_error"),
                valid=data["valid"],
                tipper=optional("tipper"),
                tipper_error=optional("tipper_error"),
                tipper_valid=optional("tipper_valid"),
                crs=str(data["crs"].item()) or None,
                metadata=json.loads(str(data["metadata_json"].item())),
                convention=convention,
            )


def merge_surveys(
    surveys: Sequence[SurveyData],
    *,
    metadata: Mapping[str, Any] | None = None,
) -> SurveyData:
    """Concatenate compatible surveys along the station axis.

    Parameters
    ----------
    surveys : sequence of SurveyData
        Non-empty surveys with identical frequencies, component order, CRS,
        impedance convention, and optional-data availability. Station names
        must be globally unique.
    metadata : mapping, optional
        Metadata for the merged object. By default a minimal provenance record
        containing ``source_survey_count`` is used; input metadata are not
        combined implicitly because equal keys may have different meanings.

    Returns
    -------
    SurveyData
        Canonical survey with concatenated station-axis arrays.

    Raises
    ------
    ValueError
        If no surveys are supplied, surveys are incompatible, optional error
        or tipper availability differs, or station names overlap.
    TypeError
        If an item is not :class:`SurveyData`.

    Examples
    --------
    >>> a = SurveyData(
    ...     np.ones((1, 2, 1), complex), [10, 1], ["A"], ["xy"], [[0, 0]]
    ... )
    >>> b = SurveyData(
    ...     np.ones((1, 2, 1), complex), [10, 1], ["B"], ["xy"], [[1, 0]]
    ... )
    >>> merged = merge_surveys([a, b])
    >>> merged.station_names, merged.shape
    (('A', 'B'), (2, 2, 1))

    Notes
    -----
    This function performs no frequency interpolation, coordinate projection,
    tensor rotation, or unit conversion. Those operations must be explicit and
    completed before merging.
    """
    items = list(surveys)
    if not items:
        raise ValueError("surveys cannot be empty.")
    reference = items[0]
    if not isinstance(reference, SurveyData):
        raise TypeError("every item must be a SurveyData instance.")
    for survey in items[1:]:
        reference.assert_compatible(survey)
    attributes = ("impedance_error", "tipper", "tipper_error")
    for name in attributes:
        availability = [getattr(survey, name) is not None for survey in items]
        if any(availability) and not all(availability):
            raise ValueError(f"all surveys must have consistent {name} availability.")
    station_names = tuple(name for survey in items for name in survey.station_names)
    if len(set(station_names)) != len(station_names):
        raise ValueError("station names must be unique across merged surveys.")

    def _concatenate(name: str):
        first = getattr(reference, name)
        if first is None:
            return None
        return np.concatenate([getattr(survey, name) for survey in items], axis=0)

    return SurveyData(
        impedance=_concatenate("impedance"),
        frequencies_hz=reference.frequencies_hz,
        station_names=station_names,
        components=reference.components,
        coordinates_m=_concatenate("coordinates_m"),
        impedance_error=_concatenate("impedance_error"),
        valid=_concatenate("valid"),
        tipper=_concatenate("tipper"),
        tipper_error=_concatenate("tipper_error"),
        tipper_valid=_concatenate("tipper_valid"),
        crs=reference.crs,
        metadata={"source_survey_count": len(items)} if metadata is None else metadata,
        convention=reference.convention,
    )

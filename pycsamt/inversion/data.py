# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Data containers for physics-based EM inversion workflows."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

from .doc import _inversion_param_docs

__all__ = ["EMData"]


@dataclass
class EMData(PyCSAMTObject, MetadataMixin):
    method: str = "mt"
    frequencies: Any = None
    times: Any = None
    rho_a: Any = None
    phase: Any = None
    values: Any = None
    errors: Any = None
    station_names: list[str] = field(default_factory=list)
    station_x: Any = None
    source: Any = field(default=None, repr=False)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.method = str(self.method).lower()
        self.frequencies = _as_float_array(self.frequencies)
        self.times = _as_float_array(self.times)
        self.rho_a = _as_float_array(self.rho_a)
        self.phase = _as_float_array(self.phase)
        self.values = _as_float_array(self.values)
        self.errors = _as_float_array(self.errors)
        self.station_x = _as_float_array(self.station_x)
        self.station_names = [str(name) for name in self.station_names]
        self.validate()

    @classmethod
    def from_dict(cls, data: dict[str, Any], **overrides: Any) -> "EMData":
        """Build :class:`EMData` from a mapping.

        Parameters
        ----------
        data : mapping
            Mapping with observation fields. Keys may include ``method``,
            ``freqs``/``frequencies``, ``periods``, ``rho_a``, ``phase``,
            ``times``, ``values``, ``errors``, ``stations``/``station_names``,
            ``station_x``, and ``metadata``.
        **overrides
            Field values that override keys read from ``data``. This is useful
            when a caller wants to force ``method`` or attach metadata while
            preserving the original observation mapping.

        Returns
        -------
        EMData
            Normalized observation container.

        Notes
        -----
        The aliases ``freqs`` and ``stations`` are normalized to
        ``frequencies`` and ``station_names``. ``periods`` are converted to
        frequency in hertz.

        Examples
        --------
        >>> from pycsamt.inversion.data import EMData
        >>> data = EMData.from_dict({
        ...     "method": "mt",
        ...     "freqs": [1.0, 10.0],
        ...     "rho_a": [100.0, 120.0],
        ...     "phase": [45.0, 47.0],
        ...     "stations": ["S01"],
        ... })
        >>> data.frequencies.tolist()
        [1.0, 10.0]
        >>> data.station_names
        ['S01']
        """
        values = dict(data)
        if "freqs" in values and "frequencies" not in values:
            values["frequencies"] = values.pop("freqs")
        if "periods" in values and "frequencies" not in values:
            periods = np.asarray(values.pop("periods"), dtype=float)
            values["frequencies"] = 1.0 / periods
        if "stations" in values and "station_names" not in values:
            values["station_names"] = values.pop("stations")
        values.update(overrides)
        return cls(**values)

    @classmethod
    def coerce(cls, data: Any, *, method: str = "mt") -> "EMData":
        """Return *data* as an :class:`EMData` instance.

        Parameters
        ----------
        data : EMData, mapping, response object, sounding object, survey object, or sequence
            Input converted to :class:`EMData`. Supported objects include
            natural-source response objects with frequency and rho/phase
            attributes; station collections of such objects; TDEM sounding
            objects exposing time gates and decay values; and survey objects
            exposing ``to_soundings()``.
        method : str, default "mt"
            Fallback method used when *data* does not declare its own method.
            TDEM sounding/survey inputs are automatically promoted to
            ``"tdem"``.

        Returns
        -------
        EMData
            Existing ``EMData`` instances are returned unchanged. Recognized
            objects are parsed into normalized arrays. Unrecognized non-null
            inputs are retained as ``source`` on an otherwise empty container so
            backend adapters can still inspect them.

        Raises
        ------
        ValueError
            If station collections have incompatible frequency counts, TDEM
            soundings have incompatible time gates, or normalized arrays fail
            validation.

        Examples
        --------
        Coerce a response-like object::

            >>> import numpy as np
            >>> from pycsamt.inversion.data import EMData
            >>> class Response:
            ...     freqs = np.array([1.0, 10.0])
            ...     rho_a = np.array([100.0, 120.0])
            ...     phase = np.array([45.0, 47.0])
            ...     station_name = "S01"
            >>> data = EMData.coerce(Response(), method="mt")
            >>> data.station_names
            ['S01']

        Coerce a small TDEM survey object exposing ``to_soundings()``::

            >>> import numpy as np
            >>> from pycsamt.inversion.data import EMData
            >>> class Sounding:
            ...     station_name = "T01"
            ...     time_gates = np.array([1e-5, 1e-4])
            ...     data = np.array([1e-8, 2e-9])
            >>> class Survey:
            ...     def to_soundings(self):
            ...         return [Sounding()]
            >>> tdem = EMData.coerce(Survey())
            >>> tdem.method
            'tdem'
        """
        if isinstance(data, cls):
            return data
        if isinstance(data, dict):
            return cls.from_dict(data, method=data.get("method", method))
        parsed = _coerce_object(data, method=method)
        if parsed is not None:
            return parsed
        return cls(method=method, source=data)

    @property
    def n_samples(self) -> int:
        """Return the number of primary samples."""
        for arr in (self.frequencies, self.times, self.rho_a, self.values):
            if arr is not None:
                return int(np.asarray(arr).shape[-1])
        return 0

    @property
    def n_stations(self) -> int:
        """Return the number of stations/soundings represented."""
        for arr in (self.rho_a, self.phase, self.values):
            if arr is None:
                continue
            arr = np.asarray(arr)
            if arr.ndim <= 1:
                return 1
            return int(arr.shape[0])
        if self.station_x is not None:
            return int(np.asarray(self.station_x).size)
        if self.station_names:
            return len(self.station_names)
        return 1

    @property
    def has_mt_response(self) -> bool:
        """Whether apparent resistivity/phase observations are available."""
        return self.frequencies is not None and (
            self.rho_a is not None or self.phase is not None
        )

    @property
    def has_tdem_response(self) -> bool:
        """Whether time-domain EM decay observations are available."""
        return self.times is not None and self.values is not None

    def validate(self) -> None:
        """Validate basic shape consistency."""
        if self.frequencies is not None and np.any(self.frequencies <= 0):
            raise ValueError("frequencies must be strictly positive.")
        if self.times is not None and np.any(self.times <= 0):
            raise ValueError("times must be strictly positive.")
        n = None
        for name in ("rho_a", "phase", "values", "errors"):
            arr = getattr(self, name)
            if arr is None:
                continue
            if arr.ndim > 2:
                raise ValueError(f"{name} must be 1-D or 2-D.")
            if n is None:
                n = int(arr.shape[-1])
            elif int(arr.shape[-1]) != n:
                raise ValueError("rho_a, phase, and errors must share sample length.")
        if self.frequencies is not None and n is not None:
            if int(self.frequencies.size) != n:
                raise ValueError("frequencies must match data sample length.")
        if self.times is not None and self.values is not None:
            if int(self.times.size) != int(self.values.shape[-1]):
                raise ValueError("times must match values sample length.")
        n_st = self.n_stations
        if self.station_x is not None and self.station_x.size not in {0, n_st}:
            raise ValueError("station_x must match number of stations.")
        if self.station_names and len(self.station_names) != n_st:
            raise ValueError("station_names must match number of stations.")


EMData.__doc__ = """
Normalized EM observation container.

``EMData`` is the backend-neutral observation object used by
:mod:`pycsamt.inversion`. It keeps natural-source and time-domain EM data in one
small shape contract so built-in, optional, and external backends can share the
same input API.

Shape convention
----------------
Frequency-domain ``rho_a`` and ``phase`` arrays may be one-dimensional for a
single station or two-dimensional with shape ``(n_stations, n_frequencies)``.
TDEM ``values`` arrays may be one-dimensional for a single sounding or
two-dimensional with shape ``(n_soundings, n_times)``. Station metadata, when
provided, must match the first dimension of these arrays.

Parameters
----------
{method}
{frequencies}
{times}
{rho_a}
{phase}
{values}
{errors}
{station_names}
{station_x}
{source}
{metadata}

Notes
-----
Use :meth:`EMData.coerce` at API boundaries. It accepts dictionaries,
response-like objects, station collections, TDEM sounding collections, and TDEM
survey objects exposing ``to_soundings()``. The original input object is kept in
``source`` for provenance and backend-specific adapters.

{examples}

{references}
""".format(
    method=_inversion_param_docs.data.method,
    frequencies=_inversion_param_docs.data.frequencies,
    times=_inversion_param_docs.data.times,
    rho_a=_inversion_param_docs.data.rho_a,
    phase=_inversion_param_docs.data.phase,
    values=_inversion_param_docs.data.values,
    errors=_inversion_param_docs.data.errors,
    station_names=_inversion_param_docs.data.station_names,
    station_x=_inversion_param_docs.data.station_x,
    source=_inversion_param_docs.data.source,
    metadata=_inversion_param_docs.data.metadata,
    examples=_inversion_param_docs.data.examples,
    references=_inversion_param_docs.data.references,
)


def _as_float_array(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value, dtype=float)
    return arr.reshape(1) if arr.ndim == 0 else arr


def _coerce_object(data: Any, *, method: str) -> EMData | None:
    if data is None or isinstance(data, (str, bytes)):
        return None
    if _is_tdem_sounding(data):
        return _from_tdem_soundings([data], method=method, source=data)
    if hasattr(data, "to_soundings"):
        try:
            soundings = data.to_soundings()
        except TypeError:
            soundings = None
        if soundings:
            return _from_tdem_soundings(soundings, method="tdem", source=data)
    if _looks_like_response(data):
        return _from_response_object(data, method=method, source=data)
    if _is_iterable_collection(data):
        items = list(data)
        if not items:
            return None
        if all(_is_tdem_sounding(item) for item in items):
            return _from_tdem_soundings(items, method="tdem", source=data)
        if all(_looks_like_response(item) for item in items):
            return _from_response_collection(items, method=method, source=data)
    return None


def _from_response_object(obj: Any, *, method: str, source: Any) -> EMData:
    freqs = _frequencies(obj)
    rho, phase = _rho_phase(obj)
    errors = _attr(obj, ("errors", "error", "rho_error", "std", "standard_deviation"))
    station_names = _station_names_from_response(obj, _n_stations_from_arrays(rho, phase))
    station_x = _station_x_from_response(obj, len(station_names) or _n_stations_from_arrays(rho, phase))
    return EMData(
        method=_method_from_source(obj, method),
        frequencies=freqs,
        rho_a=_station_sample_matrix(rho, freqs),
        phase=_station_sample_matrix(phase, freqs),
        errors=_station_sample_matrix(errors, freqs),
        station_names=station_names,
        station_x=station_x,
        source=source,
        metadata=_metadata_from_source(obj, reader="response_object"),
    )


def _from_response_collection(items: list[Any], *, method: str, source: Any) -> EMData:
    freqs = _frequencies(items[0])
    rho_rows = []
    phase_rows = []
    error_rows = []
    names = []
    xs = []
    for idx, item in enumerate(items):
        item_freqs = _frequencies(item)
        if item_freqs is not None and freqs is not None:
            if np.asarray(item_freqs, dtype=float).size != np.asarray(freqs, dtype=float).size:
                raise ValueError("all station response objects must share frequency count.")
        rho, phase = _rho_phase(item)
        if rho is not None:
            rho_rows.append(np.asarray(rho, dtype=float).reshape(-1))
        if phase is not None:
            phase_rows.append(np.asarray(phase, dtype=float).reshape(-1))
        err = _attr(item, ("errors", "error", "rho_error", "std", "standard_deviation"))
        if err is not None:
            error_rows.append(np.asarray(err, dtype=float).reshape(-1))
        names.append(_station_name(item, idx))
        xs.append(_station_x_value(item, idx))
    return EMData(
        method=_method_from_source(source, method),
        frequencies=freqs,
        rho_a=np.vstack(rho_rows) if rho_rows else None,
        phase=np.vstack(phase_rows) if phase_rows else None,
        errors=np.vstack(error_rows) if error_rows else None,
        station_names=names,
        station_x=np.asarray(xs, dtype=float),
        source=source,
        metadata={"reader": "response_collection", "n_input_objects": len(items)},
    )


def _from_tdem_soundings(items: Iterable[Any], *, method: str, source: Any) -> EMData:
    soundings = list(items)
    if not soundings:
        return EMData(method="tdem", source=source)
    times = np.asarray(_attr(soundings[0], ("time_gates", "times", "time")), dtype=float)
    rows = []
    errors = []
    names = []
    xs = []
    for idx, sounding in enumerate(soundings):
        st_times = np.asarray(_attr(sounding, ("time_gates", "times", "time")), dtype=float)
        if st_times.size != times.size or not np.allclose(st_times, times):
            raise ValueError("all TEM soundings must share time gates.")
        values = _tdem_values(sounding)
        rows.append(np.asarray(values, dtype=float).reshape(-1))
        err = _attr(sounding, ("error", "errors", "std", "standard_deviation"))
        if err is not None:
            errors.append(np.asarray(err, dtype=float).reshape(-1))
        names.append(_station_name(sounding, idx))
        xs.append(_station_x_value(sounding, idx))
    return EMData(
        method="tdem" if method == "mt" else method,
        times=times,
        values=np.vstack(rows) if len(rows) > 1 else rows[0],
        errors=np.vstack(errors) if errors else None,
        station_names=names,
        station_x=np.asarray(xs, dtype=float),
        source=source,
        metadata={"reader": "tdem_soundings", "n_soundings": len(soundings)},
    )


def _is_iterable_collection(data: Any) -> bool:
    return isinstance(data, Iterable) and not isinstance(data, (str, bytes, dict, np.ndarray))


def _is_tdem_sounding(obj: Any) -> bool:
    return _attr(obj, ("time_gates", "times", "time")) is not None and (
        _attr(obj, ("data", "values", "dBz_dt", "dbdt")) is not None
        or callable(getattr(obj, "dBdt", None))
    )


def _looks_like_response(obj: Any) -> bool:
    return _frequencies(obj) is not None and any(
        _attr(obj, names) is not None
        for names in (
            ("rho_a", "rhoa", "app_res", "apparent_resistivity", "resistivity"),
            ("phase", "phi"),
            ("rho_a_te", "rho_te"),
            ("phase_te",),
            ("rho_a_tm", "rho_tm"),
            ("phase_tm",),
        )
    )


def _frequencies(obj: Any) -> np.ndarray | None:
    freqs = _attr(obj, ("frequencies", "freqs", "freq", "frequency"))
    if freqs is not None:
        return np.asarray(freqs, dtype=float)
    periods = _attr(obj, ("periods", "period"))
    if periods is not None:
        return 1.0 / np.asarray(periods, dtype=float)
    return None


def _rho_phase(obj: Any) -> tuple[Any, Any]:
    component = str(_metadata_from_source(obj).get("component", "te")).lower()
    rho = _attr(obj, ("rho_a", "rhoa", "app_res", "apparent_resistivity"))
    phase = _attr(obj, ("phase", "phi"))
    if rho is None:
        rho = _attr(obj, (f"rho_a_{component}", f"rho_{component}", "rho_a_te", "rho_te"))
    if phase is None:
        phase = _attr(obj, (f"phase_{component}", "phase_te"))
    return rho, phase


def _station_sample_matrix(values: Any, freqs: Any) -> np.ndarray | None:
    if values is None:
        return None
    arr = np.asarray(values, dtype=float)
    if arr.ndim <= 1:
        return arr
    n_freq = 0 if freqs is None else int(np.asarray(freqs).size)
    if n_freq and arr.shape[0] == n_freq and arr.shape[-1] != n_freq:
        return arr.T
    return arr


def _n_stations_from_arrays(*arrays: Any) -> int:
    for arr in arrays:
        if arr is None:
            continue
        arr = np.asarray(arr)
        if arr.ndim > 1:
            return int(min(arr.shape[0], arr.shape[1]))
    return 1


def _station_names_from_response(obj: Any, n_stations: int) -> list[str]:
    raw = _attr(obj, ("station_names", "stations", "site_names", "names"))
    if raw is not None:
        return [str(name) for name in np.asarray(raw).reshape(-1).tolist()]
    if n_stations <= 1:
        return [_station_name(obj, 0)]
    return [f"S{i:03d}" for i in range(n_stations)]


def _station_x_from_response(obj: Any, n_stations: int) -> np.ndarray | None:
    raw = _attr(obj, ("station_x", "stations_x", "x_stations", "x", "easting"))
    if raw is None:
        return None if n_stations <= 1 else np.arange(n_stations, dtype=float)
    arr = np.asarray(raw, dtype=float)
    if arr.ndim == 0:
        return arr.reshape(1)
    return arr


def _station_name(obj: Any, idx: int) -> str:
    raw = _attr(obj, ("station_name", "station", "name", "id", "site", "site_id"))
    if raw is None or isinstance(raw, (list, tuple, np.ndarray, dict)):
        return f"S{idx:03d}"
    return str(raw)


def _station_x_value(obj: Any, idx: int) -> float:
    raw = _attr(obj, ("station_x", "x", "easting", "longitude", "lon"))
    if raw is None:
        return float(idx)
    arr = np.asarray(raw, dtype=float)
    return float(arr.reshape(-1)[0])


def _tdem_values(obj: Any) -> np.ndarray:
    dbdt = getattr(obj, "dBdt", None)
    if callable(dbdt):
        try:
            return np.asarray(dbdt(), dtype=float)
        except Exception:
            pass
    return np.asarray(_attr(obj, ("dBz_dt", "dbdt", "data", "values")), dtype=float)


def _method_from_source(obj: Any, fallback: str) -> str:
    method = _attr(obj, ("method", "survey_type", "data_type"))
    if method is None:
        return fallback
    method = str(method).lower()
    if method in {"dbdt", "dhdt", "voltage", "normalized_voltage"}:
        return "tdem"
    return method


def _metadata_from_source(obj: Any, **extra: Any) -> dict[str, Any]:
    meta = _attr(obj, ("metadata", "meta"))
    if callable(meta):
        try:
            meta = meta()
        except Exception:
            meta = None
    out = dict(meta or {}) if isinstance(meta, dict) else {}
    out.update(extra)
    return out


def _attr(obj: Any, names: tuple[str, ...]) -> Any:
    for name in names:
        if isinstance(obj, dict) and name in obj:
            return obj[name]
        if hasattr(obj, name):
            value = getattr(obj, name)
            return value() if callable(value) and name in {"metadata", "meta"} else value
    return None

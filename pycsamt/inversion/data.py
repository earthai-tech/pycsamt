# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Data containers for physics-based EM inversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = ["EMData"]


@dataclass
class EMData(PyCSAMTObject, MetadataMixin):
    """Normalized EM observation container.

    Parameters
    ----------
    method : str
        Survey method tag, for example ``"mt"``, ``"amt"``, ``"emap"``,
        or ``"tdem"``.
    frequencies : array-like, optional
        Frequency samples in Hz for MT/AMT/EMAP data.
    times : array-like, optional
        Time samples in seconds for TDEM data.
    rho_a, phase : array-like, optional
        Apparent resistivity in ohm-m and impedance phase in degrees.
    values : array-like, optional
        Generic observed values, commonly used for TDEM decay data.
    errors : array-like or None
        Absolute data errors. If omitted, workflows may derive errors from
        their configured error floor.
    station_names : list of str
        Optional station labels.
    station_x : array-like
        Optional station profile coordinates in metres.
    source : object
        Original data object or path, retained for backend adapters.
    metadata : dict
        Free-form provenance and processing metadata.
    """

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
        """Build from a mapping, accepting common alias names."""
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
        """Return *data* as an :class:`EMData` instance."""
        if isinstance(data, cls):
            return data
        if isinstance(data, dict):
            return cls.from_dict(data, method=data.get("method", method))
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


def _as_float_array(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    arr = np.asarray(value, dtype=float)
    return arr.reshape(1) if arr.ndim == 0 else arr

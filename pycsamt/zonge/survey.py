# -*- coding: utf-8 -*-
#  Author : LKouadio <etanoyau@gmail.com>
#  License: LGPL-3.0-or-later
"""
Low-level *survey-layout* helpers used by Zonge/AVG readers.

* :class:`Station`   – receiver locations & spacing utilities
"""
from __future__ import annotations

import numpy as np
from typing import Sequence, Dict, Any, List

from .base import AVGComponentBase
from .utils import chunk_by_frequency
from ..exceptions import StationError

__all__ = ["Station"]


class Station(AVGComponentBase):
    """
    One-dimensional *survey-line* geometry container.

    A **station** is the projection of each receiver site onto the
    survey line (usually measured in metres along strike).

    Parameters
    ----------
    data : Sequence[float] | np.ndarray | None
        Raw *station* distances as recorded in the field.  If the same
        station is repeated *n_freq* times the array may contain
        duplicates (``[0, 0, 0, 25, 25, 25, …]``).
    unit : {'m', 'km', 'ft'}, default ``'m'``
        Distance unit used in *data*.  Internally everything is stored
        in **metres**.
    names : Sequence[str] | None, default *None*
        Optional custom station labels.  Must match the number of
        **unique** station positions.
    normalize : bool, default *False*
        Shift the origin so the first station is at ``0.0`` m.
    """

    _CONV_FT2M = 1. / 3.280839895
    _CONV_KM2M = 1e3


    def __init__(
        self,
        data:      Sequence[float] | np.ndarray | None = None,
        *,
        unit:      str                      = "m",
        names:     Sequence[str] | None     = None,
        normalize: bool                     = False,
        **kws:     Any,
        ) -> None:

        super().__init__(**kws)

        self._raw:   np.ndarray | None = None       # pristine distances (m)
        self.unit   = unit.lower()
        self.normalize  = normalize

        # geometry fields
        self.values: np.ndarray | None = None       # unique positions (m)
        self.names:  List[str] | None  = None
        self.loc:    Dict[str, np.ndarray] = {}

        self._names_in: Sequence[str] | None = names

        if data is not None:
            self.read(data, unit=unit, names=names, normalize=normalize)

    def read(
            self,
            data:      Sequence[float] | np.ndarray,
            *,
            unit:      str                  | None = None,
            names:     Sequence[str] | None = None,
            normalize: bool | None          = None,
        ) -> None:
        """Validate *data* and populate :pyattr:`loc`, :pyattr:`values`."""


        unit      = (unit or self.unit).lower()
        normalize = self.normalize if normalize is None else normalize
        arr_m     = self._to_metre(np.asarray(data, dtype=float).ravel(), unit)

        if arr_m.size == 0:
            raise StationError("empty station array")


        pos, counts       = np.unique(arr_m, return_counts=True)
        if np.any(counts != counts[0]):
            raise StationError("inconsistent frequency count per station")

        if normalize:
            pos   = pos - pos.min()
            arr_m = arr_m - arr_m.min()

        self._raw   = arr_m
        self.values = pos

        if names is not None:
            if len(names) != len(pos):
                raise StationError("`names` length does not match "
                                   "unique station count")
            use_names: List[str] = list(names)
        else:
            use_names = [f"S{idx:02d}" for idx in range(len(pos))]
        self.names = use_names

        # dictionary mapping 
        chunks         = chunk_by_frequency(arr_m, counts[0])
        self.loc       = dict(zip(use_names, chunks))

    def write(self) -> List[str]:
        """
        Serialise to a simple ``station=<value>`` list.  This placeholder
        is mainly meant for round-trips during refactoring.
        """
        if self._raw is None:
            return []
        return [f"station={v:.3f}" for v in self._raw]

    @property
    def min(self) -> float | None:
        """Smallest station distance (metres)."""
        return None if self.values is None else float(self.values.min())

    @property
    def max(self) -> float | None:
        """Largest station distance (metres)."""
        return None if self.values is None else float(self.values.max())

    @property
    def length(self) -> float | None:
        """Total line length (metres)."""
        return None if self.values is None else self.max - self.min

    def _to_metre(self, arr: np.ndarray, unit: str) -> np.ndarray:
        """Convert *arr* to metres depending on *unit* flag."""
        if unit == "m":
            return arr
        if unit == "km":
            return arr * self._CONV_KM2M
        if unit == "ft":
            return arr * self._CONV_FT2M
        raise StationError(f"unrecognised distance unit '{unit}'")


    def __len__(self) -> int:
        return 0 if self.values is None else self.values.size

    def __repr__(self) -> str:
        if self.values is None:
            return f"{self.__class__.__name__}(empty)"
        span = f"{self.min:g} – {self.max:g} m"
        return (f"{self.__class__.__name__}(n={len(self)}, "
                f"span={span}, names={bool(self._names_in)})")

# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0
"""
*Variation* metrics shared by every “quality-control” column
(percent-error, 100 × σ, …).

Only the abstract :class:`VariationBase` is defined here; concrete
sub-classes (e.g. *EmagVar*, *RhoVar*) will live in the same file and
inherit from it.
"""
from __future__ import annotations

from typing import Sequence, Dict, Any, List
import json 
import numpy as np

from .base  import AVGComponentBase
from .utils import chunk_by_frequency, number_stations


_RAD2DEG = 180.0 / np.pi
_MRAD    = 1e3

__all__= ['PcEmag','PcRho','PcHmag', 'PcHmag', 'SPhz','SHphz', 'SEphz'] 
 
class VariationBase(AVGComponentBase):
    """
    Abstract 1-D container for *per-sample* variation metrics.

    Parameters
    ----------
    data : Sequence[float] | np.ndarray | None
        Flattened values – one per *frequency × station* pair.
    n_freq : int | None, optional
        Number of discrete frequencies in the survey line.
    n_stations : int | None, optional
        Number of receiver stations.
    to_degree : bool, default ``False``
        Convert *data* from **milliradians** to **degrees** on load.

    Notes
    -----
    The class automatically computes extra descriptors
    (:pyattr:`mean`, :pyattr:`median`, :pyattr:`count`) once data are
    loaded.
    """

    def __init__(
            self,
            data:       Sequence[float] | np.ndarray | None = None,
            n_freq:     int  | None = None,
            n_stations: int  | None = None,
            to_degree:  bool        = False,
            **kwargs:   Any,
        ) -> None:

        super().__init__(**kwargs)
        self._raw           : np.ndarray | None = None
        self._to_degree     : bool              = to_degree
        self.n_freq         : int  | None       = n_freq
        self.n_stations     : int  | None       = n_stations

        # derived — populated by :meth:`read`
        self.values         : np.ndarray | None = None
        self.v_min          : float | None      = None
        self.v_max          : float | None      = None
        self.loc            : Dict[str, np.ndarray] = {}

        # quick stats
        self.mean           : float | None      = None
        self.median         : float | None      = None
        self.count          : int   | None      = None

        if data is not None:
            self.read(
                data,
                n_freq=n_freq,
                n_stations=n_stations,
            )

    def read(
        self,
        data: Sequence[float] | np.ndarray,
        *,
        n_freq:     int | None = None,
        n_stations: int | None = None,
        ) -> None:
        """
        Parse *data* and populate all derived attributes.

        Raises
        ------
        ValueError
            If shape constraints are violated.
        """
        vec = np.asarray(data, dtype=float).ravel()
        if self._to_degree:
            vec *= _RAD2DEG / 1.0e3                         # mrad → deg

        self.n_freq     = n_freq     or self.n_freq
        self.n_stations = n_stations or self.n_stations
        if self.n_freq is None or self.n_stations is None:
            raise ValueError("`n_freq` and `n_stations` must be given")

        expected = self.n_freq * self.n_stations
        if vec.size != expected:
            raise ValueError(
                f"Data length {vec.size} ≠ n_freq×n_stations "
                f"({self.n_freq}×{self.n_stations}={expected})"
            )

        self._raw   = vec
        self.values = np.unique(vec)
        self.v_min  = float(self.values.min())
        self.v_max  = float(self.values.max())

        blocks          = chunk_by_frequency(vec, self.n_freq)
        stn_ids, _      = number_stations(self.n_stations, self.n_freq)
        self.loc        = dict(zip(stn_ids, blocks))

        # quick stats
        self.mean   = float(vec.mean())
        self.median = float(np.median(vec))
        self.count  = int(vec.size)


    def write(self) -> List[str]:
        """
        Serialise to a compact JSON line.

        Returns
        -------
        list[str]
            A single-element list ready for file output.
        """
        if self._raw is None:
            return []

        payload = {
            "n_freq"    : self.n_freq,
            "n_stations": self.n_stations,
            "to_degree" : self._to_degree,
            "min"       : self.v_min,
            "max"       : self.v_max,
            "mean"      : self.mean,
            "median"    : self.median,
            "count"     : self.count,
            "data"      : self._raw.tolist(),
        }
        return [json.dumps(payload)]


    def __len__(self) -> int:                        # noqa: D401
        """Number of stored samples."""
        return 0 if self._raw is None else self._raw.size

    def __repr__(self) -> str:
        span = (self.v_min, self.v_max) if self.v_min is not None   \
               else ('?', '?')
        return (f"{self.__class__.__name__}(len={len(self)}, "
                f"range={span}, mean={self.mean}, "
                f"median={self.median})")


class _PcVarBase(VariationBase):
    """
    Base-class for the three *percent-error* columns
    (``%Emag``, ``%Hmag``, ``%Rho``).

    Only the *presentation* layer (label, *write*, *__str__*)
    is specialised here – all heavy lifting lives in
    :class:`~pycsamt.zonge.base.VariationBase`.
    """

    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
            *,
            label:     str                  = "E",
        ) -> None:

        self.label = label.upper()
        super().__init__(data, n_freq, n_station)


    def write(self) -> list[str]:                      # pragma: no cover
        return [f"pc{self.label.lower()}mag={v:g}" for v in self._raw]

    def __str__(self) -> str:
        rng = f"{self.v_min:g}–{self.v_max:g}%"
        return f"pc{self.label}mag({rng}, stations={len(self.loc)})"


class PcEmag(_PcVarBase):
    """Percent error on *E-field* magnitude."""
    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
        ) -> None:
        super().__init__(data, n_freq, n_station, label="E")


class PcHmag(_PcVarBase):
    """Percent error on *H-field* magnitude."""
    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
        ) -> None:
        super().__init__(data, n_freq, n_station, label="H")


class PcRho(_PcVarBase):
    """Percent error on apparent resistivity."""
    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
        ) -> None:
        super().__init__(data, n_freq, n_station, label="R")


class _StdVarBase(VariationBase):
    """
    Base-class for *phase* standard-deviation columns
    (``sPhz``, ``sEphz``, ``sHphz``).

    Parameters
    ----------
    data        : sequence, optional
        Raw values as read from the AVG file (radians).
    n_freq      : int, optional
        Number of frequencies per station.
    n_station   : int, optional
        Number of stations.
    label       : str
        ``''`` ⇒ *Z* tensor, ``'E'`` ⇒ E-field, ``'H'`` ⇒ H-field.
    to_degree   : bool, default *False*
        Convert output to **degrees** instead of milliradians.
    """

    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
            *,
            label:     str,
            to_degree: bool                 = False,
        ) -> None:

        self.label     = label
        self.to_degree = to_degree
        super().__init__(data, n_freq, n_station)


    def _postprocess(self, arr: np.ndarray) -> np.ndarray:   # noqa: D401
        """Convert *radians* → *mrad* (or deg) and return copy."""
        arr_mrad = arr * _MRAD
        if self.to_degree:
            return arr_mrad * _RAD2DEG / _MRAD
        return arr_mrad

    def write(self) -> List[str]:                # pragma: no cover
        key  = f"s{self.label.lower()}phz"
        unit = "deg" if self.to_degree else "mrad"
        return [f"{key}({unit})={v:g}" for v in self._raw]

    def __str__(self) -> str:
        unit = "deg" if self.to_degree else "mrad"
        rng  = f"{self.v_min:g}–{self.v_max:g} {unit}"
        return f"s{self.label}phz({rng}, stations={len(self.loc)})"


class SPhz(_StdVarBase):
    """100 · σ(phase) on the *Z* tensor."""
    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
            *,
            to_degree: bool                 = False,
        ) -> None:
        super().__init__(data, n_freq, n_station,
                         label="", to_degree=to_degree)


class SEphz(_StdVarBase):
    """100 · σ(phase) on the *E-field* column."""
    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
            *,
            to_degree: bool                 = False,
        ) -> None:
        super().__init__(data, n_freq, n_station,
                         label="E", to_degree=to_degree)


class SHphz(_StdVarBase):
    """100 · σ(phase) on the *H/B-field* column."""
    def __init__(
            self,
            data:      Sequence[Any] | None = None,
            n_freq:    int | None           = None,
            n_station: int | None           = None,
            *,
            to_degree: bool                 = False,
        ) -> None:
        super().__init__(data, n_freq, n_station,
                         label="H", to_degree=to_degree)

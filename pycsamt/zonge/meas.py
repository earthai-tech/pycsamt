# -*- coding: utf-8 -*-
#  Author : LKouadio <etanoyau@gmail.com>
#  License: LGPL-3.0

from __future__ import annotations

from typing import Sequence, List, Any
import numpy as np
               
from ..exceptions import InputError, FrequencyError
from .base import AVGComponentBase
from .var import VariationBase  


__all__ = ["CompMeas", "Amps", "Frequency"]



class CompMeas(AVGComponentBase):
    """
    Simple **enumeration wrapper** for the four classical CSAMT
    component pairs.

    The object merely stores the *label* and validates user input so
    downstream code can rely on a consistent name.

    Valid labels
    ------------
    ``{'ExHy', 'EyHx', 'EyHy', 'ExHx'}``
    """

    _VALID: set[str] = {"ExHy", "EyHx", "EyHy", "ExHx"}

    def __init__(self,
                 name: str | None = None,
                 **kws: Any) -> None:
        super().__init__(**kws)
        self.name = None
        if name is not None:
            self.read(name)

    def read(self, label: str, **kws: Any) -> None:           # noqa: D401
        """Validate *label* and store it."""
        lab = str(label).strip()
        if lab not in self._VALID:
            raise InputError(f"unrecognised component '{lab}' – "
                             f"expect one of {sorted(self._VALID)}")
        self.name = lab

    def write(self) -> List[str]:
        """Serialise as a single ``comp=<label>`` line."""
        return [] if self.name is None else [f"comp={self.name}"]


    def __str__(self) -> str:
        return self.name or "<unset>"


class Amps(VariationBase):
    """
    Transmitter **current amplitude** per frequency / station.

    The physical unit is **ampere** (A).  The layout follows *kind-1/2*
    ordering – i.e. a flat vector with *n_freq × n_stations* entries.

    Examples
    --------
    >>> amps = Amps([1.2, 1.2, 0.9, 0.9], n_freq=2, n_stations=2)
    >>> amps.min, amps.max
    (0.9, 1.2)
    >>> amps.loc['S01']
    array([0.9, 0.9])
    """

    def __init__(self,
                 data: Sequence[float] | np.ndarray | None = None,
                 *,
                 n_freq:     int | None = None,
                 n_stations: int | None = None,
                 **kws:      Any) -> None:

        super().__init__(
            data,
            n_freq=n_freq,
            n_stations=n_stations,
            to_degree=False,
            **kws
    )

    # VariationBase implements `read`; we only adapt the writer.
    def write(self) -> List[str]:
        """Return ``amps=<value>`` records suitable for round-tripping."""
        if self._raw is None:
            return []
        return [f"amps={v:g}" for v in self._raw]

    # Optional: tighter string representation

    def __str__(self) -> str:
        if self._raw is None:
            return "Amps(<empty>)"
        span = f"{self.v_min:g}–{self.v_max:g} A"
        return f"Amps(n={len(self)}, span={span})"


class Frequency(VariationBase):
    """
    Absolute transmit **frequency** at each station (*Hz*).

    The data layout is identical to the other 1-D containers: a flat
    vector with *n_freq × n_stations* elements – first all frequencies
    for *S00*, then *S01*, and so on.

    Notes
    -----
    * Values **must** be strictly positive.
    * Convenience :py:meth:`logspace` returns a canonical,
      decade-spaced array useful for interpolation / plotting.
    """

    def __init__(
        self,
        data: Sequence[float] | np.ndarray | None = None,
        *,
        n_freq:     int | None = None,
        n_stations: int | None = None,
        **kws:      Any
    ) -> None:

        super().__init__(
            data,
            n_freq=n_freq,
            n_stations=n_stations,
            to_degree=False,
            **kws)

    def read(self,
             data:       Sequence[float] | np.ndarray,
             *,
             n_freq:     int | None = None,
             n_stations: int | None = None) -> None:
        """See :py:meth:`VariationBase.read` (positivity enforced)."""
        super().read(data, n_freq=n_freq, n_stations=n_stations)
        if np.any(self._raw <= 0):
            raise FrequencyError("frequency values must be > 0 Hz")

    def logspace(self,
                 decade_start: int,
                 decade_stop:  int,
                 n_points:     int) -> np.ndarray:
        """
        Return **log-spaced** frequency grid.

        Parameters
        ----------
        decade_start, decade_stop : int
            10-exponent of first / last frequency (e.g. ``0 -> 1 Hz``).
        n_points : int
            Number of points in the output vector.

        Examples
        --------
        >>> Frequency().logspace(0, 4, 17)[:4]
        array([ 1.     ,  1.77828,  3.16228,  5.62341])
        """
        if n_points < 2:
            raise ValueError("`n_points` must be >= 2")
        grid = np.logspace(decade_start, decade_stop, n_points,
                           endpoint=True, base=10.0)
        return grid

    def write(self) -> List[str]:
        """Serialise as ``freq=<Hz>`` lines (kind-2 flavour)."""
        return [] if self._raw is None else [f"freq={v:g}" for v in self._raw]

    def __str__(self) -> str:
        if self._raw is None:
            return "Frequency(<empty>)"
        span = f"{self.v_min:g}–{self.v_max:g} Hz"
        return f"Frequency(n={len(self)}, span={span})"


# Square-wave transmitter current (A)
_MISSING = object()        # sentinel
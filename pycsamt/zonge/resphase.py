# -*- coding: utf-8 -*-
#  Author : LKouadio <etanoyau@gmail.com>
#  License: LGPL-3.0

from __future__ import annotations

from typing import Sequence, Any, Dict, List
import numpy as np

from ..exceptions  import ResistivityError, PhaseError
from .var  import VariationBase         
from .utils   import number_stations, chunk_by_frequency
from .tensor import TensorFactory, ImpedanceTensor

__all__= ["Resistivity", "Phase"] 


class Phase(VariationBase):
    """
    Absolute *Z-tensor* phase (ϕ) per station / frequency.

    Parameters
    ----------
    data : Sequence[float] | np.ndarray | None
        Flat ``mrad`` vector – one element per *(freq, station)* pair.
    n_freq, n_stations : int | None, optional
        Grid dimensions.  Required on first :py:meth:`read`.
    to_degree : bool, default *False*
        Return phase in **degree** instead of the native **mrad**.

    Notes
    -----
    The input is **milli-radians** as produced by Zonge’s AVG files.
    Internally we convert to **radians** before delegating to
    :class:`VariationBase`; the optional degree conversion is then
    handled there.
    """

    def __init__(self,
                 data: Sequence[float] | np.ndarray | None = None,
                 *,
                 n_freq:     int | None = None,
                 n_stations: int | None = None,
                 to_degree:  bool       = False,
                 **kws:      Any) -> None:

        super().__init__(None,
                         n_freq=n_freq,
                         n_stations=n_stations,
                         to_degree=to_degree,
                         **kws)
        if data is not None:
            self.read(data,
                      n_freq=n_freq,
                      n_stations=n_stations)


    def read(self,
             data:  Sequence[float] | np.ndarray,
             *,
             n_freq:     int | None = None,
             n_stations: int | None = None) -> None:
        """
        Load phase data (*mrad*) and populate internal caches.

        All validation is eventually performed by the parent
        :py:meth:`VariationBase.read`.
        """
        # mrad → rad
        vec_rad = np.asarray(data, dtype=float).ravel() * 1.0e-3
        super().read(vec_rad,n_freq=n_freq, n_stations=n_stations)

    def write(self) -> List[str]:
        """Serialise as ``phase=<val>`` rows (native units)."""
        if self._raw is None:                        # pragma: no cover
            return []
        # restore mrad regardless of display unit
        out = self._raw * 1.0e3
        return [f"phase={v:g}" for v in out]


    def as_tensor(self) -> "ImpedanceTensor":
        """
        Return the phase grid as an :class:`~ImpedanceTensor`.

        Only the *phase* block is populated – the constructor happily
        accepts *None* for the other stacks (*z*, ρ, errors…).
        """
        if self._raw is None:
            raise PhaseError("Nothing loaded - call `read()` first")

        if self.n_freq is None or self.n_stations is None:
            raise PhaseError("`n_freq` / `n_stations` unknown")

        # reshape flat (f ⋅ s) → (f, s)
        phase_fs = self._raw.reshape(self.n_freq, self.n_stations)

        phase_dict = {
            "xx": phase_fs,           # NB: same slab for all four slots
            "xy": phase_fs,
            "yx": phase_fs,
            "yy": phase_fs,
        }

        return TensorFactory.build(phase=phase_dict)

    @classmethod
    def from_tensor(cls,
                    tensor: "ImpedanceTensor",
                    *,
                    to_degree: bool = False
                    ) -> "Phase":
        """
        Build a :class:`Phase` instance out of an existing
        :class:`~ImpedanceTensor`.

        Parameters
        ----------
        tensor : ImpedanceTensor
            Source object – must contain a *phase* block.
        to_degree : bool, default *False*
            Whether the resulting :class:`Phase` object should operate
            in degree rather than milli-radian.
        """
        if tensor.phase is None:
            raise PhaseError("Input tensor carries no phase information")

        n_freq, n_stn = tensor.phase.shape[:2]

        # take the *xx* component – all four are identical for our
        # simple AVG format; ravel back to flat (f ⋅ s)
        data_mrad = tensor.phase[:, :, 0, 0].ravel() * 1.0e3

        return cls(data_mrad,
                   n_freq=n_freq,
                   n_stations=n_stn,
                   to_degree=to_degree)

    def __str__(self) -> str:
        if self._raw is None:
            return "Phase(<empty>)"
        unit = "°" if self._to_degree else "mrad"
        span = f"{self.v_min:g}–{self.v_max:g} {unit}"
        return f"Phase(n={len(self)}, span={span})"
    
class Resistivity(VariationBase):
    """
    Apparent-resistivity (ρₐ) container – Ohm·m.

    The class stores the *raw* field estimate produced by Zonge’s AVG as
    well as, optionally, the **SRes** column generated by the ASTATIC
    utility.  Both arrays share the same *(n_freq × n_stations)* layout
    and are exposed through the :pyattr:`loc` / :pyattr:`sres` maps.

    Parameters
    ----------
    data : Sequence[float] | np.ndarray | None
        Flat vector, one value per frequency–station pair.
    n_freq, n_stations : int | None, optional
        Grid dimensions.  Required the first time :py:meth:`read`
        (or :py:meth:`set_sres`) is called.
    sres : Sequence[float] | np.ndarray | None, optional
        **Astatic-corrected** ρₐ produced by Zonge’s *ASTATIC* program.
        Must share the same length / ordering as *data*.
    **kwargs
        Forwarded to :class:`VariationBase` (currently unused).

    Attributes
    ----------
    values : np.ndarray
        Unique sorted resistivity values (Ohm·m).
    v_min, v_max : float
        Global extrema of *data*.
    loc : dict[str, np.ndarray]
        ``{station-id: ρₐ-vector}``.
    sres : dict[str, np.ndarray] | None
        Same mapping for the *SRes* column (``None`` if absent).
    """

    def __init__(
            self,
            data:       Sequence[float] | np.ndarray | None = None,
            *,
            n_freq:     int | None = None,
            n_stations: int | None = None,
            sres:       Sequence[float] | np.ndarray | None = None,
            **kwargs:   Any,
        ) -> None:

        super().__init__(
            data, n_freq=n_freq, n_stations=n_stations,
            **kwargs
        )

        self._sres_raw: np.ndarray | None          = None
        self.sres     : Dict[str, np.ndarray] | None = None

        if sres is not None:
            self.set_sres(sres, n_freq=n_freq, n_stations=n_stations)

    def set_sres(
        self,
        sres:       Sequence[float] | np.ndarray,
        *,
        n_freq:     int | None = None,
        n_stations: int | None = None,
        ) -> None:
        """
        Attach an *ASTATIC* resistivity column.

        Raises
        ------
        ResistivityError
            When the supplied vector length is inconsistent with the
            current *(n_freq, n_stations)* grid.
        """
        arr = np.asarray(sres, dtype=float).ravel()

        nf = n_freq     if n_freq     is not None else self.n_freq
        ns = n_stations if n_stations is not None else self.n_stations
        if nf is None or ns is None:
            raise ResistivityError("`n_freq` and `n_stations` are required")
        if arr.size != nf * ns:
            raise ResistivityError(
                "SRes length mismatch "
                f"(got {arr.size}, expected {nf * ns})")

        self._sres_raw = arr
        ids, _         = number_stations(ns, nf)
        chunks         = chunk_by_frequency(arr, nf)
        self.sres      = dict(zip(ids, chunks))


    def write(self) -> List[str]:
        """Serialise as simple ``rho=<val>`` (and ``sres=<val>``) rows."""
        rows = [f"rho={v:g}" for v in (self._raw or [])]
        if self._sres_raw is not None:
            rows += [f"sres={v:g}" for v in self._sres_raw]
        return rows

    def as_tensor(self) -> "ImpedanceTensor":
        """
        Pack the current ρₐ grid into an :class:`~ImpedanceTensor`.

        Only the *rho* block is filled; other stacks are left ``None``.
        """
        if self._raw is None:
            raise ResistivityError("Nothing loaded – call `read()` first")
        if self.n_freq is None or self.n_stations is None:
            raise ResistivityError("`n_freq` / `n_stations` unknown")

        rho_fs = self._raw.reshape(self.n_freq, self.n_stations)

        rho_dict = {comp: rho_fs for comp in ("xx", "xy", "yx", "yy")}

        return TensorFactory.build(rho=rho_dict)


    @classmethod
    def from_tensor(cls,
                    tensor: "ImpedanceTensor") -> "Resistivity":
        """
        Create a :class:`Resistivity` from an :class:`ImpedanceTensor`
        that already holds a ρₐ block.

        Only the canonical *rho* slab is extracted; *SRes* is not
        populated here – call :py:meth:`set_sres` afterwards if needed.
        """
        if tensor.rho is None:
            raise ResistivityError("Input tensor carries no rho block")

        n_freq, n_stn = tensor.rho.shape[:2]
        data = tensor.rho[:, :, 0, 0].ravel()   # flatten xx-component

        return cls(data,
                   n_freq=n_freq,
                   n_stations=n_stn)

    def __str__(self) -> str:
        span = f"{self.v_min:g}–{self.v_max:g} Ω·m" if self.v_min else "n/a"
        suffix = ", +SRes" if self._sres_raw is not None else ""
        return (f"Resistivity(n={len(self)}, span={span}, "
                f"stations={self.n_stations}, freqs={self.n_freq}{suffix})")

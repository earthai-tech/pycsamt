# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Scalar apparent-field estimates: Resistivity and Phase.

Both classes inherit from :class:`TensorBase` so their values can
be rearranged into 2×2 tensor-like grids per frequency (and per
station), according to the component label found in the table
(e.g., ``ExHy`` goes to slot (0, 1)).

They accept *modern* CSAVGW-like tables and common *legacy*
aliases (AMTAVG/MTEdit).  Minimal metadata (units) is preserved
and written back in a compact CSV block.

Usage
-----

>>> r = Resistivity.from_avg((df, meta))
>>> Z, f, st = r.to_tensor(var="rho", align="union")
>>> ds = r.to_xarray()

>>> p = Phase.from_avg((df, meta))
>>> p.convert_unit("deg")
>>> Zφ, f, st = p.to_tensor(var="phase")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ..exceptions import AvgDataError
from .tensor import TensorBase
from .utils import ( 
    _find_col, 
    _to_num, 
    _norm_comp, 
    to_xarray as _to_xr, 
)

__all__ = ["Resistivity", "Phase"]


@dataclass
class Resistivity(TensorBase):
    r"""
    Apparent resistivity (:math:`\rho_a`) per component.

    This is a *scalar* variable (one number per row), but the
    class inherits the tensor-alignment logic so values land in
    the appropriate :math:`2\times 2` slot given by ``comp``.

    Recognized source columns (legacy + modern)
    -------------------------------------------
    * ``ARes.mag`` (modern)
    * ``Resistivity`` (legacy tables)
    * ``rho`` / ``Rho`` / ``ares.mag`` (variants)

    Canonical internal column: ``'rho'``.

    Attributes
    ----------
    VAR_NAME : str
        Canonical column in :pyattr:`frame` (``"rho"``).
    ALIASES : tuple[str, ...]
        Ordered alias list used during :meth:`read`.
    UNIT_ATTR : str
        Dataset attribute used for units (``"Unit.Rho"``).
    """

    VAR_NAME: str = field(init=False, default="rho")
    ALIASES: Tuple[str, ...] = field(
        init=False,
        default=(
            "ARes.mag", "Resistivity", "rho", "Rho", "ares.mag",
            "ARes", "ARes_mag",
        ),
    )
    UNIT_ATTR: str = field(init=False, default="Unit.Rho")

    # API 

    def read(  # noqa: D401
        self,
        source: Union[
            pd.DataFrame, Sequence[float], np.ndarray, pd.Series
        ],
        meta: Optional[Mapping[str, Any]] = None,
        **kws: Any,
    ) -> None:
        """
        Parse *source* and produce a compact frame with
        ``station, freq, comp, rho``.  Coordinates are injected
        conservatively when absent (``station=NaN``,
        ``comp='ExHy'``).
        """
        self._meta = dict(meta or {})
        # ensure a stable unit hint
        self._meta.setdefault(self.UNIT_ATTR, "ohm·m")

        # vector-like → create a minimal tidy frame
        if isinstance(source, (list, tuple, np.ndarray, pd.Series)):
            vec = pd.to_numeric(pd.Series(source), errors="coerce")
            df = pd.DataFrame({self.VAR_NAME: vec})
            df["station"] = kws.get("station", np.nan)
            df["freq"] = kws.get("freq", np.nan)
            df["comp"] = _norm_comp(kws.get("comp", "ExHy"))
            self._frame = df[["station", "freq", "comp", self.VAR_NAME]]
            return

        if not isinstance(source, pd.DataFrame):
            raise TypeError(
                "Resistivity.read expects DataFrame or vector-like."
            )

        df = source.copy()
        col = _find_col(df, self.ALIASES)
        if col is None:
            raise AvgDataError(
                "Resistivity: none of aliases "
                f"{self.ALIASES!r} present in table."
            )

        df = df.rename(columns={col: self.VAR_NAME})

        # coords (inject defaults)
        if "comp" not in df.columns:
            df["comp"] = "ExHy"
        if "station" not in df.columns:
            df["station"] = np.nan
        if "freq" not in df.columns:
            raise AvgDataError("Resistivity requires a 'freq' column.")

        # normalize dtypes
        df["station"] = pd.to_numeric(df["station"], errors="ignore")
        df["freq"] = pd.to_numeric(df["freq"], errors="coerce")
        df["comp"] = df["comp"].map(_norm_comp)
        df[self.VAR_NAME] = df[self.VAR_NAME].map(_to_num)

        self._frame = df[["station", "freq", "comp", self.VAR_NAME]]

    def write(
        self,
        *,
        float_fmt: str = "%.6g",
        na_rep: str = "",
    ) -> Sequence[str]:
        """
        Serialise to a compact CSV block with a small meta
        preamble (units + timestamp).
        """
        meta = {self.UNIT_ATTR: self._meta.get(self.UNIT_ATTR, "ohm·m")}
        tmp = self.__class__()  # ephemeral for helper reuse
        tmp._frame = self._frame.copy()
        tmp._meta = meta
        return tmp._write_csv_block(
            cols=["station", "freq", "comp", self.VAR_NAME],
            title=r"$Resistivity Block",
            float_fmt=float_fmt,
            na_rep=na_rep,
            include_meta=True,
            stamp=True,
        )

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq", "comp"),
        attrs: Optional[Mapping[str, Any]] = None,
    ):
        """
        Convert to an :class:`xarray.Dataset` with one data
        variable (``rho``).  Attributes include a stable unit
        hint under ``Unit.Rho``.
        """
        merged = dict(self._meta)
        merged.setdefault(self.UNIT_ATTR, "ohm·m")
        if attrs:
            merged.update(attrs)
        return _to_xr(
            self._frame.copy(),
            coords=coords,
            data_vars=[self.VAR_NAME],
            attrs=merged,
        )

@dataclass
class Phase(TensorBase):
    r"""
    Impedance phase (:math:`\varphi`) per component.

    Values are typically reported in **milliradians** for modern
    CSAVGW tables (``Unit.Phase='mrad'``).  Use
    :meth:`convert_unit` to switch to degrees when desired.

    Recognized source columns (legacy + modern)
    -------------------------------------------
    * ``Z.phz`` (modern)
    * ``Phase`` (legacy)
    * ``phase`` / ``ZPHZ`` (variants)

    Canonical internal column: ``'phase'``.
    """

    VAR_NAME: str = field(init=False, default="phase")
    ALIASES: Tuple[str, ...] = field(
        init=False,
        default=("Z.phz", "z.phz", "Phase", "phase"),
    )
    UNIT_ATTR: str = field(init=False, default="Unit.Phase")

    # API #
    def read(  # noqa: D401
        self,
        source: Union[
            pd.DataFrame, Sequence[float], np.ndarray, pd.Series
        ],
        meta: Optional[Mapping[str, Any]] = None,
        **kws: Any,
    ) -> None:
        """
        Parse *source* and produce a compact frame with
        ``station, freq, comp, phase``.  Coordinates are injected
        conservatively when absent (``station=NaN``,
        ``comp='ExHy'``).  ``Unit.Phase`` defaults to ``'mrad'``.
        """
        self._meta = dict(meta or {})
        self._meta.setdefault(self.UNIT_ATTR, "mrad")

        # vector-like path
        if isinstance(source, (list, tuple, np.ndarray, pd.Series)):
            vec = pd.to_numeric(pd.Series(source), errors="coerce")
            df = pd.DataFrame({self.VAR_NAME: vec})
            df["station"] = kws.get("station", np.nan)
            df["freq"] = kws.get("freq", np.nan)
            df["comp"] = _norm_comp(kws.get("comp", "ExHy"))
            self._frame = df[["station", "freq", "comp", self.VAR_NAME]]
            return

        if not isinstance(source, pd.DataFrame):
            raise TypeError(
                "Phase.read expects DataFrame or vector-like."
            )

        df = source.copy()
        col = _find_col(df, self.ALIASES)
        if col is None:
            raise AvgDataError(
                "Phase: none of aliases "
                f"{self.ALIASES!r} present in table."
            )

        df = df.rename(columns={col: self.VAR_NAME})

        if "comp" not in df.columns:
            df["comp"] = "ExHy"
        if "station" not in df.columns:
            df["station"] = np.nan
        if "freq" not in df.columns:
            raise AvgDataError("Phase requires a 'freq' column.")

        df["station"] = pd.to_numeric(df["station"], errors="ignore")
        df["freq"] = pd.to_numeric(df["freq"], errors="coerce")
        df["comp"] = df["comp"].map(_norm_comp)
        df[self.VAR_NAME] = df[self.VAR_NAME].map(_to_num)

        self._frame = df[["station", "freq", "comp", self.VAR_NAME]]

    def write(
        self,
        *,
        float_fmt: str = "%.6g",
        na_rep: str = "",
    ) -> Sequence[str]:
        """
        Serialise to a compact CSV block with a meta preamble
        carrying ``Unit.Phase`` and a UTC timestamp.
        """
        meta = {self.UNIT_ATTR: self._meta.get(self.UNIT_ATTR, "mrad")}
        tmp = self.__class__()  # ephemeral
        tmp._frame = self._frame.copy()
        tmp._meta = meta
        return tmp._write_csv_block(
            cols=["station", "freq", "comp", self.VAR_NAME],
            title=r"$Phase Block",
            float_fmt=float_fmt,
            na_rep=na_rep,
            include_meta=True,
            stamp=True,
        )

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq", "comp"),
        attrs: Optional[Mapping[str, Any]] = None,
    ):
        """
        Convert to an :class:`xarray.Dataset` with one data
        variable (``phase``).  Attributes include a stable unit
        hint under ``Unit.Phase``.
        """
        merged = dict(self._meta)
        merged.setdefault(self.UNIT_ATTR, "mrad")
        if attrs:
            merged.update(attrs)
        return _to_xr(
            self._frame.copy(),
            coords=coords,
            data_vars=[self.VAR_NAME],
            attrs=merged,
        )

    def convert_unit(self, target: str = "mrad") -> None:
        r"""
        Convert the phase values **in place** between
        :math:`\mathrm{mrad}` and :math:`\mathrm{deg}` and update
        ``Unit.Phase`` in :pyattr:`meta`.

        .. math::
           1~\mathrm{mrad} =
           \frac{180}{\pi \cdot 1000}~\mathrm{deg}
        """
        cur = str(self._meta.get(self.UNIT_ATTR, "mrad")).lower()
        tgt = str(target).lower()
        if cur == tgt:
            return
        if tgt not in {"mrad", "deg"}:
            raise ValueError("target must be 'mrad' or 'deg'")

        col = self.VAR_NAME
        if col not in self._frame.columns:
            return

        x = pd.to_numeric(self._frame[col], errors="coerce")

        if cur == "mrad" and tgt == "deg":
            factor = 180.0 / (np.pi * 1000.0)
        elif cur == "deg" and tgt == "mrad":
            factor = (np.pi / 180.0) * 1000.0
        else:
            raise ValueError(f"unsupported conversion {cur} → {tgt}")

        self._frame[col] = x * factor
        self._meta[self.UNIT_ATTR] = tgt

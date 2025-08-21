# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
TensorBase – generic 2×2 impedance-like tensor adapter.

This mixin-like base sits on top of AVGComponentBase. It provides
helpers to move between tidy per-component tables and 3D/4D tensor
blocks suitable for numerical work:

    (freq, 2, 2)                – single-station view
    (station, freq, 2, 2)       – multi-station view

Supported component labels:
    - MT style: Zxx, Zxy, Zyx, Zyy
    - CSAMT style: ExHx, ExHy, EyHx, EyHy
"""

from __future__ import annotations

from typing import ( 
    Any, 
    Dict, 
    Iterable, 
    Mapping, 
    Optional, 
    Sequence, 
    Tuple, 
    Union
)
import numpy as np
import pandas as pd

from ..utils.deps import ensure_pkg 
from ..exceptions import AvgDataError, TensorError
from .base import AVGComponentBase

# --------------------------------------------------------------------------- #
# Component <-> matrix position maps
#   row axis:  E-field: Ex (0), Ey (1)
#   col axis:  H-field: Hx (0), Hy (1)
# --------------------------------------------------------------------------- #

__all__ = ["TensorBase"]

_COMP_POS: Dict[str, Tuple[int, int]] = {
    # MT naming
    "ZXX": (0, 0), "ZXY": (0, 1),
    "ZYX": (1, 0), "ZYY": (1, 1),

    # CSAMT naming
    "EXHX": (0, 0), "EXHY": (0, 1),
    "EYHX": (1, 0), "EYHY": (1, 1),
}

_E_AXIS = np.array(["Ex", "Ey"])
_H_AXIS = np.array(["Hx", "Hy"])


class TensorBase(AVGComponentBase):
    r"""Adds impedance-like tensor helpers to a component.

    This class acts as a mixin, providing methods to transform
    data between a tidy DataFrame format (one measurement per
    row) and a dense, multi-dimensional tensor format suitable
    for numerical computations.

    It is agnostic to the actual physical quantity being
    reshaped; the `var` parameter in its methods specifies which
    column from the internal DataFrame to use for the tensor's
    values.

    Notes
    -----
    Subclasses must provide a tidy `_frame` attribute containing
    at least the columns ``['freq', 'comp']`` and optionally
    ``'station'``.

    The tensor axes are consistently ordered:
    - 3D (single station): ``(frequency, E-field, H-field)``
    - 4D (multi-station): ``(station, frequency, E, H)``

    The E-field and H-field axes are of size 2, corresponding to
    the x and y components.

    Methods
    -------
    to_tensor(var, station=None, ...)
        Converts a data column into a NumPy ndarray with a shape
        of ``(..., 2, 2)``.
    from_tensor(tensor, freqs, var, stations=None, ...)
        Reconstructs a tidy DataFrame from a NumPy tensor.
    to_xarray_tensor(var, station=None, ...)
        Converts a data column into a labeled `xarray.DataArray`.

    See Also
    --------
    Z : A key subclass that uses these tensor operations.
    Resistivity : Another subclass that benefits from this mixin.
    Phase : A third subclass that uses this mixin.
    """

    @staticmethod
    def _ensure_columns(df: pd.DataFrame) -> None:
        need_any = {"freq", "comp"}
        missing = [c for c in need_any if c not in df.columns]
        if missing:
            raise AvgDataError(f"missing required columns: {missing}")

    @staticmethod
    def _prepare_table(
        df: pd.DataFrame,
        *,
        var: str,
        agg: str | None = "mean",
    ) -> pd.DataFrame:
        """
        Validate and reduce duplicates for (station,freq,comp).
        """
        TensorBase._ensure_columns(df)
        if var not in df.columns:
            raise AvgDataError(
                f"column '{var}' not found; available: {list(df.columns)}"
            )
        work = df.copy()

        # Normalize comp tokens and drop rows with unknown comps
        work["__comp_norm__"] = work["comp"].map(_norm_comp)
        work = work[work["__comp_norm__"].notna()].copy()

        # Best-effort numeric coercion for freq/station
        work["freq"] = pd.to_numeric(work["freq"], errors="coerce")
        if "station" in work.columns:
            work["station"] = pd.to_numeric(
                work["station"], errors="coerce")
        # groupby keys present in table
        keys = ["freq", "__comp_norm__"]
        if "station" in work.columns:
            keys = ["station"] + keys

        # Aggregate duplicates if needed
        if agg:
            gb = work.groupby(keys, sort=True, dropna=False)
            work = gb[var].agg(agg).to_frame(var).reset_index()
        return work

    def to_tensor(
        self,
        *,
        var: str,
        station: Optional[Union[int, float]] = None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        sort_freq: bool = True,
        align: str = "union",
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Convert per-component values into a 2×2 tensor per frequency.

        Returns
        -------
        tensor : np.ndarray
            If *station* is provided → shape (n_freq, 2, 2).
            Else (multi-station) → shape (n_station, n_freq, 2, 2).
        freqs : np.ndarray
            Sorted unique frequencies used in the tensor grid.
        stations : np.ndarray
            Stations used (size 0 for single-station request).
        """
        df = self._frame
        work = self._prepare_table(df, var=var, agg=agg)

        # Figure out station axis
        if station is not None or "station" not in work.columns:
            # Single-station path
            if station is not None:
                # tolerant numeric compare
                st_mask = np.isclose(
                    pd.to_numeric(work.get("station", np.nan), errors="coerce"),
                    float(station), equal_nan=False,
                ) if "station" in work.columns else np.ones(len(work), bool)
                ws = work.loc[st_mask]
            else:
                ws = work

            freqs = np.unique(ws["freq"].to_numpy())
            if sort_freq:
                freqs = np.sort(freqs)

            T = np.full((freqs.size, 2, 2), fill_value, dtype=float)

            # Fill
            if not ws.empty:
                # Map comp → (i,j)
                for _, row in ws.iterrows():
                    f = row["freq"]
                    c = row["__comp_norm__"]
                    val = row[var]
                    if pd.isna(val) or c is None or pd.isna(f):
                        continue
                    i_f = int(np.searchsorted(freqs, f))
                    i, j = _COMP_POS[c]
                    T[i_f, i, j] = float(val)

            return T, freqs, np.array([])

        # Multi-station path
        stations = _station_array(work["station"].unique())
        # union vs intersection of frequencies across stations
        if align not in {"union", "intersection"}:
            raise ValueError(
                "align must be 'union' or 'intersection'")

        if align == "union":
            freqs = np.unique(work["freq"].to_numpy())
        else:
            # intersection
            freqs = None
            for st in stations:
                mask = np.isclose(
                    pd.to_numeric(work["station"], errors="coerce"),
                    float(st),
                    equal_nan=False,
                )
                f_st = np.unique(work.loc[mask, "freq"].to_numpy())
                freqs = f_st if freqs is None else np.intersect1d(freqs, f_st)
            if freqs is None:
                freqs = np.array([])

        if sort_freq:
            freqs = np.sort(freqs)
        
        T = np.full((stations.size, freqs.size, 2, 2), 
                    fill_value, dtype=float)
        
        for si, st in enumerate(stations):
            mask = np.isclose(
                pd.to_numeric(work["station"], errors="coerce"), 
                float(st), 
                equal_nan=False
            )
            ws = work.loc[mask]
        
            for _, row in ws.iterrows():
                f = row["freq"]
                c = row["__comp_norm__"]
                val = row[var]
                if pd.isna(val) or c is None or pd.isna(f):
                    continue
        
                # find insertion point in the (sorted) freqs grid
                fi = int(np.searchsorted(freqs, f))
                # guard: if f is not exactly on 
                # the grid, skip (intersection case)
                if fi >= freqs.size or freqs[fi] != f:
                    continue
        
                i, j = _COMP_POS[c]
                T[si, fi, i, j] = float(val)

        return T, freqs, stations

    @staticmethod
    def from_tensor(
        tensor: np.ndarray,
        freqs: Sequence[float],
        *,
        var: str,
        stations: Optional[Sequence[Union[int, float, str]]] = None,
        comp_style: str = "mt",
    ) -> pd.DataFrame:
        """
        Reconstruct a tidy frame from a (…×2×2) tensor.

        Parameters
        ----------
        tensor
            Either (n_freq, 2, 2) or (n_station, n_freq, 2, 2).
        freqs
            Frequencies corresponding to axis 0 (or 1).
        var
            Column name to emit for the tensor values.
        stations
            If provided and tensor is 4-D, labels for station axis.
        comp_style
            'mt' → Zxx/Zxy/Zyx/Zyy ; 'csamt' → ExHx/ExHy/EyHx/EyHy
        """
        arr = np.asarray(tensor)
        if arr.ndim == 3:
            # (n_freq, 2, 2) → single-station
            s_axis = None
            f_axis = 0
        elif arr.ndim == 4:
            s_axis, f_axis = 0, 1 # noqa
        else:
            raise TensorError("tensor must be 3D or 4D")

        if arr.shape[-2:] != (2, 2):
            raise TensorError("last two dims must be (2,2)")

        # Choose component label set
        if comp_style.lower().startswith("mt"):
            comps = np.array(["Zxx", "Zxy", "Zyx", "Zyy"])
        else:
            comps = np.array(["ExHx", "ExHy", "EyHx", "EyHy"])

        # Build rows
        rows = []
        if s_axis is None:
            for fi, f in enumerate(freqs):
                block = arr[fi]
                vals = [block[0, 0], block[0, 1], block[1, 0], block[1, 1]]
                for comp, val in zip(comps, vals):
                    rows.append(
                        {"station": np.nan, "freq": float(f),
                         "comp": comp, var: float(val)}
                    )
        else:
            if stations is None:
                stations = list(range(arr.shape[0]))
            for si, st in enumerate(stations):
                for fi, f in enumerate(freqs):
                    block = arr[si, fi]
                    vals = [block[0, 0], block[0, 1], block[1, 0], block[1, 1]]
                    for comp, val in zip(comps, vals):
                        rows.append(
                            {"station": st, "freq": float(f),
                             "comp": comp, var: float(val)}
                        )

        return pd.DataFrame.from_records(rows)

    @ensure_pkg (
        'xarray', 
        extra="xarray is required for to_xarray_tensor()"
    )
    def to_xarray_tensor(
        self,
        *,
        var: str,
        station: Optional[Union[int, float]] = None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        attrs: Optional[Mapping[str, Any]] = None,
    ):
        """
        Return a 3-D or 4-D xarray.DataArray with dims:
            single-station   → (freq, e, h)
            multi-station    → (station, freq, e, h)
        """
        import xarray as xr  # type: ignore

        T, freqs, stations = self.to_tensor(
            var=var, station=station, agg=agg, fill_value=fill_value
        )

        # coords
        e = _E_AXIS
        h = _H_AXIS
        if stations.size == 0:
            da = xr.DataArray(
                T, dims=("freq", "e", "h"),
                coords={"freq": freqs, "e": e, "h": h},
                attrs=dict(attrs or {}),
                name=var,
            )
        else:
            da = xr.DataArray(
                T, dims=("station", "freq", "e", "h"),
                coords={"station": stations, "freq": freqs,
                        "e": e, "h": h},
                attrs=dict(attrs or {}),
                name=var,
            )
        return da

    def read(
        self,
        source: pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
        **kws: Any,
    ) -> None:
        if not isinstance(source, pd.DataFrame):
            raise TypeError("TensorBase.read expects a DataFrame.")
        df = source.copy()

        # normalise component label if present
        if "comp" in df.columns:
            df["comp"] = df["comp"].map(lambda c: _norm_comp(c))

        self._frame = df
        self._meta = dict(meta or {})

    def write(self) -> Sequence[str]:
        # Minimal, mostly for debugging; not used by tests
        if self._frame.empty:
            return ["\\ $_TensorBase", ""]
        return self._write_csv_block(
            cols=list(self._frame.columns),
            title="$_TensorBase",
            include_meta=False,
            stamp=False,
        )

    def __str__(self) -> str:
        r, c = self.shape
        cols = ", ".join(self._frame.columns[:6])
        tail = "…" if self._frame.shape[1] > 6 else ""
        return f"TensorBase[{r}×{c}] cols=[{cols}{tail}]"

    __repr__ = __str__


def _norm_comp(label: Any) -> Optional[str]:
    """
    Normalize a component label (‘ExHy’, ‘Zxy’, etc.) to an
    uppercase token present in _COMP_POS. Return None if unknown.
    """
    if label is None:
        return None
    s = str(label).strip().upper()
    # Fast path
    if s in _COMP_POS:
        return s
    # Try to strip non-alnum characters just in case
    s2 = "".join(ch for ch in s if ch.isalnum())
    return s2 if s2 in _COMP_POS else None

def _station_array(values: Iterable[Any]) -> np.ndarray:
    """
    Normalize station coordinate array for indexing. Keep numeric
    if possible, else fall back to strings.
    """
    vals = pd.Series(values)
    num = pd.to_numeric(vals, errors="coerce")
    if num.notna().all():
        return num.to_numpy()
    return vals.astype(str).to_numpy()

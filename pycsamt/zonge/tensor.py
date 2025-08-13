# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com> (structure & package)
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
from dataclasses import dataclass #, field
import numpy as np
import pandas as pd

from .base import AVGComponentBase

from ..utils.deps import ensure_pkg 
from ..exceptions import AvgDataError
# from ..constants import PI   
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


class TensorBase(AVGComponentBase):
    """
    Add impedance-like tensor helpers on top of AVGComponentBase.

    Subclasses must provide a tidy frame with at least:
        ['freq', 'comp'] and (optionally) 'station'.

    You can then call:
        to_tensor(var='rho', station=...)        # → np.ndarray
        from_tensor(tensor, stations, freqs, ...)# → DataFrame
        to_xarray_tensor(var='zabs')             # → xr.DataArray

    The base is agnostic to “what” the tensor measures; `var`
    is simply the column name you want to fold into 2×2 shape.
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
            raise ValueError("align must be 'union' or 'intersection'")

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
            s_axis, f_axis = 0, 1
        else:
            raise AvgDataError("tensor must be 3D or 4D")

        if arr.shape[-2:] != (2, 2):
            raise AvgDataError("last two dims must be (2,2)")

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

class TensorBuildError (Exception):
    """ Tensor base exceptions"""
    pass 

@dataclass(slots=True)
class ImpedanceTensor:
    """Container holding *all* 2×2 stacks for a station / site.

    Attributes
    ----------
    z : ndarray, complex
        Shape ``(n_freq, 2, 2)``.
    rho, phase : ndarray, float
        Same shape as *z*.  Phase is **degrees**.
    z_err, rho_err, phase_err : ndarray | None
        Optional 1-σ uncertainties.  When *None* no error-prop
        was provided for the corresponding field.
    """
    z:          np.ndarray
    rho:        np.ndarray
    phase:      np.ndarray
    z_err:      np.ndarray | None = None
    rho_err:    np.ndarray | None = None
    phase_err:  np.ndarray | None = None


    def _slice(self, i: int, j: int, field: str = "z") -> np.ndarray:
        return getattr(self, field)[:, i, j]

    # impedance
    @property 
    def z_xx(self): return self._slice(0, 0, "z")
    @property 
    def z_xy(self): return self._slice(0, 1, "z")
    @property 
    def z_yx(self): return self._slice(1, 0, "z")
    @property 
    def z_yy(self): return self._slice(1, 1, "z")
    # resistivity
    @property 
    def rho_xx(self): return self._slice(0, 0, "rho")
    @property 
    def rho_xy(self): return self._slice(0, 1, "rho")
    @property 
    def rho_yx(self): return self._slice(1, 0, "rho")
    @property 
    def rho_yy(self): return self._slice(1, 1, "rho")
    # phase (deg)
    @property 
    def ph_xx(self): return self._slice(0, 0, "phase")
    @property 
    def ph_xy(self): return self._slice(0, 1, "phase")
    @property 
    def ph_yx(self): return self._slice(1, 0, "phase")
    @property 
    def ph_yy(self): return self._slice(1, 1, "phase")

    def __str__(self) -> str:                       # user-friendly
        blocks = ", ".join(
            k for k, v in (("Z", self.z),
                            ("ρ", self.rho),
                            ("ϕ", self.phase)) if v is not None
        )
        nf, ns = self.z.shape[:2] if self.z is not None else \
                  self.rho.shape[:2] if self.rho is not None else \
                  self.phase.shape[:2]
        return f"ImpedanceTensor({blocks}; n_freq={nf}, n_station={ns})"

    def __repr__(self) -> str:                      # dev-oriented
        cls = self.__class__.__name__
        nf, ns = (self.z or self.rho or self.phase).shape[:2]
        flags = "/".join(flag for flag in
                         ("Z"  if self.z     is not None else "-",
                          "ρ"  if self.rho   is not None else "-",
                          "ϕ"  if self.phase is not None else "-"))
        return (f"<{cls} {flags} "
                f"(n_freq={nf}, n_station={ns}) at 0x{id(self):x}>")



class TensorFactory:
    """
    Tiny helper that stacks four *xx/xy/yx/yy* component vectors into the
    canonical ``(n_freq, 1, 2, 2)`` shape expected by
    :class:`~pycsamt.zonge.Z`.
    """

    @classmethod
    def build(
        cls,
        *,
        # any **one** (or more) of the following can be supplied
        z:        Dict[str, Sequence[Any]] | None = None,
        rho:      Dict[str, Sequence[Any]] | None = None,
        phase:    Dict[str, Sequence[Any]] | None = None,
        # optional one-sigma errors
        z_err:    Dict[str, Sequence[Any]] | None = None,
        rho_err:  Dict[str, Sequence[Any]] | None = None,
        ph_err:   Dict[str, Sequence[Any]] | None = None,
        dtype_z:  Any = complex,
        dtype_f:  Any = float,
    ) -> "ImpedanceTensor":
        """
        Assemble an :class:`ImpedanceTensor` from the components that are
        **actually** available.  Any non-supplied block (and its error)
        will be set to ``None``.

        Example
        -------
        >>> t = TensorFactory.build(rho=rho_dict)   # only ρₐ known
        >>> assert t.rho is not None and t.z is None
        """
        if z is rho is phase is None:
            raise TensorBuildError("At least one of z / rho / phase is required")

        # figure out grid size from the first non-None block
        proto = next(b for b in (z, rho, phase) if b is not None)
        cls._check_dict(proto, "prototype")
        n_freq = len(next(iter(proto.values())))

        # helper to stack or return None
        def _maybe_stack(block, kind, dtype):
            if block is None:
                return None
            cls._check_dict(block, kind)
            cls._assert_equal_len(block, n_freq, kind)
            return cls._stack(block, dtype)

        z_t     = _maybe_stack(z,        "z",     dtype_z)
        rho_t   = _maybe_stack(rho,      "rho",   dtype_f)
        phase_t = _maybe_stack(phase,    "phase", dtype_f)
        ze_t    = _maybe_stack(z_err,    "z_err", dtype_f)
        re_t    = _maybe_stack(rho_err,  "rho_err", dtype_f)
        pe_t    = _maybe_stack(ph_err,   "ph_err", dtype_f)

        return ImpedanceTensor(z_t, rho_t, phase_t, ze_t, re_t, pe_t)


    _REQ_KEYS = ("xx", "xy", "yx", "yy")

    @classmethod
    def _check_dict(cls, d: Dict[str, Sequence[Any]] | None, label: str) -> None:
        if d is None:
            return
        missing = [k for k in cls._REQ_KEYS if k not in d]
        if missing:
            raise TensorBuildError(f"{label}: missing keys {missing}")

    @staticmethod
    def _assert_equal_len(d: Dict[str, Sequence[Any]], n: int, label: str) -> None:
        for k, v in d.items():
            if len(v) != n:
                raise TensorBuildError(
                    f"{label}.{k}: length {len(v)} ≠ reference {n}"
                )

    @staticmethod
    def _stack(block: Dict[str, Sequence[Any]] | None, dtype) -> np.ndarray:
        if block is None:
            return None     # type: ignore[return-value]
        xx, xy, yx, yy = (np.asarray(block[k], dtype=dtype) for k in (
            "xx", "xy", "yx", "yy"))
        return np.stack(  # (n_freq, 2, 2)
            (np.stack((xx, xy), axis=-1),
             np.stack((yx, yy), axis=-1)),
            axis=-2
        )[..., None, :, :]        # → (n_freq, 1, 2, 2)

#  Author : LKouadio <etanoyau@gmail.com>
#  License: LGPL-3.0
"""
pycsamt.zonge.meas

Lightweight, robust containers for *measurement* columns that
appear in Zonge AVG tables.  These components follow the new
table-centric design (inherit from AVGComponentBase), so they can
be fed with a tidy DataFrame and later serialize as small CSV
blocks if desired.

Exports
-------
- CompMeas : validator/normaliser for component labels ('ExHy', …)
- Amps     : transmitter current amplitude handler (A)

Both classes keep *context* columns ('station', 'freq', 'comp')
when available, which helps downstream grouping and reshaping.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import (
    Any,
)

import numpy as np
import pandas as pd

from ..exceptions import FrequencyError, InputError
from ..utils._dependency import import_optional_dependency
from .base import AVGComponentBase
from .utils import _standardise_columns
from .utils import to_xarray as _to_xr

__all__ = ["CompMeas", "Amps", "Frequency"]


_NUMERIC_NA = {"", "*", "nan", "NaN", "None", None}


class CompMeas(AVGComponentBase):
    """
    Enumeration/validator for classical CSAMT component labels.

    The class guarantees a canonical ``comp`` column with allowed
    values (case-sensitive): ``{'ExHy','ExHx','EyHy','EyHx'}``.
    A few common case variants are accepted and normalized.

    Typical usage
    -------------
    >>> cm = CompMeas.from_avg((df, meta))
    >>> cm.unique
    ['ExHy']  # for a scalar survey with one component
    """

    # canonical forms allowed by the pipeline
    _VALID: set[str] = {
        "ExHy",
        "ExHx",
        "EyHy",
        "EyHx",
        "Zxx",
        "Zxy",
        "Zyx",
        "Zyy",
    }

    # tolerant normalisation map (upper/lower → canonical)
    _NORM: dict[str, str] = {
        "EXHY": "ExHy",
        "EXHX": "ExHx",
        "EYHY": "EyHy",
        "EYHX": "EyHx",
        "exhy": "ExHy",
        "exhx": "ExHx",
        "eyhy": "EyHy",
        "eyhx": "EyHx",
        "ZXX": "Zxx",
        "ZXY": "Zxy",
        "ZYX": "Zyx",
        "ZYY": "Zyy",
        "zxx": "Zxx",
        "zxy": "Zxy",
        "zyx": "Zyx",
        "zyy": "Zyy",
    }

    required: set[str] = set()  # flexible on input
    provides: set[str] = {"comp"}  # always provides comp

    def read(
        self, source: pd.DataFrame, meta: Mapping[str, Any] | None = None
    ) -> None:
        """
        Ensure a normalised ``comp`` column exists.

        - If missing, default to **'ExHy'** (legacy kind-1 style).
        - If present, normalise values and validate membership.
        """
        # start with a context copy to keep coords if present
        df = _keep_context(source, cols=("comp",))
        if "comp" not in df.columns:
            df["comp"] = "ExHy"

        # normalise tolerant variants → canonical
        def _norm_one(v: Any) -> str:
            s = str(v).strip()
            # Use .upper() for case-insensitivity before mapping
            return self._NORM.get(s.upper(), s)

        df["comp"] = df["comp"].map(_norm_one)

        # Drop only rows whose component is missing ("nan" strings can
        # be introduced during xarray conversion). Filtering *all*
        # invalid labels here made the validation below unreachable.
        df = df[
            df["comp"].notna() & (df["comp"].astype(str).str.lower() != "nan")
        ].copy()

        # validate – bail early with a friendly message
        bad = sorted(set(df["comp"]) - self._VALID)
        if bad:
            raise InputError(
                "unrecognised component(s) "
                f"{bad} – expect one of {sorted(self._VALID)}"
            )

        self._frame = df
        self._meta = dict(meta or {})

    def write(self) -> Sequence[str]:
        """
        Serialise the component block as a compact CSV fragment.

        Only the *context* columns and ``comp`` are emitted to keep
        files readable.  Disk I/O is the caller's responsibility.
        """
        cols = [c for c in ("station", "freq", "comp") if c in self._frame]
        if not cols:
            return []
        return self._write_csv_block(
            cols=cols,
            title="Component labels",
            include_meta=False,
            stamp=False,
        )

    @property
    def unique(self) -> list[str]:
        """Sorted unique component labels."""
        if "comp" not in self._frame:
            return []
        return sorted(pd.unique(self._frame["comp"]).tolist())

    def __str__(self) -> str:
        kinds = ",".join(self.unique[:4])
        more = "…" if len(self.unique) > 4 else ""
        return f"CompMeas[{len(self._frame)}] {{{kinds}{more}}}"


@dataclass
class _AmpStats:
    """Small stats tuple for quick diagnostics."""

    vmin: float | None = None
    vmax: float | None = None
    mean: float | None = None
    median: float | None = None
    count: int = 0


class Amps(AVGComponentBase):
    """
    Transmitter **current amplitude** container (unit: **A**).

    The class normalises the ``amps`` column to numeric, tracks
    a few useful stats, and can (optionally) export itself as a
    small CSV block.  Context columns are preserved.

    Examples
    --------
    >>> amps = Amps.from_avg((df, meta))
    >>> amps.stats.mean
    12.7
    >>> ds = amps.to_xarray()   # optional grid for convenience
    """

    required: set[str] = set()  # tolerant – legacy files
    provides: set[str] = {"amps"}

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        name: str | None = None,
    ) -> None:
        super().__init__(data=data, meta=meta, name=name)
        self._stats = _AmpStats()

    def read(
        self, source: pd.DataFrame, meta: Mapping[str, Any] | None = None
    ) -> None:
        """
        Parse *source* and populate the ``amps`` column as float.

        Non-numeric entries ('*', blanks) become NaN.  Context
        columns (station/freq/comp) are kept when present.
        """
        df = _keep_context(source, cols=("amps", "Tx.Amp"))
        # tolerate both modern 'amps' and legacy 'Tx.Amp'
        if "amps" not in df.columns and "Tx.Amp" in df.columns:
            df["amps"] = df["Tx.Amp"]

        if "amps" not in df.columns:
            # provide the column – some natural-source datasets
            # deliberately omit Tx current; we keep NaNs.
            df["amps"] = np.nan

        df["amps"] = df["amps"].map(_to_float)

        self._frame = df
        self._meta = dict(meta or {})
        self._compute_stats()

        return self

    def _compute_stats(self) -> None:
        """Compute quick stats on finite ``amps`` values."""
        s = pd.to_numeric(
            self._frame.get("amps", pd.Series(dtype=float)), errors="coerce"
        )

        s = s[np.isfinite(s.values)]
        if s.empty:
            self._stats = _AmpStats(count=0)
            return
        self._stats = _AmpStats(
            vmin=float(np.min(s)),
            vmax=float(np.max(s)),
            mean=float(np.mean(s)),
            median=float(np.median(s)),
            count=int(s.size),
        )

    @property
    def stats(self) -> _AmpStats:
        """Return min/max/mean/median/count snapshot."""
        return self._stats

    def as_series(self) -> pd.Series:
        """Return the ``amps`` column as a Series (copy)."""
        return self._frame.get("amps", pd.Series(dtype=float)).copy()

    def to_frame(self) -> pd.DataFrame:  # override for clarity
        """
        Return a small table with context + ``amps`` only.
        """
        keep = [
            c for c in ("station", "freq", "comp", "amps") if c in self._frame
        ]
        return self._frame.loc[:, keep].copy()

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq", "comp"),
        attrs: dict[str, Any] | None = None,
    ):
        """
        Optional convenience: grid the column into an xarray
        dataset (dims: station × freq × comp) with a single
        data-variable called ``amps``.
        """
        df = self.to_frame()
        if df.empty:
            return None
        return _to_xr(
            df,
            coords=coords,
            data_vars=["amps"],
            attrs=attrs or {"component": "amps"},
        )

    def write(self) -> Sequence[str]:
        """
        Serialise as a compact CSV fragment.  We keep context
        columns if present so the block remains useful alone.
        """
        cols = [
            c for c in ("station", "freq", "comp", "amps") if c in self._frame
        ]
        if not cols:
            return []
        return self._write_csv_block(
            cols=cols,
            title="Tx current (amps)",
            include_meta=False,
            stamp=False,
        )

    def __str__(self) -> str:
        s = self._stats
        if s.count == 0:
            return "Amps[empty]"
        span = f"{s.vmin:g}–{s.vmax:g} A"
        return f"Amps[n={s.count}, span={span}, mean={s.mean:g}]"


class Frequency(AVGComponentBase):
    """
    Frequency axis manager (Hz) for AVG tables.

    Goals
    -----
    • Read from *legacy* and *modern* frames (column aliases handled)
    • Enforce positivity (> 0 Hz) while tolerating missing markers
    • Offer stable unique grids across stations/components
    • Provide a compact `to_xarray()` for downstream use

    Notes
    -----
    - Legacy decimals like '.5' are parsed correctly.
    - Missing values ('*', '', None) become NaN and are ignored in
      summaries. Non-positive *numeric* entries raise `FrequencyError`.
    """

    # what we require/provide as a component table
    required: set[str] = set()  # we can construct from a vector
    provides: set[str] = {"freq"}

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        meta: Mapping[str, Any] | None = None,
        *,
        name: str | None = None,
    ) -> None:
        super().__init__(data=data, meta=meta, name=name or "Frequency")
        # ensure we always carry a frequency unit tag
        self._meta.setdefault("Unit.Freq", "Hz")

    def read(
        self,
        source: pd.DataFrame | Sequence[float] | np.ndarray | pd.Series,
        meta: Mapping[str, Any] | None = None,
        **kws: Any,
    ) -> None:
        """
        Load frequency values from a tidy frame *or* a flat vector.

        If *source* is a DataFrame, we try to keep `station` / `comp`
        when present. Otherwise we inject conservative defaults:
        `station=NaN`, `comp='ExHy'`.
        """
        self._meta = dict(meta or {})
        self._meta.setdefault("Unit.Freq", "Hz")

        # vector-like → build a tiny tidy frame
        if isinstance(source, (list, tuple, np.ndarray, pd.Series)):
            vec = pd.to_numeric(
                pd.Series(source, dtype="float64"), errors="coerce"
            )
            df = pd.DataFrame({"freq": vec})
            if "station" in kws:
                df["station"] = kws["station"]
            else:
                df["station"] = np.nan
            df["comp"] = kws.get("comp", "ExHy")
            self._frame = df[["station", "freq", "comp"]]
            self._validate_positive()
            return

        # dataframe path
        if not isinstance(source, pd.DataFrame):
            raise TypeError("Frequency.read expects DataFrame or vector-like")

        df = _standardise_columns(source.copy())

        if "freq" not in df.columns:
            raise FrequencyError(
                "Canonical column 'freq' not found in table."
            )
        # tidy coords (inject when absent)
        if "comp" not in df.columns:
            df["comp"] = "ExHy"
        if "station" not in df.columns:
            df["station"] = np.nan

        # robust numeric parsing ('.5' → 0.5, '*'/' ' → NaN)
        df["freq"] = (
            df["freq"]
            .astype(str)
            .str.strip()
            .replace({"": np.nan, "*": np.nan})
        )
        df["freq"] = pd.to_numeric(df["freq"], errors="coerce")

        self._frame = df[["station", "freq", "comp"]].copy()
        self._validate_positive()

        return self

    def write(self) -> list[str]:
        """
        Emit a compact CSV block with the contextual columns that we
        manage (`station`, `freq`, `comp`), suitable for round-tripping.
        """
        if self._frame.empty:
            return []
        cols: list[str] = []
        for c in ("station", "freq", "comp"):
            if c in self._frame.columns:
                cols.append(c)
        return self._write_csv_block(
            cols=cols,
            title="$Frequency Block",
            include_meta=True,
            stamp=True,
        )

    def _validate_positive(self) -> None:
        """Raise on non-positive *numeric* frequency values."""
        s = pd.to_numeric(
            self._frame.get("freq", pd.Series(dtype=float)), errors="coerce"
        )
        bad = s.notna() & (s <= 0.0)
        if bool(bad.any()):
            n = int(bad.sum())
            raise FrequencyError(f"found {n} non-positive frequency value(s)")

    def unique(
        self,
        *,
        sort: bool = True,
        dropna: bool = True,
        rtol: float = 1e-6,
        atol: float = 1e-12,
    ) -> np.ndarray:
        """
        Unique global frequency grid with tolerance deduplication.
        """
        s = pd.to_numeric(
            self._frame.get("freq", pd.Series(dtype=float)), errors="coerce"
        )
        if dropna:
            s = s.dropna()
        arr = s.to_numpy(dtype=float)
        if arr.size == 0:
            return arr
        uniq = _unique_tol(arr, rtol=rtol, atol=atol)
        if sort:
            uniq.sort()
        return uniq

    def by_station(
        self, *, rtol: float = 1e-6, atol: float = 1e-12
    ) -> dict[float, np.ndarray]:
        """
        Per-station unique frequency grids (sorted).
        """
        out: dict[float, np.ndarray] = {}
        if self._frame.empty:
            return out
        tmp = self._frame.copy()
        tmp["station"] = pd.to_numeric(tmp["station"], errors="coerce")
        for stn, sub in tmp.groupby("station", dropna=True):
            s = pd.to_numeric(sub["freq"], errors="coerce").dropna()
            out[float(stn)] = _unique_tol(s.to_numpy(), rtol=rtol, atol=atol)
        return out

    @property
    def n_unique(self) -> int:
        """Number of unique frequencies in the table."""
        return int(self.unique().size)

    @staticmethod
    def logspace(
        decade_start: int, decade_stop: int, n_points: int
    ) -> np.ndarray:
        """
        Canonical log-spaced grid (10**start → 10**stop), inclusive.
        """
        if n_points < 2:
            raise ValueError("n_points must be >= 2")
        return np.logspace(
            decade_start, decade_stop, n_points, endpoint=True, base=10.0
        )

    def to_xarray(
        self,
        *,
        coords: Sequence[str] = ("station", "freq", "comp"),
        attrs: dict | None = None,
    ):
        """
        Convert to an xarray.Dataset. We include `freq` as a data
        variable as well, which is helpful in some consumers; the
        coordinate is still the same `freq` axis created by `coords`.
        """
        import_optional_dependency(
            "xarray",
            extra="xarray is required for to_xarray()",
            errors="raise",
        )

        import xarray as xr

        # 1) Build a base dataset that has the dims we need
        df = self._frame.copy()
        df["present"] = True  # any simple data var to force a dense grid

        attrs = dict(self._meta or {})
        attrs.setdefault("Unit.Freq", "Hz")

        ds = _to_xr(
            df,
            coords=coords,
            data_vars=["present"],
            attrs=attrs,
        )

        # 2) Grab the coordinate values for broadcasting, then
        #    drop the 1-D 'freq' coord so we can use that
        # name for a data-var.
        st = ds.coords["station"].values
        cp = ds.coords["comp"].values
        fq = ds.coords["freq"].values  # 1-D list of freqs

        ds = ds.drop_vars(
            "freq"
        )  # remove the 1-D coord variable named 'freq'

        # 3) Broadcast freq values over (station, freq, comp)
        freq3 = np.broadcast_to(
            fq[np.newaxis, :, np.newaxis],
            (st.size, fq.size, cp.size),
        )
        ds["freq_grid"] = xr.DataArray(
            freq3,
            dims=("station", "freq", "comp"),
            coords={"station": st, "comp": cp},
        )

        return ds

    def __str__(self) -> str:
        if self._frame.empty:
            return "Frequency[0×0]"
        f = pd.to_numeric(self._frame["freq"], errors="coerce")
        f = f.dropna()
        if f.empty:
            return "Frequency[n=? span=?–? Hz, unique=0]"
        return (
            f"Frequency[{len(self._frame)}×{self._frame.shape[1]}] "
            f"span={f.min():g}–{f.max():g} Hz, unique={self.n_unique}"
        )

    __repr__ = __str__


def _to_float(x: Any) -> float | np.floating | np.nan:
    """
    Best-effort numeric coercion for table columns.

    - blanks / '*' / 'nan' → np.nan
    - integral strings → float(int) (we stay in float for pandas)
    - everything else → float(...) or np.nan on failure
    """
    if x in _NUMERIC_NA:
        return np.nan
    try:
        f = float(str(x).strip())
        return f
    except Exception:
        return np.nan


def _keep_context(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    """
    Return a *copy* of ``df`` with the selected columns if present,
    preserving a helpful context subset first: station/freq/comp.
    """
    want = [c for c in ("station", "freq", "comp") if c in df.columns]
    want += [c for c in cols if c not in want]
    keep = [c for c in want if c in df.columns]
    return df.loc[:, keep].copy()


def _unique_tol(
    arr: np.ndarray, *, rtol: float = 1e-6, atol: float = 1e-12
) -> np.ndarray:
    """
    Return sorted uniques using an *isclose* tolerance, which is
    useful when frequencies come from decimal strings like '.5'
    and suffer tiny rounding differences between stations.
    """
    if arr.size == 0:
        return arr
    a = np.sort(arr.astype(float, copy=False))
    keep = [a[0]]
    for x in a[1:]:
        if not np.isclose(x, keep[-1], rtol=rtol, atol=atol):
            keep.append(x)
    return np.asarray(keep, dtype=float)

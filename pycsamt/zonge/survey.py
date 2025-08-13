# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Survey layout helpers.

- Station: robust line-geometry container that understands modern/legacy
  AVG frames, exposes unique station positions, IDs, and handy header
  ($Stn.*) derivations for round-trips.
"""
from __future__ import annotations

from typing import ( 
    Any, Dict, 
    List, 
    Mapping, Optional, 
    Sequence, Tuple, 
    Union
)
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .base import AVGComponentBase # , AVGFrame
from .utils import number_stations
from ..exceptions import StationError


__all__ = ["Station"]


def _to_float_series(s: pd.Series) -> pd.Series:
    """Coerce a Series to float, preserving NaN for non-numeric."""
    return pd.to_numeric(s, errors="coerce").astype(float)


def _median_inc(vals: np.ndarray) -> Optional[float]:
    """Median increment if grid-like, else None."""
    if vals.size < 2:
        return None
    diffs = np.diff(vals)
    inc = np.median(diffs)
    if np.allclose(diffs, inc, rtol=0.0, atol=1e-6):
        return float(inc)
    return None


@dataclass(slots=True)
class Station(AVGComponentBase):
    """
    One-dimensional survey-line geometry container.

    This component reads the *station* coordinate from a tidy AVG frame
    (usually column ``'station'``), validates and stores unique positions,
    and (optionally) normalizes the origin to zero. It also derives a
    stable set of header keys ($Stn.*) for round-tripping modern headers.

    Key features
    -----------
    • Accepts either a DataFrame (from AVG parsers) or raw arrays.
    • Handles ragged frequency coverage (no equal-count requirement).
    • Tracks units; internal helpers can work in metres when desired.
    • Derives $Stn.Beg / $Stn.Inc / $Stn.Left / $Stn.Right if grid-like.
    • Exposes dictionary mappings useful for downstream selection.

    Notes
    -----
    We **do not** modify the caller's AVG frame; we keep our own view
    in ``self._frame`` that at minimum contains a numeric ``'station'``
    column and (if relevant) a ``'station_m'`` column when the input
    length unit is feet ('ft') and conversion to metres is helpful.
    """

    # configuration
    unit: str = "m"
    normalize: bool = False
    allow_ragged: bool = True

    # derived state
    values: np.ndarray = field(default_factory=lambda: np.empty(0))
    names: List[str] = field(default_factory=list)
    index_by_value: Dict[float, np.ndarray] = field(default_factory=dict)
    index_by_name: Dict[str, np.ndarray] = field(default_factory=dict)


    def read(
        self,
        source: Union[pd.DataFrame, Sequence[float], np.ndarray],
        meta: Mapping[str, Any] | None = None,
        *,
        unit: Optional[str] = None,
        names: Optional[Sequence[str]] = None,
        normalize: Optional[bool] = None,
        allow_ragged: Optional[bool] = None,
    ) -> None:
        """
        Populate the component from *source*.

        Parameters
        ----------
        source
            DataFrame with a ``'station'`` column **or** a 1-D array of
            station values (possibly repeated across frequencies).
        meta
            Optional file metadata. If provided and ``unit`` is None, we
            try ``meta['Unit.Length']`` for ('m'|'ft').
        unit
            Length unit of *source* values. Defaults to 'm' or to
            ``meta['Unit.Length']`` when available.
        names
            Optional labels for *unique* station positions. If omitted,
            we auto-generate: 'S00', 'S01', ...
        normalize
            When True, shift the origin so the first station is 0.0.
        allow_ragged
            When False, all stations must have identical row counts
            (one per frequency). When True (default), ragged coverage
            is allowed.
        """
        # resolve options
        munit = (meta or {}).get("Unit.Length") if meta else None
        unit = (unit or munit or self.unit or "m").lower()
        normalize = self.normalize if normalize is None else normalize
        allow_ragged = (
            self.allow_ragged if allow_ragged is None else allow_ragged
        )

        # build a minimal frame with a numeric 'station' column
        if isinstance(source, pd.DataFrame):
            if "station" not in source.columns:
                raise StationError("column 'station' missing in frame")
            frame = source[["station"]].copy()
            frame["station"] = _to_float_series(frame["station"])
        else:
            arr = np.asarray(source, dtype=float).ravel()
            if arr.size == 0:
                raise StationError("empty station array")
            frame = pd.DataFrame({"station": arr})

        # drop rows where station cannot be parsed (all we can do)
        frame = frame.loc[frame["station"].notna()].reset_index(drop=True)
        if frame.empty:
            raise StationError("no valid station values")

        # normalize if requested
        if normalize:
            frame["station"] = frame["station"] - frame["station"].min()

        # optionally provide a metres view if unit is feet
        if unit == "ft":
            frame["station_m"] = frame["station"] / 3.280839895
        elif unit == "km":
            frame["station_m"] = frame["station"] * 1e3
        else:
            # mirror for symmetry; callers can ignore it
            frame["station_m"] = frame["station"]

        # sort by station coord for stable introspection
        frame.sort_values("station", kind="mergesort", inplace=True)
        frame.reset_index(drop=True, inplace=True)

        # unique positions and index maps
        uniq = frame["station"].dropna().unique()
        uniq.sort()
        self.values = uniq.astype(float)

        # build group index maps
        groups = frame.groupby("station", sort=True, dropna=False).indices
        idx_by_val: Dict[float, np.ndarray] = {}
        for k, v in groups.items():
            try:
                key = float(k)
            except Exception:
                # groupby won't give NaN key here (we dropped NaNs)
                continue
            idx_by_val[key] = np.asarray(sorted(v), dtype=int)
        self.index_by_value = idx_by_val

        # ragged check if requested
        if not allow_ragged and self.values.size > 0:
            counts = np.array([len(idx_by_val[val]) for val in self.values])
            if not np.all(counts == counts[0]):
                raise StationError("inconsistent rows per station")

        # station names (one per unique position)
        if names is not None:
            if len(names) != self.values.size:
                raise StationError(
                    "`names` length must match unique station count"
                )
            self.names = list(names)
        else:
            ids, _ = number_stations(self.values.size, 1, prefix="S")
            self.names = ids

        # name→index map (same row indices as value map)
        self.index_by_name = {
            nm: self.index_by_value[val]
            for nm, val in zip(self.names, self.values, strict=True)
        }

        # stash the frame + a small meta view in the component payload
        self._frame = frame
        self._meta = dict(meta or {})
        # remember configuration for write/pretty
        self.unit = unit
        self.normalize = normalize
        self.allow_ragged = allow_ragged

    def write(self) -> Sequence[str]:
        """
        Emit a small, human-friendly header preamble plus a CSV of the
        station coordinate. This is mainly useful for diagnostics and
        tests; full file writing should be orchestrated at a higher
        level.
        """
        if self._frame.empty:
            return []

        # Try to derive canonical $Stn.* keys for the header
        keys = self.to_keywords()

        lines: List[str] = []
        lines.append(r"\ Station Geometry")
        for k, v in keys.items():
            lines.append(f"${k}={v}")
        lines.append("")  # spacer before table

        # CSV of station and metres view (if present)
        cols = ["station"] + (
            ["station_m"] if "station_m" in self._frame.columns else [])
        csv = self._frame.loc[:, cols].to_csv(index=False, float_format="%.6g")
        lines.extend(csv.splitlines())
        return lines

    @property
    def n_unique(self) -> int:
        """Number of unique station positions."""
        return int(self.values.size)

    @property
    def span(self) -> Optional[Tuple[float, float]]:
        """(min, max) of the station coordinate; None if empty."""
        if self.values.size == 0:
            return None
        return (float(self.values.min()), float(self.values.max()))

    @property
    def increment(self) -> Optional[float]:
        """Grid increment if evenly spaced, else None."""
        return _median_inc(self.values)

    def label_map(self) -> Dict[float, str]:
        """Map station numeric value → generated label (e.g., 'S03')."""
        return {float(v): n for v, n in zip(
            self.values, self.names, strict=True)}

    def to_keywords(self) -> Dict[str, Any]:
        """
        Derive a stable set of $Stn.* keys from current geometry:

        - Stn.Beg   : first station
        - Stn.Inc   : increment (when grid-like)
        - Stn.Left  : left edge (min)
        - Stn.Right : right edge (max)

        For convenience, we also mirror GDP-station versions when
        increment is grid-like: Stn.GdpBeg / Stn.GdpInc.
        """
        out: Dict[str, Any] = {}
        if self.values.size == 0:
            return out

        left, right = float(self.values.min()), float(self.values.max())
        out["Stn.Left"] = left
        out["Stn.Right"] = right
        out["Stn.Beg"] = left

        inc = self.increment
        if inc is not None:
            out["Stn.Inc"] = inc
            out["Stn.GdpBeg"] = left
            out["Stn.GdpInc"] = inc
        return out

 
    def __str__(self) -> str:
        if self.values.size == 0:
            return "Station(empty)"
        span = self.span
        inc = self.increment
        span_txt = f"{span[0]:g}–{span[1]:g}" if span else "?"
        inc_txt = f"{inc:g}" if inc is not None else "ragged"
        return (
            f"Station(n={self.n_unique}, span={span_txt} {self.unit}, "
            f"inc={inc_txt})"
        )

    __repr__ = __str__

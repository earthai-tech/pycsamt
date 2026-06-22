# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Borehole — ground-truth data model for EM geological interpretation.

A :class:`Borehole` stores the depth-interval log collected during or
after drilling operations.  Intervals carry the true resistivity
(TRES) measured in the lab or from downhole logging, and a lithology
name.  Both are used by :class:`~pycsamt.interp.calibrate.ModelCalibrator`
to constrain 2-D inversion models.

Supported input formats
-----------------------
* **CSV** — single file with columns ``top, bottom, lithology[, resistivity]``
* **LAS 2.0** — industry well-log ASCII standard
* **dict / list** — programmatic construction

Example
-------
>>> bh = Borehole.from_csv("K1_borehole.csv", name="Bo", x=1050.0)
>>> bh.tres_at_depth(95.0)
180.0
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

__all__ = ["Interval", "Borehole"]

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Interval
# ---------------------------------------------------------------------------

@dataclass
class Interval:
    """A single depth interval in a borehole log.

    Parameters
    ----------
    top : float
        Depth to the top of the interval, metres.
    bottom : float
        Depth to the bottom of the interval, metres.
    lithology : str
        Geological formation / lithology name.
    resistivity : float or None
        True resistivity (TRES) in Ω·m (linear, **not** log₁₀).
        ``None`` when no electrical measurement is available.
    """

    top:          float
    bottom:       float
    lithology:    str
    resistivity:  Optional[float] = None

    # ------------------------------------------------------------------

    @property
    def thickness(self) -> float:
        return self.bottom - self.top

    def contains(self, z: float) -> bool:
        return self.top <= z < self.bottom

    def __post_init__(self):
        if self.bottom <= self.top:
            raise ValueError(
                f"Interval bottom ({self.bottom}) must be > top ({self.top})."
            )


# ---------------------------------------------------------------------------
# Borehole
# ---------------------------------------------------------------------------

class Borehole:
    """Borehole / well log with depth-interval data.

    Parameters
    ----------
    name : str
        Well / borehole identifier.
    x : float
        Position along the survey profile, metres.
    collar_elevation : float
        Surface elevation at the borehole collar, metres a.s.l.
        Used only when elevation-corrected exports are requested.
    intervals : list of Interval
        Depth-interval log, sorted ascending by *top*.
    """

    def __init__(
        self,
        name: str,
        x: float,
        intervals: Sequence[Interval],
        *,
        collar_elevation: float = 0.0,
    ) -> None:
        self.name = name
        self.x = float(x)
        self.collar_elevation = float(collar_elevation)
        self.intervals: List[Interval] = sorted(
            list(intervals), key=lambda iv: iv.top
        )

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def interval_at_depth(self, z: float) -> Optional[Interval]:
        """Return the interval that contains depth *z*, or ``None``."""
        for iv in self.intervals:
            if iv.contains(z):
                return iv
        return None

    def tres_at_depth(self, z: float) -> Optional[float]:
        """Return TRES (Ω·m) at depth *z*, or ``None`` if unknown."""
        iv = self.interval_at_depth(z)
        return iv.resistivity if iv is not None else None

    def lithology_at_depth(self, z: float) -> Optional[str]:
        """Return the lithology name at depth *z*, or ``None``."""
        iv = self.interval_at_depth(z)
        return iv.lithology if iv is not None else None

    def tres_column(self, z_centers: np.ndarray) -> np.ndarray:
        """Return TRES values at *z_centers* as a float array.

        Depths not covered by any interval are ``nan``.
        """
        z = np.asarray(z_centers, dtype=float)
        out = np.full(z.shape, np.nan)
        for i, zi in enumerate(z):
            v = self.tres_at_depth(float(zi))
            if v is not None:
                out[i] = v
        return out

    @property
    def max_depth(self) -> float:
        return self.intervals[-1].bottom if self.intervals else 0.0

    @property
    def min_depth(self) -> float:
        return self.intervals[0].top if self.intervals else 0.0

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    @classmethod
    def from_csv(
        cls,
        path: PathLike,
        *,
        name: Optional[str] = None,
        x: float = 0.0,
        collar_elevation: float = 0.0,
        delimiter: str = ",",
        top_col: str = "top",
        bottom_col: str = "bottom",
        lithology_col: str = "lithology",
        resistivity_col: str = "resistivity",
    ) -> "Borehole":
        """Load from a CSV file.

        Expected columns (case-insensitive header):
        ``top, bottom, lithology[, resistivity]``

        Parameters
        ----------
        path : path-like
        name : str, optional
            Defaults to the file stem.
        x : float
            Profile position.
        collar_elevation : float
        delimiter : str
        top_col, bottom_col, lithology_col, resistivity_col
            Column header names.

        Returns
        -------
        Borehole
        """
        p = Path(path)
        if name is None:
            name = p.stem

        intervals: List[Interval] = []
        with p.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter=delimiter)
            if reader.fieldnames is None:
                raise ValueError(f"CSV file has no header: {p}")
            headers_lower = {h.lower(): h for h in reader.fieldnames}

            def _col(key: str) -> str:
                return headers_lower.get(key.lower(), key)

            for row in reader:
                top = float(row[_col(top_col)])
                bot = float(row[_col(bottom_col)])
                lith = row.get(_col(lithology_col), "").strip()
                res_raw = row.get(_col(resistivity_col), "")
                res: Optional[float] = None
                if res_raw.strip() not in ("", "nan", "None", "NA"):
                    try:
                        res = float(res_raw)
                    except ValueError:
                        pass
                intervals.append(Interval(top=top, bottom=bot,
                                          lithology=lith, resistivity=res))

        return cls(name, x, intervals, collar_elevation=collar_elevation)

    @classmethod
    def from_las(
        cls,
        path: PathLike,
        *,
        x: float = 0.0,
        collar_elevation: float = 0.0,
        depth_curve: str = "DEPT",
        resistivity_curve: str = "RESD",
        lithology_curve: Optional[str] = "LITH",
        null_value: float = -9999.25,
        step: Optional[float] = None,
    ) -> "Borehole":
        """Load from a LAS 2.0 well-log file.

        Converts the continuous depth log into discrete intervals by
        grouping consecutive samples with the same lithology code.

        Parameters
        ----------
        path : path-like
        x : float
            Profile position.
        depth_curve, resistivity_curve, lithology_curve
            Curve mnemonics.  *lithology_curve* may be ``None`` to
            assign a generic label.
        null_value : float
            LAS null sentinel replaced with ``nan``.
        step : float, optional
            If provided, override the step value from the LAS header.

        Returns
        -------
        Borehole
        """
        p = Path(path)
        name = p.stem

        curves: Dict[str, List[float]] = {}
        well_name = name
        null_val = null_value

        section = None
        curve_order: List[str] = []
        data_lines: List[str] = []

        with p.open() as fh:
            for raw in fh:
                line = raw.rstrip()
                if line.startswith("~"):
                    section = line[1:].strip().upper()
                    continue
                if line.startswith("#"):
                    continue
                if section and section.startswith("W"):
                    if "WELL" in line.upper() and "." in line:
                        parts = line.split(".")
                        if len(parts) > 1:
                            val_comment = parts[1].split(":")
                            well_name = val_comment[0].strip() or well_name
                    if "NULL" in line.upper() and "." in line:
                        parts = line.split(".")
                        val = parts[1].split(":")[0].strip()
                        try:
                            null_val = float(val)
                        except ValueError:
                            pass
                elif section and section.startswith("C"):
                    if "." in line:
                        mnem = line.split(".")[0].strip().upper()
                        if mnem:
                            curve_order.append(mnem)
                elif section and section.startswith("A"):
                    if line.strip():
                        data_lines.append(line)

        for mnem in curve_order:
            curves[mnem] = []

        for dl in data_lines:
            vals = dl.split()
            for i, v in enumerate(vals):
                if i < len(curve_order):
                    try:
                        curves[curve_order[i]].append(float(v))
                    except ValueError:
                        curves[curve_order[i]].append(null_val)

        dept_key = depth_curve.upper()
        res_key = resistivity_curve.upper()
        lith_key = lithology_curve.upper() if lithology_curve else None

        if dept_key not in curves:
            raise ValueError(f"Depth curve {depth_curve!r} not found in {p}")

        depths = np.array(curves[dept_key])
        resis = np.array(curves.get(res_key, [null_val] * len(depths)))
        liths = np.array(curves.get(lith_key, [0] * len(depths))) if lith_key else None

        mask_null = np.isclose(resis, null_val, atol=1.0)
        resis = np.where(mask_null, np.nan, resis)

        # Build intervals: group consecutive same-lithology / same-depth-step
        intervals: List[Interval] = []
        if len(depths) < 2:
            return cls(well_name, x, intervals, collar_elevation=collar_elevation)

        ds = step if step is not None else float(np.median(np.diff(depths)))

        i = 0
        while i < len(depths):
            j = i + 1
            lith_i = int(liths[i]) if liths is not None else 0
            while j < len(depths):
                lith_j = int(liths[j]) if liths is not None else 0
                if lith_j != lith_i:
                    break
                j += 1
            top = float(depths[i])
            bottom = float(depths[j - 1]) + ds
            rho_chunk = resis[i:j]
            rho_mean = float(np.nanmean(rho_chunk)) if not np.all(np.isnan(rho_chunk)) else None
            intervals.append(Interval(
                top=top, bottom=bottom,
                lithology=str(lith_i),
                resistivity=rho_mean,
            ))
            i = j

        return cls(well_name, x, intervals, collar_elevation=collar_elevation)

    # ------------------------------------------------------------------
    # Serialisation helpers
    # ------------------------------------------------------------------

    def to_dataframe(self):
        """Return intervals as a :class:`pandas.DataFrame`."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError("pandas is required for Borehole.to_dataframe") from exc
        rows = [
            {
                "top": iv.top,
                "bottom": iv.bottom,
                "thickness": iv.thickness,
                "lithology": iv.lithology,
                "resistivity": iv.resistivity,
            }
            for iv in self.intervals
        ]
        return pd.DataFrame(rows)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "x": self.x,
            "collar_elevation": self.collar_elevation,
            "intervals": [
                {
                    "top": iv.top,
                    "bottom": iv.bottom,
                    "lithology": iv.lithology,
                    "resistivity": iv.resistivity,
                }
                for iv in self.intervals
            ],
        }

    def __repr__(self) -> str:
        return (
            f"Borehole({self.name!r}, x={self.x:.1f} m, "
            f"{len(self.intervals)} intervals, "
            f"depth={self.max_depth:.1f} m)"
        )

# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmData — read, write, and build the ModEM data file.

The ModEM data format is identical for 2-D and 3-D runs.  Each block
begins with a ``>`` header that identifies the component type, time
convention, units, rotation angle, origin offset, and data count.
Data rows follow the pattern::

    Period(s)  SiteCode  GG_Lat  GG_Lon  X(m)  Y(m)  Z(m)  Component  Real  Imag  Error

Supported component types
-------------------------
* ``TE_Impedance`` / ``TM_Impedance``  — 2-D MT profile data
* ``Full_Impedance``                   — 3-D full Z-tensor (ZXX/ZXY/ZYX/ZYY)
* ``Off_Diagonal_Impedance``           — 3-D off-diagonal (ZXY/ZYX)
* ``Determinant_Impedance``            — 3-D determinant
* ``Full_Vertical_Components``         — Tipper (HZX/HZY)
* ``Phase_Tensor``                     — Phase tensor components

Entry points
------------
``ModEmData.read(path)``
    Parse an existing data file.
``ModEmData.write(path)``
    Serialise to disk.
``ModEmData.from_edi(source, config)``
    Build from an EDI collection / site list.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .base   import ModEmBase
from .config import ModEmConfig

PathLike = Union[str, Path]

__all__ = ["ModEmData"]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_KNOWN_COMPONENTS = {
    "TE_Impedance":          ("TE",),
    "TM_Impedance":          ("TM",),
    "Full_Impedance":        ("ZXX", "ZXY", "ZYX", "ZYY"),
    "Off_Diagonal_Impedance":("ZXY", "ZYX"),
    "Determinant_Impedance": ("ZDet",),
    "Full_Vertical_Components": ("HZX", "HZY"),
    "Phase_Tensor":          ("PTxx", "PTxy", "PTyx", "PTyy"),
}

# Columns in the data record array
_COL_PERIOD    = 0
_COL_SITE_IDX  = 1   # 0-based index into site_names
_COL_X         = 2   # northing (m)
_COL_Y         = 3   # easting  (m)
_COL_Z         = 4   # elevation (m)
_COL_COMP_IDX  = 5   # index into component list for this block
_COL_REAL      = 6
_COL_IMAG      = 7
_COL_ERROR     = 8


# ---------------------------------------------------------------------------
# Low-level parser
# ---------------------------------------------------------------------------

def _parse_data(path: Path) -> dict:
    """Parse a ModEM data file.

    Returns
    -------
    dict with keys:
        comment, component_type, sign_convention, units, rotation_angle,
        origin, site_names, site_coords (dict name→(x,y,z)),
        periods (sorted descending), records (list of dicts per block)
    """
    with path.open("r", errors="replace") as fh:
        raw_lines = fh.readlines()

    comment_lines: list[str] = []
    blocks: list[dict] = []
    site_names:   list[str] = []
    site_coords:  dict = {}  # name → (x, y, z)
    site_idx_map: dict = {}  # name → 0-based int

    i = 0
    N = len(raw_lines)

    def _skip_blank():
        nonlocal i
        while i < N and not raw_lines[i].strip():
            i += 1

    while i < N:
        line = raw_lines[i].rstrip("\n")
        stripped = line.strip()
        i += 1

        if not stripped:
            continue

        # Comment header
        if stripped.startswith("#"):
            comment_lines.append(stripped[1:].strip())
            continue

        # Block header (> ComponentType)
        if stripped.startswith(">"):
            comp_type = stripped[1:].strip()
            # Read 5 more > lines: convention, units, rotation, origin(x y), nPeriods nSites
            block_meta: list[str] = []
            while i < N and len(block_meta) < 5:
                ln = raw_lines[i].strip()
                i += 1
                if ln.startswith(">"):
                    block_meta.append(ln[1:].strip())
                elif ln.startswith("#") or not ln:
                    continue
                else:
                    # reached data already (malformed?) — back up
                    i -= 1
                    break

            sign_conv = block_meta[0] if len(block_meta) > 0 else ""
            units_str = block_meta[1] if len(block_meta) > 1 else ""
            rot_angle = float(block_meta[2]) if len(block_meta) > 2 else 0.0
            origin_xy = [float(v) for v in block_meta[3].split()] if len(block_meta) > 3 else [0.0, 0.0]
            counts_ln = block_meta[4] if len(block_meta) > 4 else "0 0"
            counts    = [int(v) for v in counts_ln.split()]
            n_periods = counts[0] if counts else 0
            n_sites   = counts[1] if len(counts) > 1 else 0

            # Data rows
            data_rows: list[tuple] = []
            rows_read = 0
            expected  = n_periods * n_sites * len(_KNOWN_COMPONENTS.get(comp_type, ("Z",)))
            # read until we have all rows or hit next header
            while i < N and (expected == 0 or rows_read < expected):
                ln = raw_lines[i].strip()
                if not ln or ln.startswith("#") or ln.startswith(">"):
                    break
                i += 1
                parts = ln.split()
                if len(parts) < 11:
                    continue
                period  = float(parts[0])
                code    = parts[1]
                # lat, lon are parts[2,3] but we use X,Y in metres
                x_m     = float(parts[4])
                y_m     = float(parts[5])
                z_m     = float(parts[6])
                comp    = parts[7]
                real_v  = float(parts[8])
                imag_v  = float(parts[9])
                err_v   = float(parts[10])

                if code not in site_idx_map:
                    site_idx_map[code] = len(site_names)
                    site_names.append(code)
                    site_coords[code] = (x_m, y_m, z_m)

                data_rows.append((
                    period, site_idx_map[code],
                    x_m, y_m, z_m, comp,
                    real_v, imag_v, err_v,
                ))
                rows_read += 1

            blocks.append({
                "component_type": comp_type,
                "sign_convention": sign_conv,
                "units": units_str,
                "rotation_angle": rot_angle,
                "origin": origin_xy,
                "n_periods": n_periods,
                "n_sites":   n_sites,
                "rows": data_rows,
            })

    # Collect unique periods across all blocks, sorted descending
    all_periods: list[float] = []
    for blk in blocks:
        for row in blk["rows"]:
            all_periods.append(row[0])
    periods = np.unique(all_periods)[::-1]  # descending

    return {
        "comment":      "\n".join(comment_lines),
        "blocks":       blocks,
        "site_names":   site_names,
        "site_coords":  site_coords,
        "periods":      periods,
    }


# ---------------------------------------------------------------------------
# ModEmData
# ---------------------------------------------------------------------------

class ModEmData(ModEmBase):
    """ModEM data-file container.

    Attributes
    ----------
    comment : str
        Free-text provenance comment (written to ``#`` lines).
    blocks : list[dict]
        One entry per component-type block with keys
        ``component_type``, ``sign_convention``, ``units``,
        ``rotation_angle``, ``origin``, ``n_periods``, ``n_sites``,
        ``rows`` (list of 9-tuples: period, site_idx, x, y, z,
        component_str, real, imag, error).
    site_names : list[str]
        Station codes in order of first appearance.
    site_coords : dict[str, tuple]
        ``{code: (x_m, y_m, z_m)}`` — northing, easting, elevation.
    periods : np.ndarray
        Unique periods (s) across all blocks, sorted descending.
    """

    def __init__(
        self,
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.config:      ModEmConfig      = config or ModEmConfig()
        self.comment:     str              = "ModEM data file created by pycsamt"
        self.blocks:      List[dict]       = []
        self.site_names:  List[str]        = []
        self.site_coords: Dict[str, tuple] = {}
        self.periods:     np.ndarray       = np.array([])

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------
    @property
    def n_sites(self) -> int:
        return len(self.site_names)

    @property
    def n_periods(self) -> int:
        return len(self.periods)

    @property
    def offsets(self) -> np.ndarray:
        """Y-coordinate (easting, m) for each station in *site_names* order."""
        return np.array([self.site_coords.get(n, (0., 0., 0.))[1]
                         for n in self.site_names], dtype=float)

    @property
    def x_coords(self) -> np.ndarray:
        """X (northing, m) for each station in *site_names* order."""
        return np.array([self.site_coords.get(n, (0., 0., 0.))[0]
                         for n in self.site_names], dtype=float)

    @property
    def y_coords(self) -> np.ndarray:
        """Y (easting, m) for each station in *site_names* order."""
        return self.offsets

    @property
    def component_types(self) -> List[str]:
        return [b["component_type"] for b in self.blocks]

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "ModEmData":
        """Parse an existing ModEM data file.

        Parameters
        ----------
        path : path-like

        Returns
        -------
        ModEmData
        """
        p      = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"ModEM data file not found: {p}")
        parsed = _parse_data(p)
        obj    = cls(**kwargs)
        obj.comment     = parsed["comment"]
        obj.blocks      = parsed["blocks"]
        obj.site_names  = parsed["site_names"]
        obj.site_coords = parsed["site_coords"]
        obj.periods     = parsed["periods"]
        if obj.verbose:
            obj.logger.info(
                "ModEmData.read: %d sites, %d periods, %d blocks from %s",
                obj.n_sites, obj.n_periods, len(obj.blocks), p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write data to *path* in ModEM ASCII format.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        lines: list[str] = []
        lines.append(f"# {self.comment}\n")
        lines.append("# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) "
                     "Component Real Imag Error\n")

        for blk in self.blocks:
            comp  = blk["component_type"]
            sign  = blk.get("sign_convention", self.config.sign_convention)
            units = blk.get("units", self.config.units)
            rot   = blk.get("rotation_angle", 0.0)
            orig  = blk.get("origin", [0.0, 0.0])
            rows  = blk["rows"]

            # Count unique periods and sites in this block
            periods_in = sorted({r[0] for r in rows}, reverse=True)
            sites_in   = sorted({r[1] for r in rows})
            n_per      = len(periods_in)
            n_sit      = len(sites_in)

            lines.append(f"> {comp}\n")
            lines.append(f"> {sign}\n")
            lines.append(f"> {units}\n")
            lines.append(f"> {rot:.2f}\n")
            lines.append(f"> {orig[0]:.3f} {orig[1]:.3f}\n")
            lines.append(f"> {n_per} {n_sit}\n")

            for row in rows:
                period, sidx, x, y, z, comp_str, real, imag, err = row
                code = (
                    self.site_names[sidx]
                    if sidx < len(self.site_names) else f"S{sidx:04d}"
                )
                lines.append(
                    f"{period:>15.6E} {code:>8}  "
                    f"{0.000:>8.3f}  {0.000:>8.3f}  "
                    f"{x:>14.3f}  {y:>14.3f}  {z:>10.3f}  "
                    f"{comp_str:<6}"
                    f"  {real:>15.6E}  {imag:>15.6E}  {err:>15.6E}\n"
                )

        with p.open("w") as fh:
            fh.writelines(lines)
        return p

    # ------------------------------------------------------------------
    # Builder
    # ------------------------------------------------------------------
    @classmethod
    def from_edi(
        cls,
        source,
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ) -> "ModEmData":
        """Build a ModEM data file from EDI sites.

        Parameters
        ----------
        source : iterable of duck-typed site objects
            Each item must expose:

            * ``name`` (str)
            * ``coords`` (lat, lon, elev) or ``lat``/``lon``/``elev``
            * ``freq`` (array, Hz)
            * ``z`` (array, shape (n_freq, 2, 2), complex) — impedance tensor
              in the units matching *config.units*
            * ``z_err`` (array, same shape) — error estimate on Z
              (absolute, same units); or ``None`` for floor-only errors

            For 2-D data, ``z[:,0,1]`` = Z_TE, ``z[:,1,0]`` = Z_TM.

        config : ModEmConfig, optional

        Returns
        -------
        ModEmData
        """
        cfg  = config or ModEmConfig()
        obj  = cls(config=cfg, **kwargs)

        items = _normalise_source(source)
        if not items:
            raise ValueError("ModEmData.from_edi: source is empty")

        # ---- collect global unique periods (sorted descending) ----
        all_freqs: list[float] = []
        for it in items:
            freqs = np.asarray(getattr(it, "freq", []), dtype=float)
            all_freqs.extend(freqs[freqs > 0].tolist())
        if not all_freqs:
            raise ValueError("ModEmData.from_edi: no valid frequencies found")
        global_periods = _unique_periods(all_freqs, rtol=0.01)

        # ---- site coordinates ----
        names: list[str] = []
        coords: dict = {}
        for it in items:
            name = str(getattr(it, "name", None) or f"S{len(names):04d}")
            names.append(name)
            c = getattr(it, "coords", None)
            if c is not None and len(c) >= 2:
                lat, lon = float(c[0]), float(c[1])
                elev = float(c[2]) if len(c) > 2 else 0.0
            else:
                lat  = float(getattr(it, "lat",  0.0))
                lon  = float(getattr(it, "lon",  0.0))
                elev = float(getattr(it, "elev", 0.0))
            # Convert lat/lon to X(northing)/Y(easting) in metres
            x_m, y_m = _latlon_to_xy(lat, lon, lat_ref=lat, lon_ref=lon)
            coords[name] = (x_m, y_m, elev)

        # Recompute offsets relative to centroid
        if len(names) > 1:
            lats = np.array([
                float(getattr(it, "coords", [0,0,0])[0]
                      if getattr(it, "coords", None) is not None
                      else getattr(it, "lat", 0.0))
                for it in items
            ])
            lons = np.array([
                float(getattr(it, "coords", [0,0,0])[1]
                      if getattr(it, "coords", None) is not None
                      else getattr(it, "lon", 0.0))
                for it in items
            ])
            lat0, lon0 = float(lats.mean()), float(lons.mean())
            for i, (it, name) in enumerate(zip(items, names)):
                c = getattr(it, "coords", None)
                lat = float(c[0]) if c is not None else float(getattr(it, "lat", 0.0))
                lon = float(c[1]) if c is not None else float(getattr(it, "lon", 0.0))
                elev = float(c[2]) if (c is not None and len(c) > 2) else 0.0
                x_m, y_m = _latlon_to_xy(lat, lon, lat_ref=lat0, lon_ref=lon0)
                coords[name] = (x_m, y_m, elev)

        obj.site_names  = names
        obj.site_coords = coords
        obj.periods     = global_periods

        # ---- determine which component blocks to write ----
        comp_type = cfg.component_type
        if comp_type not in _KNOWN_COMPONENTS:
            raise ValueError(
                f"Unknown component_type: {comp_type!r}. "
                f"Choose from: {list(_KNOWN_COMPONENTS)}"
            )

        # ---- build data rows ----
        rows: list[tuple] = []
        for si, (it, name) in enumerate(zip(items, names)):
            freqs_site = np.asarray(getattr(it, "freq", []), dtype=float)
            z_arr      = getattr(it, "z",     None)
            z_err_arr  = getattr(it, "z_err", None)
            x_m, y_m, z_elev = coords[name]

            if z_arr is None:
                continue
            z_arr = np.asarray(z_arr, dtype=complex)

            for pi, period in enumerate(global_periods):
                freq_target = 1.0 / period
                fi = _match_freq(freqs_site, freq_target, rtol=0.01)
                if fi is None:
                    continue

                z_max_at_freq = float(np.max(np.abs(z_arr[fi, :, :])))
                for comp_str, (ri, ci) in _comp_indices(comp_type):
                    z_val  = complex(z_arr[fi, ri, ci])
                    if z_err_arr is not None:
                        z_e = abs(float(np.asarray(z_err_arr)[fi, ri, ci]))
                    else:
                        z_e = 0.0
                    # Floor relative to max |Z| at this frequency so zero-valued
                    # components (e.g. ZXX/ZYY for 1-D structure) still get a
                    # non-zero error estimate.
                    floor = max(z_max_at_freq * cfg.error_floor_z, cfg.error_floor_z_floor)
                    z_e   = max(z_e, floor)

                    rows.append((
                        period, si,
                        x_m, y_m, z_elev,
                        comp_str,
                        float(z_val.real), float(z_val.imag), float(z_e),
                    ))

        n_per_in_rows = len({r[0] for r in rows})
        n_sit_in_rows = len({r[1] for r in rows})

        obj.blocks = [{
            "component_type":  comp_type,
            "sign_convention": cfg.sign_convention,
            "units":           cfg.units,
            "rotation_angle":  0.0,
            "origin":          [0.0, 0.0],
            "n_periods":       n_per_in_rows,
            "n_sites":         n_sit_in_rows,
            "rows":            rows,
        }]

        if obj.verbose:
            obj.logger.info(
                "ModEmData.from_edi: %d sites, %d periods, component=%s",
                len(names), len(global_periods), comp_type,
            )
        return obj


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _normalise_source(source) -> list:
    if hasattr(source, "_items"):
        return list(source)
    if hasattr(source, "edic"):
        try:
            from pycsamt.site.base import Site
            return [Site(e) for e in source.edic]
        except Exception:
            return list(source.edic)
    return list(source)


def _unique_periods(all_freqs: list, rtol: float = 0.01) -> np.ndarray:
    """Return unique periods (s) sorted ascending, merging near-duplicate freqs."""
    arr = np.sort(np.unique(np.asarray(all_freqs, dtype=float)))[::-1]  # desc freqs
    if arr.size == 0:
        return np.array([], dtype=float)
    merged = [float(arr[0])]
    for f in arr[1:]:
        if abs(f - merged[-1]) / max(abs(f), abs(merged[-1]), 1e-30) > rtol:
            merged.append(float(f))
    # convert freqs → periods; ascending period = ascending 1/f means descending f,
    # so reverse so that short periods come first
    return 1.0 / np.array(merged[::-1], dtype=float)


def _match_freq(freqs: np.ndarray, target: float, rtol: float = 0.01) -> Optional[int]:
    """Return index of *target* in *freqs* within *rtol* tolerance, or None."""
    if freqs.size == 0:
        return None
    diff = np.abs(freqs - target) / max(abs(target), 1e-30)
    idx  = int(np.argmin(diff))
    return idx if diff[idx] <= rtol else None


def _comp_indices(comp_type: str) -> List[Tuple[str, Tuple[int, int]]]:
    """Return list of (component_str, (row_idx, col_idx)) for Z array."""
    _map = {
        "ZXX": (0, 0), "ZXY": (0, 1), "ZYX": (1, 0), "ZYY": (1, 1),
        "TE":  (0, 1), "TM":  (1, 0),
        "HZX": (0, 0), "HZY": (0, 1),
        "ZDet": None,
    }
    comps = _KNOWN_COMPONENTS.get(comp_type, ())
    result = []
    for c in comps:
        idx = _map.get(c)
        if idx is not None:
            result.append((c, idx))
    return result


def _latlon_to_xy(
    lat: float, lon: float,
    lat_ref: float, lon_ref: float,
) -> Tuple[float, float]:
    """Approximate lat/lon → (x_northing_m, y_easting_m) via flat-Earth."""
    R    = 6_371_000.0
    dlat = np.radians(lat - lat_ref)
    dlon = np.radians(lon - lon_ref)
    x_m  = R * dlat
    y_m  = R * dlon * np.cos(np.radians(lat_ref))
    return float(x_m), float(y_m)

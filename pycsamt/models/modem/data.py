# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Read, write, and build ModEM data files.

The module defines :class:`ModEmData`, the container for
observed or predicted ModEM response data. It supports parsing
existing ModEM ASCII files, writing populated objects, and
building data rows from EDI-like station objects.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .base import ModEmBase
from .config import ModEmConfig
from .doc import _modem_param_docs as _params

PathLike = str | Path

__all__ = ["ModEmData"]

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

_KNOWN_COMPONENTS = {
    "TE_Impedance": ("TE",),
    "TM_Impedance": ("TM",),
    "Full_Impedance": ("ZXX", "ZXY", "ZYX", "ZYY"),
    "Off_Diagonal_Impedance": ("ZXY", "ZYX"),
    "Determinant_Impedance": ("ZDet",),
    "Full_Vertical_Components": ("HZX", "HZY"),
    "Phase_Tensor": ("PTxx", "PTxy", "PTyx", "PTyy"),
}

# Columns in the data record array
_COL_PERIOD = 0
_COL_SITE_IDX = 1   # 0-based index into site_names
_COL_X = 2   # northing (m)
_COL_Y = 3   # easting  (m)
_COL_Z = 4   # elevation (m)
_COL_COMP_IDX = 5   # index into component list for this block
_COL_REAL = 6
_COL_IMAG = 7
_COL_ERROR = 8


# ----------------------------------------------------------------------
# Low-level parser
# ----------------------------------------------------------------------

def _parse_data(path: Path) -> dict:
    """Parse a ModEM data file into plain Python values."""
    with path.open("r", errors="replace") as fh:
        raw_lines = fh.readlines()

    comment_lines: list[str] = []
    blocks: list[dict] = []
    site_names: list[str] = []
    site_coords: dict = {}  # name -> (x, y, z)
    site_idx_map: dict = {}  # name -> 0-based int

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
            # Read convention, units, rotation, origin, counts.
            block_meta: list[str] = []
            while i < N and len(block_meta) < 5:
                ln = raw_lines[i].strip()
                i += 1
                if ln.startswith(">"):
                    block_meta.append(ln[1:].strip())
                elif ln.startswith("#") or not ln:
                    continue
                else:
                    # reached data already (malformed?) - back up
                    i -= 1
                    break

            sign_conv = block_meta[0] if len(block_meta) > 0 else ""
            units_str = block_meta[1] if len(block_meta) > 1 else ""
            rot_angle = (
                float(block_meta[2]) if len(block_meta) > 2 else 0.0
            )
            if len(block_meta) > 3:
                origin_xy = [float(v) for v in block_meta[3].split()]
            else:
                origin_xy = [0.0, 0.0]
            counts_ln = block_meta[4] if len(block_meta) > 4 else "0 0"
            counts = [int(v) for v in counts_ln.split()]
            n_periods = counts[0] if counts else 0
            n_sites = counts[1] if len(counts) > 1 else 0

            # Data rows
            data_rows: list[tuple] = []
            rows_read = 0
            expected = n_periods * n_sites * len(
                _KNOWN_COMPONENTS.get(comp_type, ("Z",))
            )
            # read until we have all rows or hit next header
            while i < N and (expected == 0 or rows_read < expected):
                ln = raw_lines[i].strip()
                if not ln or ln.startswith("#") or ln.startswith(">"):
                    break
                i += 1
                parts = ln.split()
                if len(parts) < 11:
                    continue
                period = float(parts[0])
                code = parts[1]
                # lat, lon are parts[2,3] but we use X,Y in metres
                x_m = float(parts[4])
                y_m = float(parts[5])
                z_m = float(parts[6])
                comp = parts[7]
                real_v = float(parts[8])
                imag_v = float(parts[9])
                err_v = float(parts[10])

                if code not in site_idx_map:
                    site_idx_map[code] = len(site_names)
                    site_names.append(code)
                    site_coords[code] = (x_m, y_m, z_m)

                data_rows.append(
                    (
                        period,
                        site_idx_map[code],
                        x_m,
                        y_m,
                        z_m,
                        comp,
                        real_v,
                        imag_v,
                        err_v,
                    )
                )
                rows_read += 1

            blocks.append(
                {
                    "component_type": comp_type,
                    "sign_convention": sign_conv,
                    "units": units_str,
                    "rotation_angle": rot_angle,
                    "origin": origin_xy,
                    "n_periods": n_periods,
                    "n_sites": n_sites,
                    "rows": data_rows,
                }
            )

    # Collect unique periods across all blocks, sorted descending
    all_periods: list[float] = []
    for blk in blocks:
        for row in blk["rows"]:
            all_periods.append(row[0])
    periods = np.unique(all_periods)[::-1]  # descending

    return {
        "comment": "\n".join(comment_lines),
        "blocks": blocks,
        "site_names": site_names,
        "site_coords": site_coords,
        "periods": periods,
    }


# ----------------------------------------------------------------------
# ModEmData
# ----------------------------------------------------------------------

class ModEmData(ModEmBase):
    def __init__(
        self,
        config: ModEmConfig | None = None,
        **kwargs,
    ):
        """Initialize an empty ModEM data container."""
        super().__init__(**kwargs)
        self.config: ModEmConfig = config or ModEmConfig()
        self.comment: str = "ModEM data file created by pycsamt"
        self.blocks: list[dict] = []
        self.site_names: list[str] = []
        self.site_coords: dict[str, tuple] = {}
        self.periods: np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------
    @property
    def n_sites(self) -> int:
        """Number of stations represented by the data object."""
        return len(self.site_names)

    @property
    def n_periods(self) -> int:
        """Number of unique periods represented by the object."""
        return len(self.periods)

    @property
    def offsets(self) -> np.ndarray:
        """Return station easting offsets in metres."""
        return np.array(
            [
                self.site_coords.get(n, (0.0, 0.0, 0.0))[1]
                for n in self.site_names
            ],
            dtype=float,
        )

    @property
    def x_coords(self) -> np.ndarray:
        """Return station northing coordinates in metres."""
        return np.array(
            [
                self.site_coords.get(n, (0.0, 0.0, 0.0))[0]
                for n in self.site_names
            ],
            dtype=float,
        )

    @property
    def y_coords(self) -> np.ndarray:
        """Return station easting coordinates in metres."""
        return self.offsets

    @property
    def component_types(self) -> list[str]:
        """Return component-type names present in the data blocks."""
        return [b["component_type"] for b in self.blocks]

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> ModEmData:
        """Parse an existing ModEM data file.

        Parameters
        ----------
        path : path-like
            Path to a ModEM ASCII data file. The file may
            contain one or more ``>`` component blocks and
            optional ``#`` comment lines.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmData`, commonly ``config``, ``verbose``,
            or ``logger``.

        Returns
        -------
        ModEmData
            Parsed data object with comment, blocks, station
            names, station coordinates, and unique periods
            populated.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.

        Examples
        --------
        >>> from pycsamt.models.modem.data import ModEmData
        >>> data = ModEmData.read("data.dat")
        >>> data.n_sites > 0
        True
        """
        p = Path(path)
        if not p.exists():
            msg = f"ModEM data file not found: {p}"
            raise FileNotFoundError(msg)
        parsed = _parse_data(p)
        obj = cls(**kwargs)
        obj.comment = parsed["comment"]
        obj.blocks = parsed["blocks"]
        obj.site_names = parsed["site_names"]
        obj.site_coords = parsed["site_coords"]
        obj.periods = parsed["periods"]
        if obj.verbose:
            obj.logger.info(
                "ModEmData.read: %d sites, %d periods, "
                "%d blocks from %s",
                obj.n_sites,
                obj.n_periods,
                len(obj.blocks),
                p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write data to ``path`` in ModEM ASCII format.

        Parameters
        ----------
        path : path-like
            Destination data file. Parent directories are
            created before writing. Existing files are
            overwritten.

        Returns
        -------
        pathlib.Path
            Path passed to the writer, converted to
            :class:`pathlib.Path`.

        Notes
        -----
        Latitude and longitude columns are currently written as
        zeros because the object stores and uses local metre
        coordinates for ModEM modelling. The local ``X(m)``,
        ``Y(m)``, and ``Z(m)`` columns carry the station
        positions used by model builders and runners.

        Examples
        --------
        >>> from pycsamt.models.modem.data import ModEmData
        >>> data = ModEmData.read("data.dat")
        >>> path = data.write("copy.dat")
        >>> path.name
        'copy.dat'
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        lines: list[str] = []
        lines.append(f"# {self.comment}\n")
        lines.append("# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) "
                     "Component Real Imag Error\n")

        for blk in self.blocks:
            comp = blk["component_type"]
            sign = blk.get(
                "sign_convention",
                self.config.sign_convention,
            )
            units = blk.get("units", self.config.units)
            rot = blk.get("rotation_angle", 0.0)
            orig = blk.get("origin", [0.0, 0.0])
            rows = blk["rows"]

            # Count unique periods and sites in this block
            periods_in = sorted({r[0] for r in rows}, reverse=True)
            sites_in = sorted({r[1] for r in rows})
            n_per = len(periods_in)
            n_sit = len(sites_in)

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
        config: ModEmConfig | None = None,
        **kwargs,
    ) -> ModEmData:
        """Build a ModEM data file from EDI sites.

        The builder accepts a collection of site-like objects
        and writes one ModEM component block according to
        ``config.component_type``. Frequencies are merged
        across stations with a relative tolerance, converted
        to periods, and stored in descending period order.
        Impedance errors are floored using both relative and
        absolute thresholds from ``config``.

        Parameters
        ----------
        source : iterable of duck-typed site objects
            Each item must expose:

            * ``name`` (str)
            * ``coords`` (lat, lon, elev) or ``lat``/``lon``/``elev``
            * ``freq`` (array, Hz)
            * ``z`` (array, shape (n_freq, 2, 2), complex):
              impedance tensor in units matching ``config.units``.
            * ``z_err`` (array, same shape): error estimate on Z
              (absolute, same units); or ``None`` for floor-only errors

            For 2-D data, ``z[:,0,1]`` = Z_TE, ``z[:,1,0]`` = Z_TM.

        config : ModEmConfig, optional
            Configuration controlling component selection, sign
            convention, units, frequency limits, and impedance
            error floors. If omitted, defaults are used.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmData`.

        Returns
        -------
        ModEmData
            Populated data object ready to be written with
            :meth:`write` or passed to model builders.

        Raises
        ------
        ValueError
            If ``source`` is empty, contains no valid positive
            frequencies, or selects an unknown component type.

        Examples
        --------
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.data import ModEmData
        >>> cfg = ModEmConfig(component_type="Full_Impedance")
        >>> data = ModEmData.from_edi(sites, config=cfg)
        >>> data.component_types
        ['Full_Impedance']
        """
        cfg = config or ModEmConfig()
        obj = cls(config=cfg, **kwargs)

        items = _normalise_source(source)
        if not items:
            raise ValueError("ModEmData.from_edi: source is empty")

        # ---- collect global unique periods (sorted descending) ----
        all_freqs: list[float] = []
        for it in items:
            freqs = np.asarray(getattr(it, "freq", []), dtype=float)
            all_freqs.extend(freqs[freqs > 0].tolist())
        if not all_freqs:
            msg = "ModEmData.from_edi: no valid frequencies found"
            raise ValueError(msg)
        global_periods = _unique_periods(all_freqs, rtol=0.01)

        # ---- site coordinates ----
        names: list[str] = []
        coords: dict = {}
        for it in items:
            fallback = f"S{len(names):04d}"
            name = str(getattr(it, "name", None) or fallback)
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
            for it, name in zip(items, names):
                c = getattr(it, "coords", None)
                if c is not None:
                    lat = float(c[0])
                    lon = float(c[1])
                    elev = float(c[2]) if len(c) > 2 else 0.0
                else:
                    lat = float(getattr(it, "lat", 0.0))
                    lon = float(getattr(it, "lon", 0.0))
                    elev = float(getattr(it, "elev", 0.0))
                x_m, y_m = _latlon_to_xy(
                    lat,
                    lon,
                    lat_ref=lat0,
                    lon_ref=lon0,
                )
                coords[name] = (x_m, y_m, elev)

        obj.site_names = names
        obj.site_coords = coords
        obj.periods = global_periods

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
            freqs_site = np.asarray(
                getattr(it, "freq", []),
                dtype=float,
            )
            z_arr = getattr(it, "z", None)
            z_err_arr = getattr(it, "z_err", None)
            x_m, y_m, z_elev = coords[name]

            if z_arr is None:
                continue
            z_arr = np.asarray(z_arr, dtype=complex)

            for period in global_periods:
                freq_target = 1.0 / period
                fi = _match_freq(freqs_site, freq_target, rtol=0.01)
                if fi is None:
                    continue

                z_max_at_freq = float(np.max(np.abs(z_arr[fi, :, :])))
                for comp_str, (ri, ci) in _comp_indices(comp_type):
                    z_val = complex(z_arr[fi, ri, ci])
                    if z_err_arr is not None:
                        z_e = abs(
                            float(np.asarray(z_err_arr)[fi, ri, ci])
                        )
                    else:
                        z_e = 0.0
                    # Floor relative to max |Z| at this frequency.
                    floor = max(
                        z_max_at_freq * cfg.error_floor_z,
                        cfg.error_floor_z_floor,
                    )
                    z_e = max(z_e, floor)

                    rows.append(
                        (
                            period,
                            si,
                            x_m,
                            y_m,
                            z_elev,
                            comp_str,
                            float(z_val.real),
                            float(z_val.imag),
                            float(z_e),
                        )
                    )

        n_per_in_rows = len({r[0] for r in rows})
        n_sit_in_rows = len({r[1] for r in rows})

        obj.blocks = [
            {
                "component_type": comp_type,
                "sign_convention": cfg.sign_convention,
                "units": cfg.units,
                "rotation_angle": 0.0,
                "origin": [0.0, 0.0],
                "n_periods": n_per_in_rows,
                "n_sites": n_sit_in_rows,
                "rows": rows,
            }
        ]

        if obj.verbose:
            obj.logger.info(
                "ModEmData.from_edi: %d sites, %d periods, "
                "component=%s",
                len(names), len(global_periods), comp_type,
            )
        return obj


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _normalise_source(source) -> list:
    """Return a list of site-like objects from supported sources."""
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
    """Return unique periods sorted descending."""
    arr = np.sort(np.unique(np.asarray(all_freqs, dtype=float)))[::-1]
    if arr.size == 0:
        return np.array([], dtype=float)
    merged = [float(arr[0])]
    for f in arr[1:]:
        denom = max(abs(f), abs(merged[-1]), 1e-30)
        if abs(f - merged[-1]) / denom > rtol:
            merged.append(float(f))
    # Convert descending frequencies to descending periods.
    return 1.0 / np.array(merged[::-1], dtype=float)


def _match_freq(
    freqs: np.ndarray,
    target: float,
    rtol: float = 0.01,
) -> int | None:
    """Return the matching frequency index within tolerance."""
    if freqs.size == 0:
        return None
    diff = np.abs(freqs - target) / max(abs(target), 1e-30)
    idx = int(np.argmin(diff))
    return idx if diff[idx] <= rtol else None


def _comp_indices(comp_type: str) -> list[tuple[str, tuple[int, int]]]:
    """Return component names and tensor indices for a Z array."""
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
    lat: float,
    lon: float,
    lat_ref: float,
    lon_ref: float,
) -> tuple[float, float]:
    """Approximate lat/lon as local northing/easting metres."""
    R = 6_371_000.0
    dlat = np.radians(lat - lat_ref)
    dlon = np.radians(lon - lon_ref)
    x_m = R * dlat
    y_m = R * dlon * np.cos(np.radians(lat_ref))
    return float(x_m), float(y_m)


ModEmData.__doc__ = rf"""
Represent observed or predicted ModEM response data.

``ModEmData`` stores the ModEM ASCII data-file structure used
by both 2-D and 3-D workflows. A file is composed of one or
more component blocks. Each block starts with ``>`` metadata
lines describing the component family, time convention, units,
rotation, origin, and counts. Data rows then store period,
station code, local coordinates, component name, complex value,
and uncertainty.

The complex impedance tensor is commonly written as

.. math::

   \mathbf{{Z}} =
   \begin{{bmatrix}}
   Z_{{xx}} & Z_{{xy}} \\
   Z_{{yx}} & Z_{{yy}}
   \end{{bmatrix}},

with each selected component stored as real and imaginary
columns. The data builder applies an error floor

.. math::

   \sigma_Z = \max(\sigma_{{src}},
                   \epsilon |Z|_{{max}},
                   \sigma_{{min}}),

where :math:`\epsilon` is ``config.error_floor_z`` and
:math:`\sigma_{min}` is ``config.error_floor_z_floor``.

Parameters
----------
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
config : ModEmConfig
    Configuration object used for sign convention, units,
    component selection, frequency limits, and error floors.
{_params.data.comment}
{_params.data.blocks}
{_params.data.site_names}
site_coords : dict[str, tuple]
    Mapping from station name to ``(x_m, y_m, z_m)`` local
    coordinates. The x coordinate is local northing, the y
    coordinate is local easting, and z is elevation in metres.
{_params.data.periods}

Notes
-----
Supported component families include:

* ``"TE_Impedance"`` and ``"TM_Impedance"`` for 2-D profile
  workflows;
* ``"Full_Impedance"`` for ``ZXX``, ``ZXY``, ``ZYX``, and
  ``ZYY``;
* ``"Off_Diagonal_Impedance"`` for ``ZXY`` and ``ZYX``;
* ``"Determinant_Impedance"`` for determinant-style data;
* ``"Full_Vertical_Components"`` for tipper-like components;
* ``"Phase_Tensor"`` for phase-tensor component names.

``from_edi`` currently builds impedance component rows from
``z`` and ``z_err`` arrays. Component families whose values
require derived calculations, such as determinant or phase
tensor responses, may need additional preprocessing before
they can be represented fully.

See Also
--------
ModEmConfig
    Supplies component, unit, sign, frequency, and error
    settings.
InputBuilder
    Builds ModEM data and the matching model/control files.
ModEmModel2D
    Uses station coordinates from data to build 2-D models.
ModEmModel3D
    Uses station coordinates from data to build 3-D models.
ModEmRunner
    Passes observed data files to the ModEM executable.

Examples
--------
Read an existing data file:

>>> from pycsamt.models.modem.data import ModEmData
>>> data = ModEmData.read("data.dat")
>>> data.n_sites > 0
True

Build data from EDI-like station objects:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> cfg = ModEmConfig(component_type="Off_Diagonal_Impedance")
>>> data = ModEmData.from_edi(sites, config=cfg)
>>> data.component_types
['Off_Diagonal_Impedance']

Write a ModEM data file:

>>> path = data.write("modem_run/data.dat")
>>> path.name
'data.dat'

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D., and Tandon,
   K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""

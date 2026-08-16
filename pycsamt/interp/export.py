# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""export — write EM interpretation results to industry formats.

Supported formats
-----------------
* **Oasis Montaj XYZ** — point-data profile readable by Geosoft Oasis Montaj
  and any XYZ-aware GIS tool.
* **LAS 2.0** — CWLS well-log ASCII standard; compatible with petrophyscial
  software (Petrel, Kingdom, IHS Markit).
* **CSV** — flat table; station, depth, log10-rho, lithology.
* **VTK** — ASCII rectilinear-grid format readable by Paraview, QGIS, and
  other 3-D viewers.
* **Golden Software Surfer** — ``to_surfer_grid`` (DSAA regular grid,
  opens directly as a finished map/section) and ``to_surfer_xyz``
  (scattered ``X Y Z`` points at the model's real, possibly non-uniform
  cell centres, gridded by Surfer itself). Both accept *any* 2-D
  inversion result — Occam2D, ModEM, the backend-neutral
  :class:`~pycsamt.inversion.results.InversionResult`, an AI agent
  result, or a raw array — via :meth:`~pycsamt.interp.ResistivityModel.from_any`.

Example
-------
>>> from pycsamt.interp import export
>>> export.to_oasis_montaj_xyz(logs, "profile_K1.xyz")  # doctest: +SKIP
>>> export.to_las(logs[0], "S17.las")  # doctest: +SKIP
>>> export.to_csv(logs, "all_stations.csv")  # doctest: +SKIP
>>> export.to_surfer_grid(model, "section.grd")  # doctest: +SKIP
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Any, Union

import numpy as np

from ._base import ResistivityModel
from ..geology.lithology import StratigraphicLog

__all__ = [
    "to_oasis_montaj_xyz",
    "to_las",
    "to_csv",
    "to_vtk",
    "to_surfer_grid",
    "to_surfer_xyz",
]

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Oasis Montaj XYZ
# ---------------------------------------------------------------------------


def to_oasis_montaj_xyz(
    logs: Sequence[StratigraphicLog],
    path: PathLike,
    *,
    y: float = 0.0,
    elevation: np.ndarray | None = None,
    log_rho: bool = True,
    channels: list[str] | None = None,
) -> Path:
    """Write pseudo-stratigraphic logs to Oasis Montaj XYZ format.

    Each station is written as a separate ``/ Line`` block.  The depth
    axis is encoded as a negative Z value (depth below surface).

    Parameters
    ----------
    logs : list of StratigraphicLog
    path : path-like
        Output file path (extension ``.xyz`` recommended).
    y : float
        Northing to assign to all points (single-profile surveys).
    elevation : ndarray (n_stations,), optional
        Surface elevation for each station in metres a.s.l.
        If provided, Z = elevation − depth.  Otherwise Z = −depth.
    log_rho : bool
        Write :math:`\\log_{10}(\\rho)` (default) or linear ρ (Ω·m).
    channels : list of str, optional
        Override column headers.  Defaults to
        ``['X', 'Y', 'Z', 'RESD', 'LITH']``.

    Returns
    -------
    Path
        Written file path.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    hdr = channels if channels else ["X", "Y", "Z", "RESD", "LITH"]

    with out.open("w") as fh:
        fh.write("/ pycsamt.interp — Oasis Montaj XYZ\n")
        fh.write(f"/ Generated: {datetime.utcnow().isoformat()}\n")
        fh.write("/ " + "  ".join(hdr) + "\n")

        for k, log in enumerate(logs):
            elev = float(elevation[k]) if elevation is not None else 0.0
            fh.write(f"/ Line {log.station_name}\n")
            for _iz, (z, rho) in enumerate(zip(log.z_centers, log.rho_log10)):
                if np.isnan(rho):
                    continue
                x_val = log.station_x
                z_val = elev - float(z)
                rho_val = float(rho) if log_rho else float(10.0**rho)
                # find lithology for this depth cell
                lith = ""
                for layer in log.layers:
                    if layer.top <= float(z) < layer.bottom:
                        lith = layer.lithology
                        break
                lith_safe = lith.replace(" ", "_")
                fh.write(
                    f"{x_val:>14.3f}  {y:>14.3f}  {z_val:>12.3f}  "
                    f"{rho_val:>10.5f}  {lith_safe}\n"
                )
            fh.write("\n")

    return out


# ---------------------------------------------------------------------------
# LAS 2.0
# ---------------------------------------------------------------------------


def to_las(
    log: StratigraphicLog,
    path: PathLike,
    *,
    well_name: str = "",
    company: str = "pycsamt",
    null_value: float = -9999.25,
    log_rho: bool = True,
) -> Path:
    """Write a single station log to LAS 2.0 format.

    Parameters
    ----------
    log : StratigraphicLog
    path : path-like
        Output ``.las`` file.
    well_name : str
        Well / station identifier; defaults to ``log.station_name``.
    company : str
    null_value : float
        LAS null sentinel for missing values.
    log_rho : bool
        Write :math:`\\log_{10}(\\rho)` (default) or linear ρ.

    Returns
    -------
    Path
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    w_name = well_name or log.station_name

    z = log.z_centers
    rho = log.rho_log10.copy()
    if not log_rho:
        rho = 10.0**rho

    # build lithology-code column
    lith_codes = np.full(len(z), 0, dtype=int)
    for i, zi in enumerate(z):
        for layer in log.layers:
            if layer.top <= float(zi) < layer.bottom:
                lith_codes[i] = hash(layer.lithology) % 1000
                break

    step = float(np.median(np.diff(z))) if len(z) > 1 else 1.0
    rho_unit = "LOG10.OHMM" if log_rho else "OHMM"

    lines = [
        "~VERSION INFORMATION",
        " VERS.                  2.0:   CWLS log ASCII Standard - VERSION 2.0",
        " WRAP.                  NO:    ONE LINE PER DEPTH STEP",
        "~WELL INFORMATION",
        f" STRT.M          {z[0]:.4f}:  START DEPTH",
        f" STOP.M          {z[-1]:.4f}:  STOP DEPTH",
        f" STEP.M          {step:.4f}:  STEP",
        f" NULL.           {null_value}:  NULL VALUE",
        f" COMP.           {company}:  COMPANY",
        f" WELL.           {w_name}:  WELL",
        f" DATE.           {datetime.utcnow().strftime('%Y-%m-%d')}:  DATE",
        "~CURVE INFORMATION",
        " DEPT.M                 :  DEPTH",
        f" RESD.{rho_unit}  :  RESISTIVITY (EM-derived)",
        " LITH.                  :  LITHOLOGY CODE (hash % 1000)",
        "~A  DEPT         RESD         LITH",
    ]

    with out.open("w") as fh:
        fh.write("\n".join(lines) + "\n")
        for zi, ri, li in zip(z, rho, lith_codes):
            ri_out = ri if not np.isnan(ri) else null_value
            fh.write(f"  {zi:10.4f}  {ri_out:12.5f}  {li:6d}\n")

    return out


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------


def to_csv(
    logs: Sequence[StratigraphicLog],
    path: PathLike,
    *,
    log_rho: bool = True,
) -> Path:
    """Write all station logs to a flat CSV file.

    Columns: ``station, x_m, depth_m, rho_log10, rho_ohm_m, lithology``

    Parameters
    ----------
    logs : list of StratigraphicLog
    path : path-like
    log_rho : bool
        Include the ``rho_log10`` column (always written); if ``False``
        only ``rho_ohm_m`` is added.

    Returns
    -------
    Path
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "station",
        "x_m",
        "depth_m",
        "rho_log10",
        "rho_ohm_m",
        "lithology",
    ]

    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for log in logs:
            for zi, rho_l in zip(log.z_centers, log.rho_log10):
                if np.isnan(rho_l):
                    continue
                lith = ""
                for layer in log.layers:
                    if layer.top <= float(zi) < layer.bottom:
                        lith = layer.lithology
                        break
                writer.writerow(
                    {
                        "station": log.station_name,
                        "x_m": f"{log.station_x:.3f}",
                        "depth_m": f"{float(zi):.3f}",
                        "rho_log10": f"{float(rho_l):.5f}",
                        "rho_ohm_m": f"{10.0 ** float(rho_l):.3f}",
                        "lithology": lith,
                    }
                )

    return out


# ---------------------------------------------------------------------------
# VTK rectilinear grid
# ---------------------------------------------------------------------------


def to_vtk(
    model: ResistivityModel,
    path: PathLike,
    *,
    log_rho: bool = True,
    field_name: str = "log10_rho",
) -> Path:
    """Write a 2-D resistivity model as a VTK rectilinear grid.

    The output is plain ASCII VTK (``RECTILINEAR_GRID``) readable by
    Paraview, QGIS (via the SimpleVTK plugin), and most 3-D viewers.

    Parameters
    ----------
    model : ResistivityModel
    path : path-like
        Output ``.vtk`` file.
    log_rho : bool
        Write :math:`\\log_{10}(\\rho)` (default) or linear ρ.
    field_name : str
        VTK scalar name.

    Returns
    -------
    Path
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    x = model.x_centers
    z = model.z_centers
    rho = model.rho_2d if log_rho else 10.0**model.rho_2d

    nx, nz = len(x), len(z)

    with out.open("w") as fh:
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write("pycsamt.interp resistivity model\n")
        fh.write("ASCII\n")
        fh.write("DATASET RECTILINEAR_GRID\n")
        fh.write(f"DIMENSIONS {nx} {nz} 1\n")

        fh.write(f"X_COORDINATES {nx} float\n")
        fh.write(" ".join(f"{v:.3f}" for v in x) + "\n")

        fh.write(f"Y_COORDINATES {nz} float\n")
        fh.write(" ".join(f"{v:.3f}" for v in z) + "\n")

        fh.write("Z_COORDINATES 1 float\n0.0\n")

        n_cells = nx * nz
        fh.write(f"POINT_DATA {n_cells}\n")
        fh.write(f"SCALARS {field_name} float 1\n")
        fh.write("LOOKUP_TABLE default\n")

        for iz in range(nz):
            for ix in range(nx):
                val = rho[iz, ix]
                fh.write(f"{val if not np.isnan(val) else -9999.0:.5f}\n")

    return out


# ---------------------------------------------------------------------------
# Golden Software Surfer
# ---------------------------------------------------------------------------

_SURFER_BLANK = 1.70141e38  # Surfer's own default blanking sentinel


def _elevation_at_x(
    x_centers: np.ndarray,
    elevation: np.ndarray | None,
    chainage: np.ndarray | None,
    fallback_chainage: np.ndarray,
) -> np.ndarray:
    """Interpolate station elevation onto *x_centers*; ``0.0`` everywhere
    when no elevation was supplied (so "Y = elevation − depth" collapses
    to plain "Y = −depth" and every Surfer writer needs only one code
    path)."""
    if elevation is None:
        return np.zeros_like(x_centers)
    chain = np.asarray(
        chainage if chainage is not None else fallback_chainage, dtype=float
    )
    elev = np.asarray(elevation, dtype=float)
    order = np.argsort(chain)
    return np.interp(x_centers, chain[order], elev[order])


def to_surfer_grid(
    model: Any,
    path: PathLike,
    *,
    log_rho: bool = True,
    nx: int | None = None,
    ny: int | None = None,
    elevation: np.ndarray | None = None,
    station_x: np.ndarray | None = None,
    chainage: np.ndarray | None = None,
    blank_value: float = _SURFER_BLANK,
) -> Path:
    """Write a 2-D resistivity model as a Golden Software Surfer DSAA grid.

    Accepts *any* 2-D inversion result via
    :meth:`~pycsamt.interp.ResistivityModel.from_any` (Occam2D, ModEM,
    the backend-neutral ``InversionResult``, an AI agent result, a raw
    ``(x_centers, z_centers, rho_2d)`` triple, or an existing
    :class:`ResistivityModel`) — not just resistivity models already in
    hand. The written file opens directly in Surfer (``File > Open``)
    as a finished grid, no further gridding step needed.

    DSAA is a *regular* grid format: X and Y must each be evenly
    spaced. The source model's cell centres generally are not (Occam2D
    depth cells expand geometrically, mesh padding is uneven, ...), so
    this resamples onto a new ``nx`` x ``ny`` regular grid first. The
    resample is exact and separable when *elevation* is omitted (a
    straightforward 1-D linear interpolation along each axis, since a
    ``ResistivityModel``'s ``rho_2d`` already sits on a true
    ``x_centers`` x ``z_centers`` tensor grid). With *elevation* given,
    each x-column's local depth-to-elevation mapping differs (real
    terrain), so the grid becomes column-sheared before resampling;
    this handles that with one linear interpolation per column rather
    than a full 2-D scattered regrid, still exact along each axis and
    scipy-free. Grid cells outside a column's real elevation coverage
    (above local terrain, or below the deepest column-local sample) are
    written as *blank_value*, Surfer's own missing-data sentinel.

    Parameters
    ----------
    model : object
        Anything :meth:`ResistivityModel.from_any` accepts. Native
        MARE2DEM results (an unstructured triangular mesh, not a
        regular grid) are not supported — regrid onto a
        ``(x_centers, z_centers, rho_2d)`` triple first. A native
        Occam2D/ModEM result is not auto-cropped to the station-carrying
        core — pass ``model.clip_to_stations()`` (after
        ``ResistivityModel.from_any``) instead of the raw result if
        Occam2D's own wide boundary-padding columns should not appear
        in the output.
    path : path-like
        Output ``.grd`` file.
    log_rho : bool
        Write :math:`\\log_{10}(\\rho)` (default) or linear ρ (Ω·m).
    nx, ny : int, optional
        Output grid resolution. Default to the source model's own
        ``n_x``/``n_z`` (a reasonable resolution match; raise for a
        finer resample).
    elevation : ndarray, optional
        Real terrain elevation (m a.s.l.) — pairs with *chainage* (or
        *station_x* when *chainage* is omitted). When given, Y becomes
        real elevation (``elevation − depth``); otherwise Y is
        ``−depth`` (flat datum at 0, matching
        :func:`to_oasis_montaj_xyz`'s convention).
    station_x : ndarray, optional
        Overrides the station positions used to resolve *model* via
        :meth:`ResistivityModel.from_any` (only meaningful for input
        forms that need it, e.g. a raw array triple).
    chainage : ndarray, optional
        Along-profile positions matching *elevation*, metres. Defaults
        to the resolved model's own ``station_x``.
    blank_value : float
        Surfer's missing-data sentinel. Default matches Surfer's own
        default (``1.70141e38``); change only if a project's Surfer
        installation uses a different one.

    Returns
    -------
    Path

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.interp import export
    >>> x = np.linspace(0.0, 900.0, 10)
    >>> z = np.array([10.0, 50.0, 120.0])
    >>> rho_log10 = np.tile(np.array([[2.0], [2.5], [3.0]]), (1, 10))
    >>> path = export.to_surfer_grid(
    ...     (x, z, rho_log10), "surfer_grid_doctest.grd", nx=5, ny=3,
    ... )
    >>> lines = path.read_text().splitlines()
    >>> lines[0]
    'DSAA'
    >>> lines[1]
    '5 3'
    >>> path.unlink()
    """
    rm = ResistivityModel.from_any(model, station_x=station_x)
    values_src = rm.rho_2d if log_rho else 10.0**rm.rho_2d

    x_c = rm.x_centers
    z_c = rm.z_centers
    elev_at_x = _elevation_at_x(x_c, elevation, chainage, rm.station_x)

    n_x = int(nx) if nx else rm.n_x
    n_y = int(ny) if ny else rm.n_z
    x_uniform = np.linspace(float(x_c.min()), float(x_c.max()), n_x)

    # Pass 1: interpolate every original depth row along x onto x_uniform.
    values_x = np.empty((rm.n_z, n_x))
    for iz in range(rm.n_z):
        values_x[iz, :] = np.interp(x_uniform, x_c, values_src[iz, :])
    elev_uniform = np.interp(x_uniform, x_c, elev_at_x)

    y_min = float(np.min(elev_at_x - z_c.max()))
    y_max = float(np.max(elev_at_x - z_c.min()))
    y_uniform = np.linspace(y_min, y_max, n_y)

    # Pass 2: for each x_uniform column, interpolate along that column's
    # own local (elevation − depth) axis onto the shared y_uniform axis,
    # blanking outside its real coverage.
    grid = np.full((n_y, n_x), blank_value)
    for ix in range(n_x):
        y_local = elev_uniform[ix] - z_c
        order = np.argsort(y_local)
        y_sorted = y_local[order]
        col_sorted = values_x[order, ix]
        in_range = (y_uniform >= y_sorted[0]) & (y_uniform <= y_sorted[-1])
        grid[in_range, ix] = np.interp(
            y_uniform[in_range], y_sorted, col_sorted
        )

    # A NaN in the source model (a masked/unresolved cell) propagates
    # through linear interpolation as NaN, not as blank_value -- and
    # Surfer's DSAA reader treats "nan" as an unparsable number, not as
    # its blanking sentinel, so this must be caught explicitly rather
    # than relying on the blank_value fill above.
    grid[np.isnan(grid)] = blank_value

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    finite = grid[grid != blank_value]
    z_min = float(finite.min()) if finite.size else 0.0
    z_max = float(finite.max()) if finite.size else 0.0

    with out.open("w") as fh:
        fh.write("DSAA\n")
        fh.write(f"{n_x} {n_y}\n")
        fh.write(f"{x_uniform.min():.6g} {x_uniform.max():.6g}\n")
        fh.write(f"{y_uniform.min():.6g} {y_uniform.max():.6g}\n")
        fh.write(f"{z_min:.6g} {z_max:.6g}\n")
        for iy in range(n_y):  # row 1 = y_min, matching the DSAA spec
            fh.write(" ".join(f"{v:.6g}" for v in grid[iy, :]) + "\n")

    return out


def to_surfer_xyz(
    model: Any,
    path: PathLike,
    *,
    log_rho: bool = True,
    elevation: np.ndarray | None = None,
    station_x: np.ndarray | None = None,
    chainage: np.ndarray | None = None,
    delimiter: str = "\t",
    header: bool = True,
) -> Path:
    """Write a 2-D resistivity model as scattered Surfer ``X Y Z`` points.

    Unlike :func:`to_surfer_grid`, this writes the model's *actual*
    cell centres exactly as they are — no resampling, no regularity
    requirement, nothing blanked — and lets Surfer's own gridding
    (``Grid > Data``) turn it into a map with whatever method and
    blanking the user chooses. Accepts the same input forms as
    :func:`to_surfer_grid` via :meth:`~pycsamt.interp.ResistivityModel.from_any`.

    Parameters
    ----------
    model : object
        Anything :meth:`ResistivityModel.from_any` accepts.
    path : path-like
        Output text file (``.xyz`` / ``.dat`` recommended).
    log_rho : bool
        Write :math:`\\log_{10}(\\rho)` (default) or linear ρ (Ω·m).
    elevation, station_x, chainage :
        Same meaning as :func:`to_surfer_grid` — with *elevation*, Y is
        real elevation (``elevation − depth``); without it, Y is
        ``−depth``.
    delimiter : str
        Column separator (default tab, Surfer's own default import
        delimiter).
    header : bool
        Write an ``X<delim>Y<delim>Z`` header row.

    Returns
    -------
    Path

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.interp import export
    >>> x = np.array([0.0, 100.0])
    >>> z = np.array([10.0, 50.0])
    >>> rho_log10 = np.array([[2.0, 2.1], [2.5, 2.6]])
    >>> path = export.to_surfer_xyz((x, z, rho_log10), "surfer_xyz_doctest.dat")
    >>> path.read_text().splitlines()
    ['X\\tY\\tZ', '0.000\\t-10.000\\t2.00000', '0.000\\t-50.000\\t2.50000', '100.000\\t-10.000\\t2.10000', '100.000\\t-50.000\\t2.60000']
    >>> path.unlink()
    """
    rm = ResistivityModel.from_any(model, station_x=station_x)
    values = rm.rho_2d if log_rho else 10.0**rm.rho_2d
    x_c = rm.x_centers
    z_c = rm.z_centers
    elev_at_x = _elevation_at_x(x_c, elevation, chainage, rm.station_x)

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as fh:
        if header:
            fh.write(delimiter.join(["X", "Y", "Z"]) + "\n")
        for ix in range(rm.n_x):
            y_col = elev_at_x[ix] - z_c
            for iz in range(rm.n_z):
                val = values[iz, ix]
                if np.isnan(val):
                    continue
                fh.write(
                    delimiter.join(
                        [f"{x_c[ix]:.3f}", f"{y_col[iz]:.3f}", f"{val:.5f}"]
                    )
                    + "\n"
                )

    return out

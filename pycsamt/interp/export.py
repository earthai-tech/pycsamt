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

Example
-------
>>> from pycsamt.interp import export
>>> export.to_oasis_montaj_xyz(logs, "profile_K1.xyz")
>>> export.to_las(logs[0], "S17.las")
>>> export.to_csv(logs, "all_stations.csv")
"""
from __future__ import annotations

import csv
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Union

import numpy as np

from ._base import ResistivityModel
from .lithology import StratigraphicLog

__all__ = [
    "to_oasis_montaj_xyz",
    "to_las",
    "to_csv",
    "to_vtk",
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
            for _iz, (z, rho) in enumerate(
                zip(log.z_centers, log.rho_log10)
            ):
                if np.isnan(rho):
                    continue
                x_val = log.station_x
                z_val = elev - float(z)
                rho_val = float(rho) if log_rho else float(10.0 ** rho)
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
        rho = 10.0 ** rho

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

    fieldnames = ["station", "x_m", "depth_m", "rho_log10", "rho_ohm_m", "lithology"]

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
                writer.writerow({
                    "station":   log.station_name,
                    "x_m":       f"{log.station_x:.3f}",
                    "depth_m":   f"{float(zi):.3f}",
                    "rho_log10": f"{float(rho_l):.5f}",
                    "rho_ohm_m": f"{10.0 ** float(rho_l):.3f}",
                    "lithology": lith,
                })

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
    rho = model.rho_2d if log_rho else 10.0 ** model.rho_2d

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

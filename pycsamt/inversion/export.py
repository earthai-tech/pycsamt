# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Export helpers for :mod:`pycsamt.inversion` results."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .results import InversionResult

PathLike = str | Path

__all__ = ["to_csv", "to_npz"]


def to_csv(result: InversionResult, path: PathLike, *, log_rho: bool = True) -> Path:
    """Export a recovered resistivity model as long-form CSV."""
    model = result.to_resistivity_model()
    rho = np.asarray(model.rho_2d, dtype=float)
    values = rho if log_rho else 10.0 ** rho
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        unit = "log10_ohm_m" if log_rho else "ohm_m"
        writer.writerow(["x_m", "z_m", f"rho_{unit}", "station"])
        for ix, x in enumerate(model.x_centers):
            station = ""
            if ix < len(model.station_names):
                station = model.station_names[ix]
            for iz, z in enumerate(model.z_centers):
                writer.writerow([float(x), float(z), float(values[iz, ix]), station])
    return out


def to_npz(result: InversionResult, path: PathLike) -> Path:
    """Export common result arrays to a compressed NumPy archive."""
    model = result.to_resistivity_model()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        rho_2d=np.asarray(model.rho_2d, dtype=float),
        x_centers=np.asarray(model.x_centers, dtype=float),
        z_centers=np.asarray(model.z_centers, dtype=float),
        station_x=np.asarray(model.station_x, dtype=float),
        station_names=np.asarray(model.station_names, dtype=str),
        method=result.method,
        dimension=result.dimension,
        backend=result.backend,
        rms=result.rms,
    )
    return out

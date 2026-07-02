# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Extract 2-D vertical sections ("curtains") from a ModEM 3-D model.

A ModEM inversion produces one 3-D resistivity volume for the whole
survey — individual survey lines are not separate model files, they
are paths of stations *through* that one volume. :func:`station_curtain`
samples the volume along an ordered list of stations to build the 2-D
(station x depth) resistivity section for one line, so multi-line
surveys inverted as a single 3-D run can still be shown as per-line
fence/depth panels.

Registration between the model grid (cell widths only, no absolute
origin) and the data file's station coordinates (metres, centred near
zero) is derived empirically: ModEM builds the horizontal grid with
symmetric padding around the station footprint, so the grid's
geometric centre always coincides with the midpoint of the station
coordinate range. That gives a closed-form offset without needing to
parse padding/core-zone boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Sequence

import numpy as np

if TYPE_CHECKING:
    from .data import ModEmData
    from .model3d import ModEmModel3D

__all__ = ["ModemSection", "station_curtain"]


@dataclass
class ModemSection:
    """One 2-D vertical resistivity section sliced from a 3-D model."""

    z: np.ndarray
    """Depth below the ground surface (m), earth layers only."""

    rho: np.ndarray
    """Linear resistivity (Ohm.m), shape ``(n_z, n_stations)``."""

    station_names: list[str] = field(default_factory=list)
    """Station names aligned with ``rho``'s second axis."""


def _cell_centers(widths: np.ndarray) -> np.ndarray:
    """Cell-centre coordinates (m) from an array of cell widths."""
    nodes = np.concatenate([[0.0], np.cumsum(widths)])
    return nodes[:-1] + widths / 2.0


def _nearest_index(centers: np.ndarray, value: float) -> int:
    idx = int(np.searchsorted(centers, value))
    if idx <= 0:
        return 0
    if idx >= centers.size:
        return centers.size - 1
    return idx if abs(centers[idx] - value) < abs(centers[idx - 1] - value) else idx - 1


def _grid_offset(widths: np.ndarray, coords: np.ndarray) -> float:
    """Shift that maps *coords* (data frame) onto *widths*' grid frame.

    ModEM pads the horizontal grid symmetrically around the station
    footprint, so the grid's total-extent midpoint always equals the
    midpoint of the station coordinate range.
    """
    finite = coords[np.isfinite(coords)]
    if finite.size == 0:
        return 0.0
    data_mid = (float(finite.min()) + float(finite.max())) / 2.0
    grid_mid = float(np.sum(widths)) / 2.0
    return grid_mid - data_mid


def station_curtain(
    model: "ModEmModel3D",
    data: "ModEmData",
    station_names: Sequence[str],
) -> ModemSection:
    """Sample *model*'s 3-D resistivity volume along *station_names*.

    Parameters
    ----------
    model : ModEmModel3D
        A read 3-D ModEM model (e.g. ``Modular_NLCG_115.rho``).
    data : ModEmData
        The matching data file, supplying ``site_coords`` (metres,
        same run) used to locate each station in the model grid.
    station_names : sequence of str
        Ordered station names for one survey line. Stations missing
        from ``data.site_coords`` are silently skipped.

    Returns
    -------
    ModemSection
        Depth axis (earth layers only, air layers excluded) and a
        ``(n_z, n_kept)`` linear-resistivity curtain, nearest-cell
        sampled — no smoothing across the model's own discretisation.
    """
    x_centers = _cell_centers(model.x_widths)
    y_centers = _cell_centers(model.y_widths)

    z_widths_earth = np.asarray(model.z_widths[model.n_air:], dtype=float)
    z_nodes = np.concatenate([[0.0], np.cumsum(z_widths_earth)])
    z_centers = z_nodes[:-1] + z_widths_earth / 2.0

    rho_earth = model.rho_linear[model.n_air:, :, :]

    offset_x = _grid_offset(model.x_widths, data.x_coords)
    offset_y = _grid_offset(model.y_widths, data.y_coords)

    names: list[str] = []
    columns: list[np.ndarray] = []
    for name in station_names:
        xyz = data.site_coords.get(name)
        if xyz is None:
            continue
        x_m, y_m = float(xyz[0]), float(xyz[1])
        ix = _nearest_index(x_centers, x_m + offset_x)
        iy = _nearest_index(y_centers, y_m + offset_y)
        columns.append(rho_earth[:, iy, ix])
        names.append(name)

    rho = (
        np.column_stack(columns)
        if columns
        else np.empty((z_centers.size, 0), dtype=float)
    )
    return ModemSection(z=z_centers, rho=rho, station_names=names)

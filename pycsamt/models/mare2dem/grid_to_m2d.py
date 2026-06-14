# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Convert a 2-D resistivity grid to a MARE2DEM model.

Port of ``m2d_gridToM2D.m``.

:func:`grid_to_mare2dem` takes a regular 2-D resistivity grid (y, z,
ρ) and writes the corresponding MARE2DEM ``.poly``, ``.resistivity``,
and ``.settings`` files to a specified output directory.

The grid is assumed to be cell-centred.  Boundary node coordinates are
derived from the cell spacings, a padding region with the log-mean
resistivity is added around the grid, and an air layer (10¹² Ω·m) is
placed above the model.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .iotools.poly import PolyFile, write_poly
from .iotools.resistivity import ResistivityFile, write_resistivity
from .iotools.settings import SettingsFile, write_settings

__all__ = ["grid_to_mare2dem"]


def _make_nodes_segs_regions(
    Y: np.ndarray,
    Z: np.ndarray,
    Rho: np.ndarray,
    padding_y: float,
    padding_z: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build mesh nodes, segments, region seeds, and resistivities.

    Returns
    -------
    nodes : (n_nodes, 2) float array – (y, z) coordinates.
    segs  : (n_segs, 3) int array – (i, j, boundary_marker).
    regions : (n_regions, 2) float array – seed points.
    res   : (n_regions,) float array – resistivity per region (linear Ω·m).
    """
    nZ_c, nY_c = Y.shape   # number of cell-centre rows and columns

    # 1D cell-centre arrays (assume uniform within rows/columns)
    y_c = Y[0, :]           # shape (nY_c,)
    z_c = Z[:, 0]           # shape (nZ_c,)

    # 1D cell-edge arrays
    def _edges(centers: np.ndarray) -> np.ndarray:
        d = np.diff(centers)
        e = np.empty(len(centers) + 1)
        e[0] = centers[0] - d[0] / 2.0
        e[1:-1] = (centers[:-1] + centers[1:]) / 2.0
        e[-1] = centers[-1] + d[-1] / 2.0
        return e

    y_e = _edges(y_c)   # length nY_c + 1
    z_e = _edges(z_c)   # length nZ_c + 1

    nY = len(y_e)   # nY_c + 1
    nZ = len(z_e)   # nZ_c + 1

    # Node index: iz * nY + iy  (0-based)
    Ye, Ze = np.meshgrid(y_e, z_e)   # both (nZ, nY)
    nodes = np.column_stack([Ye.ravel(), Ze.ravel()])  # (nZ*nY, 2)

    # node index function (1-based)
    def nidx(iz: int, iy: int) -> int:
        return iz * nY + iy + 1

    # horizontal segments: row iz, columns iy..iy+1
    segs_h = [
        [nidx(iz, iy), nidx(iz, iy + 1), 1]
        for iz in range(nZ)
        for iy in range(nY - 1)
    ]

    # vertical segments: column iy, rows iz..iz+1
    segs_v = [
        [nidx(iz, iy), nidx(iz + 1, iy), 1]
        for iz in range(nZ - 1)
        for iy in range(nY)
    ]

    segs_grid = np.array(segs_h + segs_v, dtype=int)

    min_y = nodes[:, 0].min()
    max_y = nodes[:, 0].max()
    min_z = nodes[:, 1].min()
    max_z = nodes[:, 1].max()

    # bounding-box nodes (6 extra nodes, 1-based offsets)
    box_nodes = np.array([
        [min_y - padding_y, min_z - padding_z],
        [min_y - padding_y, min_z],
        [min_y - padding_y, max_z + padding_z],
        [max_y + padding_y, max_z + padding_z],
        [max_y + padding_y, min_z],
        [max_y + padding_y, min_z - padding_z],
    ])
    all_nodes = np.vstack([box_nodes, nodes])
    n_box = 6
    # renumber grid segs
    segs_grid[:, :2] += n_box

    # boundary box segs
    box_segs = np.array([
        [1, 2, 1], [2, 3, 1], [3, 4, 1],
        [4, 5, 1], [5, 6, 1], [6, 1, 1],
    ], dtype=int)

    # connect sides to top of central grid
    top_left_grid = n_box + nidx(0, 0)           # iz=0, iy=0
    top_right_grid = n_box + nidx(0, nY - 1)     # iz=0, iy=nY-1
    connect_segs = np.array([
        [2, top_left_grid, 1],
        [5, top_right_grid, 1],
    ], dtype=int)

    all_segs = np.vstack([box_segs, connect_segs, segs_grid])

    # region seed points: cell centres (Y, Z are already cell centres)
    regions = np.column_stack([Y.ravel(), Z.ravel()])
    res = Rho.ravel().copy()

    # add air and ground padding regions
    mean_rho = 10.0 ** np.mean(np.log10(np.maximum(Rho, 1e-10)))
    extra_seeds = np.array([
        [min_y - padding_y + 0.1, min_z - padding_z + 0.1],  # air
        [min_y - padding_y + 0.1, min_z + 0.1],              # ground padding
    ])
    regions = np.vstack([regions, extra_seeds])
    res = np.concatenate([res, [1e12, mean_rho]])

    return all_nodes, all_segs, regions, res


def grid_to_mare2dem(
    Y: np.ndarray,
    Z: np.ndarray,
    Rho: np.ndarray,
    *,
    padding_y: float = 50000.0,
    padding_z: float = 50000.0,
    out_dir: str | Path = ".",
    model_name: str = "mare2dem",
    data_file: str = "mare2dem.emdata",
    target_misfit: float = 1.0,
    max_iterations: int = 100,
) -> dict[str, Path]:
    """Create a MARE2DEM model from a regular 2-D resistivity grid.

    Port of ``m2d_gridToM2D.m``.

    Parameters
    ----------
    Y : array-like, shape (nz, ny)
        y-coordinates (along-profile) of grid cell centres in metres.
        Values must vary along columns; all rows for the same column
        share the same y value.
    Z : array-like, shape (nz, ny)
        Depth coordinates of grid cell centres in metres (positive down).
        Values must vary along rows; all columns for the same row share
        the same z value.
    Rho : array-like, shape (nz, ny)
        Resistivity at each cell centre in Ω·m.
    padding_y : float, default 50000.0
        Lateral padding in metres added outside the grid.
    padding_z : float, default 50000.0
        Vertical padding in metres added below and above the grid.
    out_dir : path-like, default "."
        Output directory.
    model_name : str, default "mare2dem"
        Stem name used for all output files.
    data_file : str, default "mare2dem.emdata"
        Name of the associated data file written into the model header.
    target_misfit : float, default 1.0
        Target normalized RMS misfit.
    max_iterations : int, default 100
        Maximum inversion iterations.

    Returns
    -------
    dict[str, pathlib.Path]
        Keys ``"resistivity"``, ``"poly"``, ``"settings"``.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.grid_to_m2d import grid_to_mare2dem
    >>> y1d = np.linspace(-5000, 5000, 21)
    >>> z1d = np.linspace(0, 3000, 11)
    >>> Y, Z = np.meshgrid(y1d, z1d)
    >>> Rho = np.ones_like(Y) * 10.0      # 10 Ω·m half-space
    >>> files = grid_to_mare2dem(Y, Z, Rho, out_dir="/tmp/m2d_grid_test")
    >>> files["resistivity"].exists()
    True
    """
    Y_arr = np.asarray(Y, dtype=float)
    Z_arr = np.asarray(Z, dtype=float)
    Rho_arr = np.asarray(Rho, dtype=float)

    dest = Path(out_dir)
    dest.mkdir(parents=True, exist_ok=True)

    nodes, segs, regions, res = _make_nodes_segs_regions(
        Y_arr, Z_arr, Rho_arr, padding_y, padding_z
    )

    # ---- .poly file ----
    pf = PolyFile()
    pf.nodes = nodes
    pf.segments = segs[:, :2].astype(int)
    pf.segment_markers = segs[:, 2].astype(int)
    pf.regions = np.column_stack([
        regions,
        np.arange(1, len(regions) + 1, dtype=float),
        -np.ones(len(regions), dtype=float),
    ])
    poly_path = write_poly(pf, dest / f"{model_name}.poly")

    # ---- .resistivity file ----
    rf = ResistivityFile(
        resistivity_file=str(dest / f"{model_name}.0.resistivity"),
        poly_file=f"{model_name}.poly",
        data_file=data_file,
        settings_file="mare2dem.settings",
        anisotropy="isotropic",
        target_misfit=target_misfit,
        max_iterations=max_iterations,
        global_bounds=np.array([-1.0, 5.0]),
        roughness_penalty_method="gradient",
        yz_penalty_weights=np.array([3.0, 1.0]),
        penalty_cut_weight=0.1,
    )
    n_reg = len(res)
    rf.resistivity = res.reshape(n_reg, 1)
    rf.free_parameter = np.ones((n_reg, 1), dtype=float)
    # air region and padding: fix them
    rf.free_parameter[-2:] = 0
    rf.bounds = np.tile([-1.0, 5.0], (n_reg, 1))
    rf.prejudice = np.zeros((n_reg, 2))
    res_path = write_resistivity(rf, dest / f"{model_name}.0.resistivity")

    # ---- .settings file ----
    sf = SettingsFile()
    settings_path = write_settings(sf, dest / "mare2dem.settings")

    return {
        "resistivity": res_path,
        "poly": poly_path,
        "settings": settings_path,
    }

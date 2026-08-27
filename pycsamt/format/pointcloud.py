# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Generic point-cloud extraction from any :class:`PCSFModel` geometry kind.

Phase 7 of the PCSF format plan needs a 3-D view that works the same
way regardless of which backend produced the file — "no backend-
specific glue code in the view layer" is the phase's own definition of
done. :func:`pcsf_to_point_cloud` is that one shared extraction: it
knows how to turn each of the four geometry kinds
(:data:`~pycsamt.format.schema.GEOMETRY_KINDS`) into a flat
``(x, y, z, log10_rho)`` point cloud, so a view (the desktop 3-D panel
in this phase, potentially others later) only has to know how to
scatter-plot points, never how Occam2D/ModEM/MARE2DEM/a multiline
stack each store their own geometry.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .schema import PCSFModel

__all__ = ["PointCloud", "pcsf_to_point_cloud"]


@dataclass(frozen=True)
class PointCloud:
    """Flat point cloud ready for a 3-D scatter view.

    Attributes
    ----------
    x, y, z : ndarray, shape (n,)
        Position, metres. ``z`` is elevation-like (positive up) —
        callers plotting depth sections see negative values below the
        surface, matching the sign convention already used by
        ``pycsamt.app.web.callbacks.map3d``'s own 3-D views.
    value : ndarray, shape (n,)
        :math:`\\log_{10}(\\rho / \\Omega\\mathrm{m})`.
    label : str
        Short description of what was plotted (geometry kind + any
        subsampling applied), for a status bar / axis title.
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    value: np.ndarray
    label: str


def _log10_rho(resistivity: np.ndarray) -> np.ndarray:
    return np.log10(np.clip(np.asarray(resistivity, dtype=float), 1e-3, None))


def _finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(arrays[0].shape, dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


def _subsample(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    value: np.ndarray,
    max_points: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool]:
    n = x.size
    if n <= max_points:
        return x, y, z, value, False
    rng = np.random.default_rng(seed)
    idx = rng.choice(n, size=max_points, replace=False)
    idx.sort()
    return x[idx], y[idx], z[idx], value[idx], True


def _grid2d_points(model: PCSFModel) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    geo = model.geometry
    x2d, z2d = np.meshgrid(geo.x, geo.z)  # both (n_z, n_x)
    y2d = np.zeros_like(x2d)
    value2d = _log10_rho(model.resistivity)
    return x2d.ravel(), y2d.ravel(), -z2d.ravel(), value2d.ravel()


def _grid3d_points(model: PCSFModel) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    geo = model.geometry
    # resistivity is (n_z, n_y, n_x); build matching (n_z, n_y, n_x) coord
    # grids directly rather than transposing, so no axis-order mistake
    # can creep in between the two.
    z3d, y3d, x3d = np.meshgrid(geo.z, geo.y, geo.x, indexing="ij")
    value3d = _log10_rho(model.resistivity)
    return x3d.ravel(), y3d.ravel(), -z3d.ravel(), value3d.ravel()


def _mesh_unstructured_points(
    model: PCSFModel,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    geo = model.geometry
    centroids = geo.nodes[geo.connectivity].mean(axis=1)  # (n_tri, 2)
    if model.resistivity.shape[0] != geo.connectivity.shape[0]:
        # PCSFModel.validate() also accepts the compact, region-count
        # form of resistivity for mesh_unstructured (shape (n_regions,));
        # a point per triangle needs the per-cell expansion instead —
        # every adapter in this package already produces that (e.g.
        # mare2dem_to_pcsf's resistivity_by_region field carries the
        # compact table separately), so this only fires for a
        # hand-built model that skipped the expansion.
        raise ValueError(
            "pcsf_to_point_cloud needs per-cell mesh_unstructured "
            f"resistivity (shape ({geo.connectivity.shape[0]},)), got "
            f"shape {model.resistivity.shape} — expand "
            "resistivity_by_region onto each triangle's region id first."
        )
    value = _log10_rho(model.resistivity)
    # A MARE2DEM-style mesh lives in one vertical (x, z) plane; y is a
    # single cross-strike position (0.0 -- there is no second line to
    # offset against for a lone unstructured mesh).
    x = centroids[:, 0]
    z = centroids[:, 1]
    y = np.zeros_like(x)
    return x, y, -z, value


def _multiline_points(
    model: PCSFModel,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xs, ys, zs, values = [], [], [], []
    for line in model.geometry.lines:
        x2d, z2d = np.meshgrid(line.geometry.x, line.geometry.z)
        y2d = np.full_like(x2d, line.offset_y)
        value2d = _log10_rho(line.resistivity)
        xs.append(x2d.ravel())
        ys.append(y2d.ravel())
        zs.append(-z2d.ravel())
        values.append(value2d.ravel())
    return (
        np.concatenate(xs),
        np.concatenate(ys),
        np.concatenate(zs),
        np.concatenate(values),
    )


_EXTRACTORS = {
    "grid2d": _grid2d_points,
    "grid3d": _grid3d_points,
    "mesh_unstructured": _mesh_unstructured_points,
    "multiline": _multiline_points,
}


def pcsf_to_point_cloud(
    model: PCSFModel,
    *,
    max_points: int = 200_000,
    seed: int = 0,
) -> PointCloud:
    """Flatten any :class:`PCSFModel` geometry into one 3-D point cloud.

    Parameters
    ----------
    model : PCSFModel
        Any geometry kind.
    max_points : int, default 200_000
        Random (seeded, reproducible) subsample cap — a native
        ``grid3d``/``mesh_unstructured`` model can carry hundreds of
        thousands of cells, too many for an interactive scatter plot.
    seed : int, default 0
        Subsampling RNG seed, for a reproducible view across renders.

    Returns
    -------
    PointCloud
        Non-finite values (masked cells, log of non-positive
        resistivity) are dropped, not zeroed.

    Raises
    ------
    ValueError
        If ``model.geometry.kind`` is not one of
        :data:`~pycsamt.format.schema.GEOMETRY_KINDS`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format import Grid2DGeometry, PCSFModel
    >>> from pycsamt.format.pointcloud import pcsf_to_point_cloud
    >>> geometry = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
    >>> model = PCSFModel(geometry=geometry, resistivity=np.array([[100.0, 110.0], [50.0, 55.0]]))
    >>> cloud = pcsf_to_point_cloud(model)
    >>> cloud.x.shape
    (4,)
    """
    kind = model.geometry.kind
    extractor = _EXTRACTORS.get(kind)
    if extractor is None:
        raise ValueError(
            f"unsupported geometry kind {kind!r}; expected one of "
            f"{tuple(_EXTRACTORS)}"
        )
    x, y, z, value = extractor(model)
    keep = _finite_mask(x, y, z, value)
    x, y, z, value = x[keep], y[keep], z[keep], value[keep]

    x, y, z, value, subsampled = _subsample(x, y, z, value, max_points, seed)
    label = f"{kind} ({x.size} points{', subsampled' if subsampled else ''})"
    return PointCloud(x=x, y=y, z=z, value=value, label=label)

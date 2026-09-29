# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Horizontal slices through a 3-D inversion model (PCSF ``grid3d``).

:func:`model_depth_slice` cuts the resistivity volume at a depth (or
averages it from the surface down to that depth), crops it to the station
area -- the padding cells of a ModEM model reach hundreds of kilometres --
and, when the stations carry longitude/latitude, geo-references the cell
grid with an affine fit between the stations' local model coordinates and
their geographic positions.  It is the real model, not an interpolation of
per-station values.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ModelSlice", "model_depth_slice"]


@dataclass
class ModelSlice:
    """A horizontal cut through a 3-D model.

    ``x``/``y`` are cell-edge coordinates of shape ``(ny + 1, nx + 1)``:
    longitude/latitude when ``geo`` is True, else local model metres.
    ``rho`` has shape ``(ny, nx)`` (Ω·m, NaN outside the model).
    """

    x: np.ndarray
    y: np.ndarray
    rho: np.ndarray
    depth: float
    mode: str  # "at" | "mean"
    geo: bool


def _affine(src: np.ndarray, dst: np.ndarray) -> np.ndarray | None:
    """3x2 least-squares affine map ``[x, y, 1] @ A -> dst``."""
    ok = np.isfinite(src).all(axis=1) & np.isfinite(dst).all(axis=1)
    if ok.sum() < 3:
        return None
    A = np.column_stack([src[ok], np.ones(ok.sum())])
    if np.linalg.matrix_rank(A) < 3:  # collinear stations: no 2-D fit
        return None
    coef, *_ = np.linalg.lstsq(A, dst[ok], rcond=None)
    return coef


def model_depth_slice(
    model,
    depth: float,
    *,
    mode: str = "at",
    margin: float = 0.25,
    georeference: bool = True,
) -> ModelSlice:
    """Resistivity of a ``grid3d`` model at (or down to) *depth*.

    Parameters
    ----------
    model : PCSFModel
        A model read with :func:`pycsamt.format.text.read_pcsf_or_pcsm`.
    depth : float
        Metres below the model top (PCSF ``z`` is depth, positive down).
    mode : {"at", "mean"}
        ``"at"``: log-linear interpolation between the two cell layers
        around *depth*.  ``"mean"``: log-mean of every layer from the top
        down to *depth*.
    margin : float
        Extra area around the stations, as a fraction of their extent.
    georeference : bool
        Map the cells to longitude/latitude when the stations allow it.

    Raises
    ------
    ValueError
        Not a ``grid3d`` model, or *depth* outside the model.
    """
    g = model.geometry
    if getattr(model, "kind", None) != "grid3d" and \
            type(g).__name__ != "Grid3DGeometry":
        raise ValueError("a depth slice needs a 3-D (grid3d) model")
    rho = np.asarray(model.resistivity, dtype=float)  # (nz, ny, nx)
    z = np.asarray(g.z, dtype=float)
    if not (z.min() <= depth <= z.max()) and mode == "at":
        raise ValueError(f"{depth:g} m is outside the model "
                         f"({z.min():g}-{z.max():g} m)")
    logr = np.log10(np.where(rho > 0, rho, np.nan))
    if mode == "mean":
        sel = z <= max(depth, z.min())
        with np.errstate(invalid="ignore"):
            layer = np.nanmean(logr[sel], axis=0)
    else:
        k = int(np.clip(np.searchsorted(z, depth), 1, z.size - 1))
        w = (depth - z[k - 1]) / max(z[k] - z[k - 1], 1e-12)
        layer = (1 - w) * logr[k - 1] + w * logr[k]

    ox, oy = float(g.origin[0]), float(g.origin[1])
    xe = np.asarray(g.x_nodes, float) + ox
    ye = np.asarray(g.y_nodes, float) + oy
    # Stations share the model frame (grid x/y + origin)
    st = model.stations
    sx = np.asarray(st.x, dtype=float) if st is not None else np.array([])
    sy = np.asarray(st.y, dtype=float) if st is not None else np.array([])
    good = np.isfinite(sx) & np.isfinite(sy)
    if good.any():
        span = max(np.ptp(sx[good]), np.ptp(sy[good]), 1.0)
        lo_x, hi_x = sx[good].min() - margin * span, sx[good].max() + \
            margin * span
        lo_y, hi_y = sy[good].min() - margin * span, sy[good].max() + \
            margin * span
        ix = np.flatnonzero((xe[1:] >= lo_x) & (xe[:-1] <= hi_x))
        iy = np.flatnonzero((ye[1:] >= lo_y) & (ye[:-1] <= hi_y))
        if ix.size and iy.size:
            layer = layer[iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
            xe = xe[ix[0]:ix[-1] + 2]
            ye = ye[iy[0]:iy[-1] + 2]
    X, Y = np.meshgrid(xe, ye)
    geo = False
    if georeference and st is not None and st.lat is not None \
            and st.lon is not None:
        src = np.column_stack([sx, sy])
        dst = np.column_stack([np.asarray(st.lon, float),
                               np.asarray(st.lat, float)])
        coef = _affine(src, dst)
        if coef is not None:
            P = np.column_stack([X.ravel(), Y.ravel(), np.ones(X.size)])
            LL = P @ coef
            X = LL[:, 0].reshape(X.shape)
            Y = LL[:, 1].reshape(Y.shape)
            geo = True
    return ModelSlice(x=X, y=Y, rho=10.0 ** layer, depth=float(depth),
                      mode=mode, geo=geo)

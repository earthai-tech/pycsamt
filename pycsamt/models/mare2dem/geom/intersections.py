# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Line-segment intersection utilities.

Port of ``m2d_getIntersections.m``.

Tests whether any segment from a set of segments *a* intersects a
single query segment *b*, and returns intersection coordinates and
parametric positions.
"""

from __future__ import annotations

import numpy as np

__all__ = ["get_intersections", "do_rects_overlap"]


def do_rects_overlap(
    a_rect: np.ndarray,
    many_rects: np.ndarray,
    *,
    tol: float = 100.0 * np.finfo(float).eps,
) -> np.ndarray:
    """Return indices of rows in *many_rects* that overlap *a_rect*.

    Port of the inner ``doRectsOverlap`` function in
    ``m2d_getIntersections.m``.

    Parameters
    ----------
    a_rect : array-like, shape (4,)
        Single bounding rectangle ``[x0, x1, y0, y1]``.
    many_rects : array-like, shape (n, 4)
        Array of bounding rectangles, same format.
    tol : float
        Tolerance for the overlap test.

    Returns
    -------
    numpy.ndarray of int
        Row indices of *many_rects* whose bounding boxes overlap
        *a_rect*.
    """
    a = np.asarray(a_rect, dtype=float)
    b = np.asarray(many_rects, dtype=float)
    if b.ndim == 1:
        b = b.reshape(1, -1)

    ax1, ax2 = min(a[0], a[1]), max(a[0], a[1])
    ay1, ay2 = min(a[2], a[3]), max(a[2], a[3])

    bx1 = np.minimum(b[:, 0], b[:, 1])
    bx2 = np.maximum(b[:, 0], b[:, 1])
    by1 = np.minimum(b[:, 2], b[:, 3])
    by2 = np.maximum(b[:, 2], b[:, 3])

    no_overlap = (
        (ax1 - bx2 > tol)
        | (bx1 - ax2 > tol)
        | (ay1 - by2 > tol)
        | (by1 - ay2 > tol)
    )
    return np.where(~no_overlap)[0]


def get_intersections(
    xya: np.ndarray,
    xyb: np.ndarray,
    *,
    tol: float = 1000.0 * np.finfo(float).eps,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Find intersections between segments in *xya* and segment *xyb*.

    Port of ``m2d_getIntersections.m``.

    Parameters
    ----------
    xya : array-like, shape (n_a, 4)
        Set of line segments. Each row is ``[x0, x1, y0, y1]``.
    xyb : array-like, shape (1, 4) or (4,)
        Single query segment ``[x0, x1, y0, y1]``.
    tol : float
        Interior-point tolerance (intersection on exact endpoint is
        excluded).

    Returns
    -------
    intersect : numpy.ndarray of int
        Original row indices in *xya* where an interior intersection
        with *xyb* exists.
    xi : numpy.ndarray of float
        x coordinates of intersection points.
    yi : numpy.ndarray of float
        y coordinates of intersection points.
    pa : numpy.ndarray of float, shape (n_a,)
        Parametric position along each segment in *xya*
        (``-1`` for segments not tested).
    pb : numpy.ndarray of float, shape (n_a,)
        Parametric position along *xyb* (``-1`` for untested).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.geom.intersections import get_intersections
    >>> segs_a = np.array([[0.0, 2.0, 1.0, 1.0]])   # horizontal segment
    >>> seg_b  = np.array([1.0, 1.0, 0.0, 2.0])     # vertical segment
    >>> inter, xi, yi, pa, pb = get_intersections(segs_a, seg_b)
    >>> xi
    array([1.])
    >>> yi
    array([1.])
    """
    xya = np.asarray(xya, dtype=float)
    if xya.ndim == 1:
        xya = xya.reshape(1, -1)
    xyb = np.asarray(xyb, dtype=float).ravel()
    na_input = len(xya)

    pa_full = -np.ones(na_input)
    pb_full = -np.ones(na_input)

    # Bounding-box pre-filter
    overlap_idx = do_rects_overlap(xyb, xya, tol=tol)
    if len(overlap_idx) == 0:
        return np.array([], dtype=int), np.array([]), np.array([]), pa_full, pb_full

    sub = xya[overlap_idx]
    na = len(sub)
    xyb_rep = np.tile(xyb, (na, 1))

    dxb = xyb_rep[:, 1] - xyb_rep[:, 0]
    dxa = sub[:, 1] - sub[:, 0]
    dyb = xyb_rep[:, 3] - xyb_rep[:, 2]
    dya = sub[:, 3] - sub[:, 2]
    dy1ab = sub[:, 2] - xyb_rep[:, 2]
    dx1ab = sub[:, 0] - xyb_rep[:, 0]

    num_a = dxb * dy1ab - dyb * dx1ab
    num_b = dxa * dy1ab - dya * dx1ab
    den = dyb * dxa - dxb * dya

    with np.errstate(divide="ignore", invalid="ignore"):
        pa_sub = np.where(den != 0, num_a / den, np.inf)
        pb_sub = np.where(den != 0, num_b / den, np.inf)

    xi_sub = sub[:, 0] + dxa * pa_sub
    yi_sub = sub[:, 2] + dya * pa_sub

    interior = (
        (pa_sub > 0 + tol) & (pa_sub < 1 - tol)
        & (pb_sub > 0 + tol) & (pb_sub < 1 - tol)
    )

    pa_full[overlap_idx] = pa_sub
    pb_full[overlap_idx] = pb_sub

    inter_local = np.where(interior)[0]
    intersect = overlap_idx[inter_local]
    xi = xi_sub[inter_local]
    yi = yi_sub[inter_local]

    return intersect, xi, yi, pa_full, pb_full

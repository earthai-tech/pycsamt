# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Terrain-following coordinate transform for 2-D section plots.

A standard 2-D resistivity section stores:

    * ``x`` — along-profile distance (km)
    * ``z`` — depth below a *flat datum* at z = 0 (km, positive downward)

When real surface topography is available, the datum is no longer flat.
This module transforms the flat depth-grid into a *terrain-following*
coordinate frame in which z = 0 at every x position corresponds to the
local surface elevation:

    z_real(x, z) = elev(x) - z                  [both in km]

The result can be fed directly to :func:`matplotlib.pyplot.pcolormesh`
using a 2-D ``Z`` argument (supported since Matplotlib 3.3) so that the
mesh cells drape over the terrain.

Typical usage::

    from pycsamt.topo.drape import interp_elev, drape_section
    elev_at_x = interp_elev(chain_km, elev_km, x_centres_km)
    x_nodes, z_draped, data = drape_section(x_nodes, z_nodes, rho_2d, elev_at_x)
    ax.pcolormesh(x_nodes, z_draped, data, ...)
"""

from __future__ import annotations

import numpy as np

try:
    from scipy.interpolate import interp1d as _scipy_interp1d
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

__all__ = [
    "interp_elev",
    "drape_section",
    "mask_above_topo",
    "station_surface_z",
]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def interp_elev(
    chainage_km: np.ndarray,
    elev_km: np.ndarray,
    x_query_km: np.ndarray,
    method: str = "linear",
) -> np.ndarray:
    """Interpolate station elevations to arbitrary profile positions.

    Clamps extrapolated values to the boundary station elevations so
    the terrain surface never shoots up unexpectedly at the section edges.

    Parameters
    ----------
    chainage_km : array_like, shape (n_stations,)
        Along-profile distances of the stations (km).
    elev_km : array_like, shape (n_stations,)
        Terrain elevation at each station (km a.s.l.).
    x_query_km : array_like, shape (m,)
        Positions at which to evaluate the interpolated elevation (km).
    method : {"linear", "cubic", "nearest"}
        Interpolation method.  ``"cubic"`` requires scipy.

    Returns
    -------
    numpy.ndarray, shape (m,)
        Interpolated elevation in **km** a.s.l.
    """
    chainage_km = np.asarray(chainage_km, dtype=float).ravel()
    elev_km     = np.asarray(elev_km,     dtype=float).ravel()
    x_query_km  = np.asarray(x_query_km,  dtype=float).ravel()

    if len(chainage_km) == 0:
        return np.zeros_like(x_query_km)
    if len(chainage_km) == 1:
        return np.full_like(x_query_km, elev_km[0])

    # Remove NaN / Inf elevation values before interpolating
    valid = np.isfinite(elev_km)
    c = chainage_km[valid]
    e = elev_km[valid]
    if len(c) < 2:
        return np.full_like(x_query_km, e[0] if len(e) else 0.0)

    if _HAS_SCIPY and method == "cubic":
        fn = _scipy_interp1d(c, e, kind="cubic",
                             bounds_error=False,
                             fill_value=(e[0], e[-1]))
        return fn(x_query_km)

    # numpy linear / nearest
    if method == "nearest":
        idx = np.argmin(
            np.abs(x_query_km[:, np.newaxis] - c[np.newaxis, :]), axis=1
        )
        return e[idx]

    # default: linear with clamped extrapolation
    return np.interp(x_query_km, c, e, left=e[0], right=e[-1])


def drape_section(
    x_nodes: np.ndarray,
    z_nodes: np.ndarray,
    data: np.ndarray,
    elev_at_centres: np.ndarray,
    exaggeration: float = 1.0,
    clip_above_surface: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transform a flat depth section into terrain-following coordinates.

    Builds a 2-D ``z_draped`` array (shape ``(nz+1, nx+1)``) where each
    column is shifted so that ``z_nodes[0]`` (the surface) aligns with
    the local terrain elevation.  The result can be passed directly to
    :func:`~matplotlib.axes.Axes.pcolormesh` as the ``Y`` argument.

    Parameters
    ----------
    x_nodes : array_like, shape (nx+1,)
        Horizontal pcolormesh node positions (km).
    z_nodes : array_like, shape (nz+1,)
        Depth node positions (km, positive downward from flat datum).
        ``z_nodes[0]`` should be 0 (surface) or the shallowest depth.
    data : array_like, shape (nz, nx)
        2-D data values (e.g. log10(rho)).
    elev_at_centres : array_like, shape (nx,)
        Terrain elevation at each *cell-centre* x position (km a.s.l.).
    exaggeration : float
        Vertical exaggeration applied to both the elevation offset and
        the depth axis.  Values > 1 amplify relief for display purposes.
    clip_above_surface : bool
        If ``True``, set cells that are above the terrain surface
        (unreachable subsurface) to ``NaN``.

    Returns
    -------
    x_nodes : numpy.ndarray, shape (nx+1,)
        Unchanged horizontal node positions.
    z_draped : numpy.ndarray, shape (nz+1, nx+1)
        2-D elevation node array.  Column *j* equals
        ``elev_at_nodes[j] - z_nodes * exaggeration``.
    data_out : numpy.ndarray, shape (nz, nx)
        Input data, optionally NaN-masked above terrain.

    Notes
    -----
    The terrain-following coordinate for node ``(k, j)`` is::

        z_draped[k, j] = elev_nodes[j] - z_nodes[k] * exaggeration

    where ``elev_nodes`` is the terrain interpolated to the x *node*
    positions (``nx+1`` values) from the cell-centre values.
    """
    x_nodes          = np.asarray(x_nodes,          dtype=float)
    z_nodes          = np.asarray(z_nodes,          dtype=float)
    data             = np.asarray(data,             dtype=float)
    elev_at_centres  = np.asarray(elev_at_centres,  dtype=float)

    nx = x_nodes.shape[0] - 1
    nz = z_nodes.shape[0] - 1

    # Interpolate elevation from cell centres to node positions (nx+1)
    if nx > 0:
        x_centres = (x_nodes[:-1] + x_nodes[1:]) / 2.0
        elev_nodes = np.interp(x_nodes, x_centres, elev_at_centres,
                               left=elev_at_centres[0],
                               right=elev_at_centres[-1])
    else:
        elev_nodes = np.full(nx + 1, elev_at_centres[0] if len(elev_at_centres) else 0.0)

    # Build 2-D draped z grid  (nz+1, nx+1)
    # Broadcasting: rows = depth nodes, cols = x nodes
    z_draped = (elev_nodes[np.newaxis, :]                    # (1, nx+1)
                - z_nodes[:, np.newaxis] * exaggeration)     # (nz+1, 1)

    data_out = data.copy()
    if clip_above_surface and nx > 0 and nz > 0:
        data_out = mask_above_topo(x_nodes, z_nodes, data_out,
                                   elev_at_centres, exaggeration)

    return x_nodes, z_draped, data_out


def mask_above_topo(
    x_nodes: np.ndarray,
    z_nodes: np.ndarray,
    data: np.ndarray,
    elev_at_centres: np.ndarray,
    exaggeration: float = 1.0,
) -> np.ndarray:
    """Set data cells that lie above the terrain surface to NaN.

    A cell at column *j* and depth-row *k* has absolute elevation::

        cell_elev = elev_at_centres[j] - z_centre[k] * exaggeration

    where ``z_centre[k] = (z_nodes[k] + z_nodes[k+1]) / 2``.

    A cell is *above the surface* when ``cell_elev > elev_at_centres[j]``,
    which simplifies to ``z_centre[k] < 0``.  For standard meshes
    (all z ≥ 0) this never occurs; the masking is only relevant for
    meshes whose z-origin is below the deepest station (uncommon).

    The more practically useful masking is for *varying terrain*: station
    A sits at 800 m, station B sits at 200 m.  A cell at depth 500 m
    below A has absolute elevation 300 m — it is above the surface at B.
    This function masks such cells so they do not appear in the plot.

    Parameters
    ----------
    x_nodes : (nx+1,) node positions
    z_nodes : (nz+1,) depth nodes (km, positive down)
    data : (nz, nx) data array
    elev_at_centres : (nx,) elevation at cell centres (km)
    exaggeration : float

    Returns
    -------
    numpy.ndarray, shape (nz, nx)
        Data with cells above the terrain set to NaN.
    """
    z_nodes         = np.asarray(z_nodes,          dtype=float)
    elev_at_centres = np.asarray(elev_at_centres,  dtype=float)
    data_out        = np.asarray(data,              dtype=float).copy()

    nz = z_nodes.shape[0] - 1
    nx = data_out.shape[1]

    z_centres = (z_nodes[:-1] + z_nodes[1:]) / 2.0   # (nz,)
    max_elev  = np.max(elev_at_centres) if len(elev_at_centres) else 0.0

    for j in range(nx):
        surf = elev_at_centres[j] if j < len(elev_at_centres) else 0.0
        for k in range(nz):
            cell_elev = surf - z_centres[k] * exaggeration
            if cell_elev > max_elev:
                data_out[k, j] = np.nan
    return data_out


def station_surface_z(
    chainage_km: np.ndarray,
    elev_km: np.ndarray,
    station_x_km: np.ndarray,
    exaggeration: float = 1.0,
) -> np.ndarray:
    """Return the terrain-draped z-coordinate for station marker positions.

    In the draped coordinate frame the station sits at the surface
    elevation, not at z = 0.

    Parameters
    ----------
    chainage_km : (n,) station chainage values
    elev_km : (n,) elevation at each station (km)
    station_x_km : (m,) x positions where markers should be drawn
    exaggeration : float

    Returns
    -------
    numpy.ndarray, shape (m,)
        z-coordinates (km) at which to place station markers.
    """
    elev_interp = interp_elev(chainage_km, elev_km, station_x_km)
    # In the draped frame surface is at elev; z_nodes[0] = 0 maps to elev.
    # The marker sits exactly at the surface, so its draped z = elev.
    # With exaggeration the datum shifts but the surface stays at elev.
    return elev_interp * exaggeration if exaggeration != 1.0 else elev_interp

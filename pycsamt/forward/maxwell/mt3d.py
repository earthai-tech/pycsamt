# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Research-only, small-grid 3-D MT finite-difference adapter.

**Status: research-only, not for production use.** Per
``docs/source/development/adr/AI-INVERSION-M6-3D-ADR.md`` (the M6
architecture-decision record),
an in-house 3-D Maxwell solver failed the feasibility gate for
production training: direct sparse solves do not scale to realistic
3-D mesh sizes (empirically confirmed there; extrapolates to tens of
minutes per single solve at production resolution). This module exists
for small-grid research and unit validation only. Production 3-D work
should use :mod:`pycsamt.forward.maxwell.modem3d` (a trusted external
backend) once available.

Physics
-------
:class:`MT3DAdapter` discretizes the frequency-domain equation

.. math::

    \\nabla \\times \\nabla \\times E + i \\omega \\mu_0 \\sigma E = 0

on a Cartesian Yee (staggered edge/face) grid, with E on cell edges
and H on cell faces. Cell widths may be non-uniform per axis (a
padded/graded tensor mesh, like
:mod:`pycsamt.forward.maxwell.mt2d` already uses) — see "Non-uniform
mesh support" below. The discrete curl operators (:func:`_curl_e2h`,
:func:`_curl_h2e`) are built and verified independently of the
physics: applying :func:`_curl_e2h` to a uniform field gives exactly
zero, to a linear field gives the exact analytic curl, and
``face_divergence @ curl_e2h`` is the exact zero matrix (the
topological identity ``div(curl(E)) = 0``), on both uniform and
non-uniform grids — see ``pycsamt/forward/tests/test_maxwell_mt3d.py``.
The overall sign convention (``H = -curl(E) / (i omega mu0)``,
assembled as
``(curl_h2e @ curl_e2h + i*omega*mu0*diag(sigma_edge)) @ E = 0``) was
anchored empirically against the analytic half-space limit, not
assumed from a textbook derivation, because this codebase's
depth-increases-downward coordinate convention does not necessarily
share a textbook right-handed frame's sign. Two independent horizontal
polarizations are solved per frequency (boundary ``Ex`` driven and
boundary ``Ey`` driven); receivers combine both to recover the full
impedance tensor.

Non-uniform mesh support
-------------------------
``supports_nonuniform_mesh=True``. Earlier versions of this module
supported uniform cell spacing only, and the *documented* consequence
was that :func:`~pycsamt.forward.maxwell.benchmarks.layered_earth_benchmark`
failed by 30-45% even after mesh refinement. Investigating that
failure (see git history / the AI-inversion project memory for the
session this was diagnosed in) found it was **not** a physics defect
in the boundary-condition approximation itself: the per-column
decay in :func:`_column_decay` is exact for a laterally uniform
half-space regardless of domain size (that is why the half-space
benchmark always passed), and for a genuinely layered earth its error
stays negligible wherever the field has already decayed close to zero
by the domain edge — the same "boundary far enough away" argument
:mod:`pycsamt.forward.maxwell.mt2d` documents for its own boundary
treatment. The real problem was that a *uniform*-only mesh cannot
reach several skin depths of lateral/vertical extent **and** resolve
a few-hundred-metre layer interface within the ``maximum_cells``
budget at the same time, because 3-D cell count scales as the cube of
resolution — refining resolution while keeping the same small domain
does not help, which is exactly the "refining the mesh does not
reduce them" symptom that was previously (incorrectly) attributed to
the boundary approximation. Generalizing the curl operators
(:func:`_curl_e2h`, :func:`_curl_h2e`) and edge-conductivity averaging
(:func:`_sigma_on_edges`) to per-axis, non-uniform cell widths lets a
padded mesh (fine cells near the receiver/structure, geometrically
growing cells outward) reach the same physical extent at a fraction
of the cell cost, matching what :mod:`pycsamt.forward.maxwell.mt2d`
already does. :func:`_curl_h2e` in particular must divide by the
*dual*-grid spacing (the average of the two neighbouring primal cell
widths) rather than either cell's own width — the two coincide on a
uniform mesh, which is why the uniform-only implementation could get
away without this distinction.

Deliberate scope reductions (all enforced, not just documented)
-----------------------------------------------------------------

* **Small grids only.** ``maximum_cells`` (default 6,000, overridable)
  rejects larger problems outright rather than silently taking an
  impractical amount of time, operationalizing the ADR's feasibility
  finding directly in code. A padded, non-uniform mesh (see above)
  substantially relaxes what this budget can achieve, but the ADR's
  underlying conclusion — that a direct solve does not scale to
  realistic production 3-D mesh sizes — is unaffected: this module
  remains research/small-grid only.
* **Direct sparse solve only** (:func:`scipy.sparse.linalg.spsolve`);
  no iterative solver or preconditioner, per the ADR's conclusion that
  a real one is a research problem in itself.
* **Surface receivers only**, within the mesh's horizontal extent,
  checked in :meth:`MT3DAdapter.assess`.
* **Boundary conditions are a per-column exponential approximation**
  (:func:`_column_decay`, the same level of rigor as
  :mod:`pycsamt.forward.em2d`'s own ``_ey_1d_profile``): each boundary
  edge uses the exponential decay implied by its nearest cell column's
  *single local layer*, not a full multi-layer recursion. As described
  above, this is exact for a half-space and negligible-error elsewhere
  once the domain is wide/deep enough, not a source of the previously
  measured layered-earth bias. ``Ez`` is fixed at zero on every
  boundary (plane-wave incidence has no driven vertical field far from
  the domain interior).
* Isotropic conductivity, vacuum permeability, ``exp(+iwt)`` only, no
  inactive-cell/topography support — same restrictions as
  :class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter`.

Measured accuracy
------------------
On a padded, non-uniform 16x16x16 (4,096-cell) research-scale grid —
fine (150 m / 120 m) core cells near the receiver and shallow
interfaces, geometrically padded out to ~30 km laterally and ~8 km
vertically, comfortably beyond the deepest layer's skin depth at the
lowest benchmark frequency — this adapter passes the default
:class:`~pycsamt.forward.maxwell.benchmarks.BenchmarkThresholds` for
**both** :func:`~pycsamt.forward.maxwell.benchmarks.half_space_benchmark`
(~1.9% normalized RMS, ~1.3 degree phase error) **and**
:func:`~pycsamt.forward.maxwell.benchmarks.layered_earth_benchmark`
(~3.5% normalized RMS, ~4.5% amplitude error, <1 degree phase error);
see ``pycsamt/forward/tests/test_maxwell_mt3d.py``.
:data:`_VERIFIED_BENCHMARKS` now includes both.
"""

from __future__ import annotations

import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .adapters import AdapterPolicy, BaseMaxwellAdapter
from .backends import BackendCapabilities, CompatibilityReport
from .contracts import ForwardResult, MaxwellProblem, SolverDiagnostics

__all__ = ["MT3DAdapter"]

_MU0 = 4.0e-7 * np.pi
_SURFACE_TOLERANCE_M = 1e-6
_PERMEABILITY_RTOL = 1e-9
_DEFAULT_MAX_CELLS = 6_000

# Both benchmarks genuinely pass pycsamt/forward/tests/test_maxwell_mt3d.py's
# real benchmark run, on the padded/non-uniform calibrated mesh there
# (see "Measured accuracy" above) -- not silently rounded up.
_VERIFIED_BENCHMARKS: tuple[str, ...] = ("half-space", "layered-earth")


def _edge_shapes(nx: int, ny: int, nz: int) -> dict[str, tuple[int, ...]]:
    return {
        "ex": (nx, ny + 1, nz + 1),
        "ey": (nx + 1, ny, nz + 1),
        "ez": (nx + 1, ny + 1, nz),
        "hx": (nx + 1, ny, nz),
        "hy": (nx, ny + 1, nz),
        "hz": (nx, ny, nz + 1),
    }


def _dual_widths(h: np.ndarray) -> np.ndarray:
    """Return node-centred dual-grid spacings, length ``len(h) + 1``.

    Entry ``i`` (for ``1 <= i <= len(h) - 1``) is the distance between
    the centres of primal cells ``i - 1`` and ``i``, i.e.
    ``0.5 * (h[i - 1] + h[i])``. Entries 0 and ``len(h)`` are unused
    (those nodes are mesh boundaries, handled separately) and left as
    zero so an accidental use fails loudly (division by zero) rather
    than silently.
    """
    n = len(h)
    dual = np.zeros(n + 1)
    dual[1:n] = 0.5 * (h[:-1] + h[1:])
    return dual


def _curl_e2h(
    nx: int,
    ny: int,
    nz: int,
    hx: np.ndarray,
    hy: np.ndarray,
    hz: np.ndarray,
) -> sparse.csr_matrix:
    """Discrete curl E -> H on a (possibly non-uniform) tensor grid.

    ``hx``/``hy``/``hz`` are the primal cell widths along each axis
    (length ``nx``/``ny``/``nz``). Each finite difference here spans
    exactly one primal cell, so it is divided by that cell's own
    width at the relevant index -- no averaging needed.
    """
    s = _edge_shapes(nx, ny, nz)
    n_ex, n_ey, n_ez = (int(np.prod(s[k])) for k in ("ex", "ey", "ez"))
    n_hx, n_hy, n_hz = (int(np.prod(s[k])) for k in ("hx", "hy", "hz"))
    off_ey, off_ez = n_ex, n_ex + n_ey
    off_hy, off_hz = n_hx, n_hx + n_hy

    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []

    def add(r: int, c: int, v: float) -> None:
        rows.append(r)
        cols.append(c)
        vals.append(v)

    def ex_idx(i: int, j: int, k: int) -> int:
        return int(np.ravel_multi_index((i, j, k), s["ex"]))

    def ey_idx(i: int, j: int, k: int) -> int:
        return int(np.ravel_multi_index((i, j, k), s["ey"]))

    def ez_idx(i: int, j: int, k: int) -> int:
        return int(np.ravel_multi_index((i, j, k), s["ez"]))

    for i in range(nx + 1):
        for j in range(ny):
            for k in range(nz):
                r = int(np.ravel_multi_index((i, j, k), s["hx"]))
                add(r, off_ez + ez_idx(i, j + 1, k), 1.0 / hy[j])
                add(r, off_ez + ez_idx(i, j, k), -1.0 / hy[j])
                add(r, off_ey + ey_idx(i, j, k + 1), -1.0 / hz[k])
                add(r, off_ey + ey_idx(i, j, k), 1.0 / hz[k])

    for i in range(nx):
        for j in range(ny + 1):
            for k in range(nz):
                r = off_hy + int(np.ravel_multi_index((i, j, k), s["hy"]))
                add(r, ex_idx(i, j, k + 1), 1.0 / hz[k])
                add(r, ex_idx(i, j, k), -1.0 / hz[k])
                add(r, off_ez + ez_idx(i + 1, j, k), -1.0 / hx[i])
                add(r, off_ez + ez_idx(i, j, k), 1.0 / hx[i])

    for i in range(nx):
        for j in range(ny):
            for k in range(nz + 1):
                r = off_hz + int(np.ravel_multi_index((i, j, k), s["hz"]))
                add(r, off_ey + ey_idx(i + 1, j, k), 1.0 / hx[i])
                add(r, off_ey + ey_idx(i, j, k), -1.0 / hx[i])
                add(r, ex_idx(i, j + 1, k), -1.0 / hy[j])
                add(r, ex_idx(i, j, k), 1.0 / hy[j])

    n_e = n_ex + n_ey + n_ez
    n_h = n_hx + n_hy + n_hz
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n_h, n_e))


def _curl_h2e(
    nx: int,
    ny: int,
    nz: int,
    hx: np.ndarray,
    hy: np.ndarray,
    hz: np.ndarray,
) -> sparse.csr_matrix:
    """Discrete curl H -> E on a (possibly non-uniform) tensor grid.

    Each finite difference here connects two H values that straddle
    an E edge from *opposite* primal cells (e.g. ``Hz`` at cell
    centres ``j - 1`` and ``j``), so the correct denominator is the
    *dual*-grid spacing between those two cell centres --
    ``0.5 * (h[j - 1] + h[j])`` -- not either cell's own width. Using
    the primal width here (as a uniform-only implementation legitimately
    can, since dual and primal spacing coincide when all cells are the
    same size) silently double-counts curvature on a padded mesh.
    """
    s = _edge_shapes(nx, ny, nz)
    n_ex, n_ey, n_ez = (int(np.prod(s[k])) for k in ("ex", "ey", "ez"))
    n_hx, n_hy, n_hz = (int(np.prod(s[k])) for k in ("hx", "hy", "hz"))
    off_ey, off_ez = n_ex, n_ex + n_ey
    off_hy, off_hz = n_hx, n_hx + n_hy

    dual_x, dual_y, dual_z = (
        _dual_widths(hx),
        _dual_widths(hy),
        _dual_widths(hz),
    )

    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []

    def add(r: int, c: int, v: float) -> None:
        rows.append(r)
        cols.append(c)
        vals.append(v)

    def hx_idx(i: int, j: int, k: int) -> int:
        return int(np.ravel_multi_index((i, j, k), s["hx"]))

    def hy_idx(i: int, j: int, k: int) -> int:
        return int(np.ravel_multi_index((i, j, k), s["hy"]))

    def hz_idx(i: int, j: int, k: int) -> int:
        return int(np.ravel_multi_index((i, j, k), s["hz"]))

    for i in range(nx):
        for j in range(1, ny):
            for k in range(1, nz):
                r = int(np.ravel_multi_index((i, j, k), s["ex"]))
                add(r, off_hz + hz_idx(i, j, k), 1.0 / dual_y[j])
                add(r, off_hz + hz_idx(i, j - 1, k), -1.0 / dual_y[j])
                add(r, off_hy + hy_idx(i, j, k), -1.0 / dual_z[k])
                add(r, off_hy + hy_idx(i, j, k - 1), 1.0 / dual_z[k])

    for i in range(1, nx):
        for j in range(ny):
            for k in range(1, nz):
                r = off_ey + int(np.ravel_multi_index((i, j, k), s["ey"]))
                add(r, hx_idx(i, j, k), 1.0 / dual_z[k])
                add(r, hx_idx(i, j, k - 1), -1.0 / dual_z[k])
                add(r, off_hz + hz_idx(i, j, k), -1.0 / dual_x[i])
                add(r, off_hz + hz_idx(i - 1, j, k), 1.0 / dual_x[i])

    for i in range(1, nx):
        for j in range(1, ny):
            for k in range(nz):
                r = off_ez + int(np.ravel_multi_index((i, j, k), s["ez"]))
                add(r, off_hy + hy_idx(i, j, k), 1.0 / dual_x[i])
                add(r, off_hy + hy_idx(i - 1, j, k), -1.0 / dual_x[i])
                add(r, hx_idx(i, j, k), -1.0 / dual_y[j])
                add(r, hx_idx(i, j - 1, k), 1.0 / dual_y[j])

    n_e = n_ex + n_ey + n_ez
    n_h = n_hx + n_hy + n_hz
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n_e, n_h))


def _boundary_masks(
    nx: int, ny: int, nz: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = _edge_shapes(nx, ny, nz)
    ex_bnd = np.zeros(s["ex"], dtype=bool)
    ex_bnd[:, 0, :] = ex_bnd[:, -1, :] = True
    ex_bnd[:, :, 0] = ex_bnd[:, :, -1] = True
    ey_bnd = np.zeros(s["ey"], dtype=bool)
    ey_bnd[0, :, :] = ey_bnd[-1, :, :] = True
    ey_bnd[:, :, 0] = ey_bnd[:, :, -1] = True
    ez_bnd = np.zeros(s["ez"], dtype=bool)
    ez_bnd[0, :, :] = ez_bnd[-1, :, :] = True
    ez_bnd[:, 0, :] = ez_bnd[:, -1, :] = True
    return ex_bnd, ey_bnd, ez_bnd


def _sigma_on_edges(
    sigma_xyz: np.ndarray,
    hx: np.ndarray,
    hy: np.ndarray,
    hz: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average cell conductivity onto edges, weighted by transverse area.

    Each interior edge is surrounded by 4 cells in the plane
    perpendicular to it; a plain 1/4-1/4-1/4-1/4 average (correct only
    when all 4 cells have equal cross-section) silently mis-weights a
    padded/non-uniform mesh in favour of whichever neighbour happens to
    come first in memory. Weighting by each cell's transverse
    cross-sectional area is the same convention already used by
    :func:`pycsamt.forward.em2d._assemble_te`'s nodal average.
    """
    sx = np.zeros((nx, ny + 1, nz + 1))
    wy_l, wy_r = hy[: ny - 1, None], hy[1:ny, None]
    wz_l, wz_r = hz[None, : nz - 1], hz[None, 1:nz]
    w_tl, w_tr = wy_l * wz_l, wy_r * wz_l
    w_bl, w_br = wy_l * wz_r, wy_r * wz_r
    w_sum = w_tl + w_tr + w_bl + w_br
    sx[:, 1:ny, 1:nz] = (
        w_tl * sigma_xyz[:, 0 : ny - 1, 0 : nz - 1]
        + w_tr * sigma_xyz[:, 1:ny, 0 : nz - 1]
        + w_bl * sigma_xyz[:, 0 : ny - 1, 1:nz]
        + w_br * sigma_xyz[:, 1:ny, 1:nz]
    ) / w_sum

    sy = np.zeros((nx + 1, ny, nz + 1))
    wx_l, wx_r = hx[: nx - 1, None, None], hx[1:nx, None, None]
    wz_l, wz_r = hz[None, None, : nz - 1], hz[None, None, 1:nz]
    w_tl, w_tr = wx_l * wz_l, wx_r * wz_l
    w_bl, w_br = wx_l * wz_r, wx_r * wz_r
    w_sum = w_tl + w_tr + w_bl + w_br
    sy[1:nx, :, 1:nz] = (
        w_tl * sigma_xyz[0 : nx - 1, :, 0 : nz - 1]
        + w_tr * sigma_xyz[1:nx, :, 0 : nz - 1]
        + w_bl * sigma_xyz[0 : nx - 1, :, 1:nz]
        + w_br * sigma_xyz[1:nx, :, 1:nz]
    ) / w_sum

    sz = np.zeros((nx + 1, ny + 1, nz))
    wx_l, wx_r = hx[: nx - 1, None, None], hx[1:nx, None, None]
    wy_l, wy_r = hy[None, : ny - 1, None], hy[None, 1:ny, None]
    w_tl, w_tr = wx_l * wy_l, wx_r * wy_l
    w_bl, w_br = wx_l * wy_r, wx_r * wy_r
    w_sum = w_tl + w_tr + w_bl + w_br
    sz[1:nx, 1:ny, :] = (
        w_tl * sigma_xyz[0 : nx - 1, 0 : ny - 1, :]
        + w_tr * sigma_xyz[1:nx, 0 : ny - 1, :]
        + w_bl * sigma_xyz[0 : nx - 1, 1:ny, :]
        + w_br * sigma_xyz[1:nx, 1:ny, :]
    ) / w_sum
    return sx, sy, sz


def _column_decay(
    sigma_profile: np.ndarray,
    dz_profile: np.ndarray,
    z_query: np.ndarray,
    omega: float,
) -> np.ndarray:
    wavenumber = np.sqrt(1j * omega * _MU0 * sigma_profile)
    top = np.concatenate([[0.0], np.cumsum(dz_profile)])
    layer = np.clip(
        np.searchsorted(top, z_query, side="right") - 1,
        0,
        len(sigma_profile) - 1,
    )
    return np.exp(-wavenumber[layer] * z_query)


def _nearest_column(
    sigma_xyz: np.ndarray, i: int, j: int, nx: int, ny: int
) -> np.ndarray:
    return sigma_xyz[min(max(i, 0), nx - 1), min(max(j, 0), ny - 1), :]


def _boundary_fields(
    driven: str,
    nx: int,
    ny: int,
    nz: int,
    hz: np.ndarray,
    sigma_xyz: np.ndarray,
    omega: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dirichlet boundary values for one polarization ("x" or "y")."""
    s = _edge_shapes(nx, ny, nz)
    dz_profile = hz[:-1] if nz > 1 else np.zeros(0)
    z_ex = np.concatenate([[0.0], np.cumsum(hz)])

    bc_ex = np.zeros(s["ex"], dtype=complex)
    bc_ey = np.zeros(s["ey"], dtype=complex)
    bc_ez = np.zeros(s["ez"], dtype=complex)  # Ez is always 0 on boundary

    target = bc_ex if driven == "x" else bc_ey
    n_i, n_j = target.shape[0], target.shape[1]
    for i in range(n_i):
        for j in range(n_j):
            column = _nearest_column(sigma_xyz, i, j, nx, ny)
            target[i, j, :] = _column_decay(column, dz_profile, z_ex, omega)
    return bc_ex, bc_ey, bc_ez


def _assemble(
    nx: int,
    ny: int,
    nz: int,
    hx: np.ndarray,
    hy: np.ndarray,
    hz: np.ndarray,
    sigma_xyz: np.ndarray,
    omega: float,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, np.ndarray]:
    c1 = _curl_e2h(nx, ny, nz, hx, hy, hz)
    c2 = _curl_h2e(nx, ny, nz, hx, hy, hz)
    sx, sy, sz = _sigma_on_edges(sigma_xyz, hx, hy, hz, nx, ny, nz)
    sigma_edge = np.concatenate([sx.ravel(), sy.ravel(), sz.ravel()])

    a = (c2 @ c1) + sparse.diags(1j * omega * _MU0 * sigma_edge)

    ex_bnd, ey_bnd, ez_bnd = _boundary_masks(nx, ny, nz)
    boundary = np.concatenate([ex_bnd.ravel(), ey_bnd.ravel(), ez_bnd.ravel()])
    idx = np.flatnonzero(boundary)
    a = a.tocsr()
    for r in idx:
        a.data[a.indptr[r] : a.indptr[r + 1]] = 0.0
    a = a.tolil()
    for r in idx:
        a[r, r] = 1.0
    a = a.tocsr()
    a.eliminate_zeros()
    return a, c1, boundary


def _solve_polarization(
    a: sparse.csr_matrix,
    c1: sparse.csr_matrix,
    boundary: np.ndarray,
    bc_ex: np.ndarray,
    bc_ey: np.ndarray,
    bc_ez: np.ndarray,
    omega: float,
) -> tuple[np.ndarray, float]:
    b = np.zeros(a.shape[0], dtype=complex)
    bc_flat = np.concatenate([bc_ex.ravel(), bc_ey.ravel(), bc_ez.ravel()])
    idx = np.flatnonzero(boundary)
    b[idx] = bc_flat[idx]
    x = spsolve(a.tocsc(), b)
    denom = np.linalg.norm(b)
    residual = float(
        np.linalg.norm(a @ x - b) / (denom if denom > 0 else 1e-300)
    )
    h = -(c1 @ x) / (1j * omega * _MU0)
    return np.concatenate([x, h]), residual


def _bilinear(
    values: np.ndarray,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    x_query: float,
    y_query: float,
) -> complex:
    ix = np.clip(
        np.searchsorted(x_centers, x_query) - 1, 0, len(x_centers) - 2
    )
    iy = np.clip(
        np.searchsorted(y_centers, y_query) - 1, 0, len(y_centers) - 2
    )
    x0, x1 = x_centers[ix], x_centers[ix + 1]
    y0, y1 = y_centers[iy], y_centers[iy + 1]
    tx = 0.0 if x1 == x0 else (x_query - x0) / (x1 - x0)
    ty = 0.0 if y1 == y0 else (y_query - y0) / (y1 - y0)
    v00, v01 = values[ix, iy], values[ix, iy + 1]
    v10, v11 = values[ix + 1, iy], values[ix + 1, iy + 1]
    return (
        v00 * (1 - tx) * (1 - ty)
        + v10 * tx * (1 - ty)
        + v01 * (1 - tx) * ty
        + v11 * tx * ty
    )


_COMPONENT_INDEX = {"zxx": (0, 0), "zxy": (0, 1), "zyx": (1, 0), "zyy": (1, 1)}


class MT3DAdapter(BaseMaxwellAdapter):
    """Research-only 3-D MT adapter (see module docstring for scope).

    Parameters
    ----------
    version : str, default="1.0-research"
        Adapter version reported in every
        :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.
    policy : AdapterPolicy or None, optional
        Solver-independent result acceptance policy.
    max_cells : int, default=6000
        Safety ceiling on total mesh cells
        (:attr:`~pycsamt.forward.maxwell.backends.BackendCapabilities.maximum_cells`).
        A direct sparse solve becomes impractically slow well before
        typical production 3-D mesh sizes; see
        ``docs/source/development/adr/AI-INVERSION-M6-3D-ADR.md``. Raise this
        only if you have
        confirmed the resulting solve time is acceptable for your use.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.maxwell import (
    ...     MaxwellMesh,
    ...     MaxwellProblem,
    ...     ReceiverSet,
    ... )
    >>> mesh = MaxwellMesh(
    ...     np.linspace(0, 4000, 9),
    ...     np.linspace(0, 3000, 11),
    ...     np.linspace(0, 4000, 9),
    ... )
    >>> problem = MaxwellProblem(
    ...     mesh,
    ...     np.full(mesh.shape, 1.0 / 100.0),
    ...     [1.0],
    ...     ReceiverSet([[2000.0, 2000.0, 0.0]], ["S00"]),
    ...     ("zxy", "zyx"),
    ... )
    >>> result = MT3DAdapter().solve(problem)
    >>> result.shape
    (1, 1, 2)
    """

    def __init__(
        self,
        *,
        version: str = "1.0-research",
        policy: AdapterPolicy | None = None,
        max_cells: int = _DEFAULT_MAX_CELLS,
    ) -> None:
        capabilities = BackendCapabilities(
            name="mt3d",
            version=version,
            dimensions=(3,),
            components=("zxx", "zxy", "zyx", "zyy"),
            time_conventions=("exp(+iwt)",),
            supports_nonuniform_mesh=True,
            supports_inactive_cells=False,
            supports_topography=False,
            supports_anisotropy=False,
            maximum_cells=max_cells,
            verified_benchmarks=_VERIFIED_BENCHMARKS,
        )
        super().__init__(capabilities, policy)

    def assess(self, problem: MaxwellProblem) -> CompatibilityReport:
        """Assess a problem, adding this solver's research-only checks.

        Parameters
        ----------
        problem : MaxwellProblem
            Candidate simulation problem.

        Returns
        -------
        CompatibilityReport
            The generic capability report from
            :meth:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter.assess`,
            extended with surface-receiver, horizontal-bounds, and
            vacuum-permeability checks.

        Examples
        --------
        See :class:`MT3DAdapter` for a complete solve example.
        """
        base = super().assess(problem)
        errors = list(base.errors)
        # The coordinate/permeability checks below assume a 3-D
        # receiver/mesh shape; a wrong-dimension problem is already
        # rejected by the base capability report, so skip them rather
        # than index into a 2-D receiver array.
        if problem.mesh.dimension == 3:
            coordinates = problem.receivers.coordinates_m
            depth = coordinates[:, 2]
            if np.any(np.abs(depth) > _SURFACE_TOLERANCE_M):
                errors.append(
                    "mt3d only evaluates receivers at the surface (z=0)"
                )
            x, y = coordinates[:, 0], coordinates[:, 1]
            x_lo, x_hi = problem.mesh.x_edges_m[0], problem.mesh.x_edges_m[-1]
            y_lo, y_hi = problem.mesh.y_edges_m[0], problem.mesh.y_edges_m[-1]
            if np.any((x < x_lo) | (x > x_hi) | (y < y_lo) | (y > y_hi)):
                errors.append(
                    "receiver x/y coordinates must lie within the mesh "
                    f"[{x_lo:g}, {x_hi:g}] x [{y_lo:g}, {y_hi:g}] m"
                )
        if not np.isclose(
            problem.magnetic_permeability_h_m,
            _MU0,
            rtol=_PERMEABILITY_RTOL,
        ):
            errors.append(
                f"mt3d assumes vacuum magnetic permeability ({_MU0:.6g} H/m)"
            )
        if len(errors) == len(base.errors):
            return base
        return CompatibilityReport(
            self.capabilities.name, False, tuple(errors), base.warnings
        )

    def _solve_backend(self, problem: MaxwellProblem) -> ForwardResult:
        mesh = problem.mesh
        nz, ny, nx = mesh.shape
        widths = mesh.cell_widths_m
        hx, hy, hz = widths["x"], widths["y"], widths["z"]
        sigma_xyz = np.transpose(problem.conductivity_s_m, (2, 1, 0))
        x_centers = mesh.cell_centres_m["x"]
        y_centers = mesh.cell_centres_m["y"]
        x_nodes = mesh.x_edges_m
        y_nodes = mesh.y_edges_m

        n_freq = len(problem.frequencies_hz)
        n_station = problem.receivers.count
        n_component = len(problem.components)
        impedance = np.empty((n_station, n_freq, n_component), dtype=complex)
        converged = np.empty((n_freq, 2), dtype=bool)
        residual = np.empty((n_freq, 2))
        start = time.monotonic()

        n_ex = nx * (ny + 1) * (nz + 1)
        n_ey = (nx + 1) * ny * (nz + 1)
        n_hx = (nx + 1) * ny * nz
        n_hy = nx * (ny + 1) * nz

        for fi, frequency in enumerate(problem.frequencies_hz):
            omega = 2.0 * np.pi * float(frequency)
            a, c1, boundary = _assemble(
                nx, ny, nz, hx, hy, hz, sigma_xyz, omega
            )
            hx_by_pol = []
            hy_by_pol = []
            for pi, driven in enumerate(("x", "y")):
                bc_ex, bc_ey, bc_ez = _boundary_fields(
                    driven, nx, ny, nz, hz, sigma_xyz, omega
                )
                solution, res = _solve_polarization(
                    a, c1, boundary, bc_ex, bc_ey, bc_ez, omega
                )
                residual[fi, pi] = res
                converged[fi, pi] = np.isfinite(res)
                h = solution[n_ex + n_ey + ((nx + 1) * (ny + 1) * nz) :]
                hx_field = h[:n_hx].reshape((nx + 1, ny, nz))[:, :, 0]
                hy_field = h[n_hx : n_hx + n_hy].reshape((nx, ny + 1, nz))[
                    :, :, 0
                ]
                hx_by_pol.append(hx_field)
                hy_by_pol.append(hy_field)

            # E at the surface is exactly the prescribed incident field
            # (1 for the driven polarization's own axis, 0 otherwise),
            # since the top boundary is uniform across x, y by
            # construction; no interpolation is needed for E.
            e_mat = np.eye(2, dtype=complex)
            for si in range(n_station):
                x_q, y_q = (
                    problem.receivers.coordinates_m[si, 0],
                    problem.receivers.coordinates_m[si, 1],
                )
                h_mat = np.empty((2, 2), dtype=complex)
                for pi in range(2):
                    # Hx lives at (x-node, y-center); Hy at (x-center,
                    # y-node) -- the two Yee dual-grid H locations.
                    h_mat[0, pi] = _bilinear(
                        hx_by_pol[pi], x_nodes, y_centers, x_q, y_q
                    )
                    h_mat[1, pi] = _bilinear(
                        hy_by_pol[pi], x_centers, y_nodes, x_q, y_q
                    )
                z_full = e_mat @ np.linalg.inv(h_mat)
                for ci, name in enumerate(problem.components):
                    row, col = _COMPONENT_INDEX[name]
                    impedance[si, fi, ci] = z_full[row, col]

        runtime_s = time.monotonic() - start
        diagnostics = SolverDiagnostics(
            converged,
            np.zeros(residual.shape, dtype=int),
            np.where(np.isfinite(residual), residual, 0.0),
            runtime_s,
        )
        return ForwardResult(
            problem.problem_hash,
            problem.frequencies_hz,
            problem.receivers.names,
            problem.components,
            impedance,
            None,
            self.capabilities.name,
            self.capabilities.version,
            diagnostics,
        )

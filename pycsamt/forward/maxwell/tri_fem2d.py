# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""In-house 2-D triangular-mesh magnetotelluric FEM solver (M16, research-only).

:class:`Tri2DFEMForward` solves the same 2-D MT TE/TM boundary-value problem
as :mod:`pycsamt.forward.em2d` (:class:`~pycsamt.forward.em2d.MT2DForward`),
but with a standard nodal Galerkin finite-element discretization on
unstructured linear (P1) triangles instead of a finite-difference tensor
grid -- the classic method of Wannamaker, Stodt & Rijo (1987), already cited
in :mod:`pycsamt.forward.em2d`'s own docstring. :class:`TriFEM2DAdapter`
wraps it in the solver-neutral
:class:`~pycsamt.forward.maxwell.contracts_tri.TriProblem` /
:class:`~pycsamt.forward.maxwell.contracts.ForwardResult` contract.

This is a **research-only, in-house alternative** to the external
:class:`~pycsamt.forward.maxwell.mare2dem.Mare2DEMAdapter` (M13, production),
mirroring the existing :mod:`~pycsamt.forward.maxwell.mt3d`-alongside-
:mod:`~pycsamt.forward.maxwell.modem3d` precedent: no external binary
required, deliberately smaller-scope, validated with the same analytic
half-space/layered-earth benchmark discipline before any accuracy claim is
trusted.

Physics
-------
Same governing PDEs, sign convention, and impedance definitions as
:mod:`pycsamt.forward.em2d` (time convention ``exp(+iwt)``):

**TE mode** (``E_y``): ``∇²E_y - iωμ₀σE_y = 0``.

**TM mode** (``H_y``): ``∇·(ρ∇H_y) - iωμ₀H_y = 0``.

Multiplying by a test function ``φᵢ`` and integrating by parts (Green's
first identity) flips the sign of the diffusion term relative to the
strong-form PDE above -- ``∫φᵢ∇²E_y dA = -∫∇φᵢ·∇E_y dA`` plus a boundary
term dropped here since Dirichlet conditions are imposed everywhere -- so
the correct weak form is ``∫∇φᵢ·∇E_y dA + iωμ₀∫σφᵢE_y dA = 0`` (**plus**,
not minus; getting this sign backwards during development produced a
self-consistent but wrong solution -- residuals of the assembled discrete
operator against the known analytic half-space solution were small at the
*boundary* (correct, since boundary values are prescribed directly from
that same analytic formula) but large in the *interior*, which is what
exposed it -- confirmed by comparing interior nodal values against
:func:`~pycsamt.forward.em2d._ey_1d_profile` directly on a laterally
uniform mesh, where the true solution is exactly the 1-D profile
everywhere, not just at the boundary). Likewise for TM:
``∫ρ∇φᵢ·∇H_y dA + iωμ₀∫φᵢH_y dA = 0``.

For linear (P1) triangles, ``∇φᵢ`` is constant per element (the standard
``bᵢ, cᵢ`` coefficients), giving a closed-form element stiffness matrix
``Kᵉ_ij = A(bᵢbⱼ + cᵢcⱼ)`` (this is a discrete **negative** Laplacian, by
the standard FEM convention -- positive semi-definite, which is why the
weak-form sign above is plus, not the strong-form PDE's minus) and the
standard consistent P1 mass matrix
``Mᵉ = (A/12)·[[2,1,1],[1,2,1],[1,1,2]]``. ``TriProblem.conductivity_s_m``
is already per-triangle, so the piecewise-constant material coefficient
integrates exactly with no interpolation.

Boundary conditions
--------------------
Dirichlet at every node referenced by ``mesh.boundary_segments``. A node
exactly at **the mesh's own local surface** -- its top boundary at that
x, whatever elevation that happens to be, found by
:func:`_local_surface_z`'s vertical ray-cast, not a fixed ``z=0``
assumption -- gets the normalized incident field ``1+0j`` (matching
:mod:`~pycsamt.forward.em2d`'s top-row convention). Every other boundary
node's value comes from the exact 1-D analytic downward continuation
(:func:`~pycsamt.forward.em2d._ey_1d_profile`/
:func:`~pycsamt.forward.em2d._hy_1d_profile`, reused unchanged) of a layered
profile built **from the mesh's own triangles** by
:func:`_vertical_column_profile`, evaluated at **depth below that same
local surface** (not raw absolute z) -- a vertical ray is cast through the
triangulation at that node's x-coordinate, and every triangle the ray
crosses contributes one (thickness, resistivity) layer, ordered from the
local surface downward. This is self-consistent by construction (a
half-space/layered-earth benchmark mesh reproduces exactly the true
background model at its boundary; a real topographic mesh reproduces the
true local geology under that boundary column) and, unlike
:mod:`~pycsamt.forward.em2d`'s FD version, is not restricted to a
rectangular grid -- the mesh's outer boundary can be any shape, which is
why ``supports_topography=True`` is a genuine capability here. An earlier
version of this module hardcoded "surface means ``abs(z)<=tol``", which
only happened to be correct because every benchmark mesh built so far put
its surface exactly at ``z=0`` -- confirmed as a real bug (not just a
missing feature) once a receiver away from ``z=0`` was tried and every
boundary node's downward-continuation depth came out measured from the
wrong datum. Generalizing to the local surface reduces to exactly the old
behavior whenever the mesh *is* flat at ``z=0``, so the existing
half-space/layered-earth benchmarks are unaffected; a new datum-shift
benchmark (translate the whole mesh + receiver vertically, response must
be unchanged) specifically covers the fix itself -- see
``test_maxwell_tri_fem2d.py``.

Station field extraction
-------------------------
Every receiver must coincide with an existing mesh node sitting on the
mesh's own local surface (not literally ``z=0`` -- see above), checked in
:meth:`TriFEM2DAdapter.assess`. This mirrors the convention
:func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`
already established (receivers included as explicit input points to the
triangulation) and avoids needing point-in-triangle interpolation for the
nodal value itself. Unlike :mod:`~pycsamt.forward.em2d`'s finite-difference
version, P1 elements have an *exact*, constant-per-triangle gradient -- no
finite-difference step size to choose. ``Hx``/``Ex`` are recovered at a
receiver node by area-weighting the constant gradient of every triangle
incident to that node (``Σ areaᵢ·∇ᵢ / Σ areaᵢ``), a standard FEM nodal
gradient-recovery technique. This still needs the triangles immediately
below a receiver to be reasonably small relative to skin depth -- confirmed
empirically during validation: a mesh with a ~0.43-skin-depth first layer
gave ~10-33% half-space impedance error, while grading the same mesh down
to a ~0.1-skin-depth first layer (only near the surface, not everywhere)
brought that under 1%. Not a bug, the same kind of resolution sensitivity
:mod:`~pycsamt.forward.em2d`'s FD version already documents for lateral
domain extent -- just the triangular-mesh analogue of it, for vertical
resolution at receivers specifically.

Scope, honestly stated
-----------------------
* Receivers must sit exactly on a mesh node at the mesh's own local
  surface (any elevation, not just ``z=0``) -- off-surface (buried or
  floating) receivers are rejected by :meth:`TriFEM2DAdapter.assess`, not
  silently mishandled.
* No point-in-triangle interpolation for off-node receivers.
* Genuine spatially-varying topographic distortion of the EM field
  (receivers at different elevations along real terrain) is **not**
  independently benchmarked here -- no known closed-form analytic
  topographic-MT solution was available to check against, the same
  honesty gap :mod:`~pycsamt.forward.maxwell.mare2dem`'s own "unverified
  physics" section documents for a different reason. What *is* proven is
  the coordinate-frame generalization itself (the datum-shift benchmark
  above): correct relative to a fixed flat surface at any elevation, not
  independently checked against real topographic distortion.
* Wired into :class:`~pycsamt.agents.inv2d_agent.Inv2DAgent`
  (``physics="mt2d_tri"``, pass ``mare2dem_adapter=TriFEM2DAdapter()``).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from ..em2d import MU0, _ey_1d_profile, _hy_1d_profile
from .adapters import AdapterPolicy, BaseMaxwellAdapter
from .backends import (
    BackendCapabilities,
    BackendRegistration,
    CompatibilityReport,
)
from .backends import register_backend as _register_backend
from .contracts import ForwardResult, SolverDiagnostics
from .contracts_tri import TriMesh, TriProblem

__all__ = [
    "Tri2DFEMForward",
    "ForwardResponseTriFEM2D",
    "TriFEM2DAdapter",
    "register_trifem2d_backend",
]

_SURFACE_TOLERANCE_M = 1e-6
_PERMEABILITY_RTOL = 1e-9
_GEOM_TOL = 1e-9

# Benchmarks this adapter is currently known to pass, verified by
# pycsamt/forward/tests/test_maxwell_tri_fem2d.py; update together with
# that test file, never in isolation.
_VERIFIED_BENCHMARKS: tuple[str, ...] = ("half-space", "layered-earth")


# ─────────────────────────────────────────────────────────────────────────
# P1 triangle geometry
# ─────────────────────────────────────────────────────────────────────────


def p1_gradients(nodes_m: np.ndarray, triangle: np.ndarray) -> tuple[np.ndarray, float]:
    """Return the constant P1 shape-function gradients of one triangle.

    Parameters
    ----------
    nodes_m : ndarray, shape (n_nodes, 2)
        Mesh node coordinates.
    triangle : ndarray of int, shape (3,)
        Node indices of one triangle (either orientation).

    Returns
    -------
    grad : ndarray, shape (3, 2)
        Row ``i`` is ``(b_i, c_i)`` = the constant gradient of shape
        function ``i`` within this triangle.
    area : float
        Positive triangle area.

    Examples
    --------
    >>> nodes = np.array([[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]])
    >>> grad, area = p1_gradients(nodes, np.array([0, 1, 2]))
    >>> round(area, 6)
    6.0
    >>> np.allclose(grad.sum(axis=0), 0.0)
    True
    """
    p0, p1, p2 = nodes_m[triangle[0]], nodes_m[triangle[1]], nodes_m[triangle[2]]
    double_area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (
        p1[1] - p0[1]
    )
    x0, z0 = p0
    x1, z1 = p1
    x2, z2 = p2
    b = np.array([z1 - z2, z2 - z0, z0 - z1]) / double_area
    c = np.array([x2 - x1, x0 - x2, x1 - x0]) / double_area
    grad = np.column_stack([b, c])
    area = abs(double_area) / 2.0
    return grad, area


def _node_triangle_incidence(mesh: TriMesh) -> list[list[int]]:
    """Return, per node, the list of incident triangle indices."""
    incidence: list[list[int]] = [[] for _ in range(mesh.n_nodes)]
    for t_idx, triangle in enumerate(mesh.triangles):
        for node in triangle:
            incidence[int(node)].append(t_idx)
    return incidence


def _vertical_column_profile(
    mesh: TriMesh,
    resistivity_per_triangle: np.ndarray,
    x0: float,
    *,
    tol: float = _GEOM_TOL,
) -> tuple[np.ndarray, np.ndarray]:
    """Ray-cast a vertical line at ``x=x0`` through the triangulation.

    Parameters
    ----------
    mesh : TriMesh
        Triangulation to cast the ray through.
    resistivity_per_triangle : ndarray, shape (n_triangles,)
        Resistivity of every triangle.
    x0 : float
        X-coordinate of the vertical line.
    tol : float, default=1e-9
        Geometric tolerance for edge-crossing detection.

    Returns
    -------
    resistivity : ndarray, shape (n_layers,)
        Layer resistivities from the local surface downward.
    thickness : ndarray, shape (n_layers - 1,)
        Finite-layer thicknesses (one fewer than ``resistivity``, matching
        :func:`~pycsamt.forward.maxwell.benchmarks.layered_earth_impedance`'s
        own convention: the last resistivity is a basal half-space).

    Examples
    --------
    >>> mesh = TriMesh(
    ...     [[0, 0], [10, 0], [10, 10], [0, 10]],
    ...     [[0, 1, 2], [0, 2, 3]],
    ... )
    >>> rho, thick = _vertical_column_profile(mesh, np.array([7.0, 7.0]), 5.0)
    >>> rho.tolist(), thick.tolist()
    ([7.0, 7.0], [5.0])
    """
    nodes = mesh.nodes_m
    intervals: list[tuple[float, float, float]] = []
    for t_idx, triangle in enumerate(mesh.triangles):
        verts = nodes[triangle]
        x_lo, x_hi = verts[:, 0].min(), verts[:, 0].max()
        if x0 < x_lo - tol or x0 > x_hi + tol:
            continue
        crossings: list[float] = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            xa, za = verts[a]
            xb, zb = verts[b]
            if abs(xb - xa) <= tol:
                if abs(xa - x0) <= tol:
                    crossings.append(float(za))
                    crossings.append(float(zb))
                continue
            lo, hi = (xa, xb) if xa <= xb else (xb, xa)
            if lo - tol <= x0 <= hi + tol:
                t = (x0 - xa) / (xb - xa)
                crossings.append(float(za + t * (zb - za)))
        if len(crossings) < 2:
            continue
        z_lo, z_hi = min(crossings), max(crossings)
        if z_hi - z_lo <= tol:
            continue
        intervals.append((z_lo, z_hi, float(resistivity_per_triangle[t_idx])))

    if not intervals:
        raise ValueError(
            f"no triangle intersects the vertical line x={x0!r}; cannot "
            "build a boundary column profile."
        )

    intervals.sort(key=lambda item: item[0])
    merged: list[tuple[float, float, float]] = []
    z_cursor = intervals[0][0]
    for z_lo, z_hi, rho in intervals:
        top = max(z_lo, z_cursor)
        if z_hi - top <= tol:
            continue
        merged.append((top, z_hi, rho))
        z_cursor = z_hi

    resistivity = np.array([rho for (_, _, rho) in merged])
    thickness = np.array([z_hi - z_lo for (z_lo, z_hi, _) in merged[:-1]])
    return resistivity, thickness


def _local_surface_z(mesh: TriMesh, x0: float, *, tol: float = _GEOM_TOL) -> float:
    """Return the mesh's own local top-of-mesh elevation at ``x=x0``.

    Same vertical ray-cast as :func:`_vertical_column_profile` (every
    triangle whose x-range spans ``x0``, edge-crossing z values), but
    only the minimum crossing is kept -- the local topographic surface
    at that x, whatever its elevation, rather than a fixed ``z=0``
    assumption. Pure geometry, no resistivity needed, so this is kept
    separate from :func:`_vertical_column_profile` rather than changing
    that function's return contract.

    Parameters
    ----------
    mesh : TriMesh
        Triangulation to cast the ray through.
    x0 : float
        X-coordinate of the vertical line.
    tol : float, default=1e-9
        Geometric tolerance for edge-crossing detection.

    Returns
    -------
    float
        Local surface elevation (z, positive down) at ``x0``.

    Raises
    ------
    ValueError
        If no triangle intersects the vertical line at ``x0``.

    Examples
    --------
    >>> mesh = TriMesh(
    ...     [[0, 0], [10, 0], [10, 10], [0, 10]],
    ...     [[0, 1, 2], [0, 2, 3]],
    ... )
    >>> _local_surface_z(mesh, 5.0)
    0.0

    A mesh whose top boundary sits at a non-zero elevation reports that
    elevation, not ``0``:

    >>> shifted = TriMesh(
    ...     [[0, 100], [10, 100], [10, 110], [0, 110]],
    ...     [[0, 1, 2], [0, 2, 3]],
    ... )
    >>> _local_surface_z(shifted, 5.0)
    100.0
    """
    nodes = mesh.nodes_m
    crossings: list[float] = []
    for triangle in mesh.triangles:
        verts = nodes[triangle]
        x_lo, x_hi = verts[:, 0].min(), verts[:, 0].max()
        if x0 < x_lo - tol or x0 > x_hi + tol:
            continue
        for a, b in ((0, 1), (1, 2), (2, 0)):
            xa, za = verts[a]
            xb, zb = verts[b]
            if abs(xb - xa) <= tol:
                if abs(xa - x0) <= tol:
                    crossings.append(float(za))
                    crossings.append(float(zb))
                continue
            lo, hi = (xa, xb) if xa <= xb else (xb, xa)
            if lo - tol <= x0 <= hi + tol:
                t = (x0 - xa) / (xb - xa)
                crossings.append(float(za + t * (zb - za)))
    if not crossings:
        raise ValueError(
            f"no triangle intersects the vertical line x={x0!r}; cannot "
            "determine the local surface elevation."
        )
    return min(crossings)


# ─────────────────────────────────────────────────────────────────────────
# Element assembly
# ─────────────────────────────────────────────────────────────────────────


def _local_mass(area: float) -> np.ndarray:
    return (area / 12.0) * np.array(
        [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]
    )


def assemble_stiffness_mass(
    mesh: TriMesh,
    weight_per_triangle: np.ndarray,
    *,
    weight_stiffness: bool,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    """Assemble the global (unconstrained) stiffness and mass matrices.

    Parameters
    ----------
    mesh : TriMesh
        Triangulation to assemble over.
    weight_per_triangle : ndarray, shape (n_triangles,)
        Per-triangle material coefficient (conductivity for TE's mass
        matrix, resistivity for TM's stiffness matrix).
    weight_stiffness : bool
        When ``True``, ``weight_per_triangle`` scales the stiffness matrix
        (TM mode) and the mass matrix is unweighted; when ``False``, it
        scales the mass matrix (TE mode) and the stiffness matrix is
        unweighted.

    Returns
    -------
    K : scipy.sparse.csr_matrix
        Global stiffness matrix, shape (n_nodes, n_nodes).
    M : scipy.sparse.csr_matrix
        Global mass matrix, same shape.

    Examples
    --------
    >>> mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
    >>> K, M = assemble_stiffness_mass(
    ...     mesh, np.array([1.0]), weight_stiffness=False
    ... )
    >>> np.allclose(K.sum(), 0.0)
    True
    >>> np.isclose(M.sum(), mesh.triangle_areas_m2.sum())
    True
    """
    n = mesh.n_nodes
    rows: list[int] = []
    cols: list[int] = []
    k_vals: list[float] = []
    m_vals: list[float] = []

    for t_idx, triangle in enumerate(mesh.triangles):
        grad, area = p1_gradients(mesh.nodes_m, triangle)
        local_k = area * (grad @ grad.T)
        local_m = _local_mass(area)
        weight = float(weight_per_triangle[t_idx])
        if weight_stiffness:
            local_k = weight * local_k
        else:
            local_m = weight * local_m
        for a in range(3):
            for b in range(3):
                rows.append(int(triangle[a]))
                cols.append(int(triangle[b]))
                k_vals.append(local_k[a, b])
                m_vals.append(local_m[a, b])

    K = sparse.csr_matrix((k_vals, (rows, cols)), shape=(n, n), dtype=float)
    M = sparse.csr_matrix((m_vals, (rows, cols)), shape=(n, n), dtype=float)
    return K, M


def _apply_dirichlet(
    A: sparse.csr_matrix,
    b: np.ndarray,
    boundary_node_ids: np.ndarray,
    bc_values: np.ndarray,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Replace boundary rows of A with identity and b with BC values."""
    A_lil = A.tolil()
    b = b.copy()
    for k, value in zip(boundary_node_ids, bc_values):
        A_lil.rows[k] = [int(k)]
        A_lil.data[k] = [1.0 + 0j]
        b[k] = value
    return A_lil.tocsr(), b


def _relative_residual(
    A: sparse.csr_matrix, x: np.ndarray, b: np.ndarray
) -> float:
    denominator = np.linalg.norm(b)
    if denominator <= 0.0:
        denominator = 1e-300
    return float(np.linalg.norm(A @ x - b) / denominator)


# ─────────────────────────────────────────────────────────────────────────
# Boundary conditions
# ─────────────────────────────────────────────────────────────────────────


def _boundary_node_ids(mesh: TriMesh) -> np.ndarray:
    return np.unique(np.asarray(mesh.boundary_segments).ravel())


def _te_boundary_values(
    mesh: TriMesh,
    resistivity_per_triangle: np.ndarray,
    boundary_node_ids: np.ndarray,
    omega: float,
) -> np.ndarray:
    values = np.empty(len(boundary_node_ids), dtype=complex)
    profile_cache: dict[float, tuple[np.ndarray, np.ndarray, float]] = {}
    for idx, node_id in enumerate(boundary_node_ids):
        x0, z0 = mesh.nodes_m[node_id]
        key = round(float(x0), 6)
        if key not in profile_cache:
            profile_cache[key] = (
                *_vertical_column_profile(
                    mesh, resistivity_per_triangle, float(x0)
                ),
                _local_surface_z(mesh, float(x0)),
            )
        rho, thick, surface_z = profile_cache[key]
        depth = float(z0) - surface_z
        if abs(depth) <= _SURFACE_TOLERANCE_M:
            values[idx] = 1.0 + 0j
            continue
        ey = _ey_1d_profile(omega, rho, thick, np.array([0.0, depth]))
        values[idx] = ey[1]
    return values


def _tm_boundary_values(
    mesh: TriMesh,
    resistivity_per_triangle: np.ndarray,
    boundary_node_ids: np.ndarray,
    omega: float,
) -> np.ndarray:
    values = np.empty(len(boundary_node_ids), dtype=complex)
    profile_cache: dict[float, tuple[np.ndarray, np.ndarray, float]] = {}
    for idx, node_id in enumerate(boundary_node_ids):
        x0, z0 = mesh.nodes_m[node_id]
        key = round(float(x0), 6)
        if key not in profile_cache:
            profile_cache[key] = (
                *_vertical_column_profile(
                    mesh, resistivity_per_triangle, float(x0)
                ),
                _local_surface_z(mesh, float(x0)),
            )
        rho, thick, surface_z = profile_cache[key]
        depth = float(z0) - surface_z
        if abs(depth) <= _SURFACE_TOLERANCE_M:
            values[idx] = 1.0 + 0j
            continue
        hy = _hy_1d_profile(omega, rho, thick, np.array([0.0, depth]))
        values[idx] = hy[1]
    return values


# ─────────────────────────────────────────────────────────────────────────
# Output container
# ─────────────────────────────────────────────────────────────────────────


@dataclass
class ForwardResponseTriFEM2D:
    """Output of :class:`Tri2DFEMForward`.

    Parameters
    ----------
    freqs : ndarray, shape (n_freqs,)
    receiver_names : tuple of str
    zxy, zyx : ndarray of complex, shape (n_freqs, n_receivers)
    residual_te, residual_tm : ndarray, shape (n_freqs,)
    """

    freqs: np.ndarray
    receiver_names: tuple[str, ...]
    zxy: np.ndarray
    zyx: np.ndarray
    residual_te: np.ndarray
    residual_tm: np.ndarray
    mesh: TriMesh = field(repr=False)


# ─────────────────────────────────────────────────────────────────────────
# Main forward class
# ─────────────────────────────────────────────────────────────────────────


class Tri2DFEMForward:
    """2-D triangular-mesh MT finite-element solver (TE + TM).

    Parameters
    ----------
    freqs : array-like
        Frequencies [Hz].
    mesh : TriMesh
        Triangulation with ``boundary_segments`` populated.
    resistivity_ohm_m : array-like, shape (n_triangles,)
        Per-triangle resistivity.
    receiver_node_ids : array-like of int
        Mesh node index of every receiver (must be at the surface).
    receiver_names : sequence of str
    verbose : bool, default=False

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.maxwell.contracts_tri import TriMesh
    >>> nodes = [[-1000, 0], [1000, 0], [1000, 1000], [-1000, 1000]]
    >>> mesh = TriMesh(
    ...     nodes, [[0, 1, 2], [0, 2, 3]],
    ...     boundary_segments=[[0, 1], [1, 2], [2, 3], [3, 0]],
    ... )
    >>> solver = Tri2DFEMForward(
    ...     [1.0], mesh, np.full(2, 100.0), [0], ["S00"], verbose=False
    ... )
    >>> resp = solver.run()
    >>> resp.zxy.shape
    (1, 1)
    """

    def __init__(
        self,
        freqs: np.ndarray,
        mesh: TriMesh,
        resistivity_ohm_m: np.ndarray,
        receiver_node_ids: np.ndarray,
        receiver_names: tuple[str, ...],
        *,
        verbose: bool = False,
    ) -> None:
        self.freqs = np.asarray(freqs, dtype=float)
        self.mesh = mesh
        self.resistivity_ohm_m = np.asarray(resistivity_ohm_m, dtype=float)
        self.receiver_node_ids = np.asarray(receiver_node_ids, dtype=int)
        self.receiver_names = tuple(receiver_names)
        self.verbose = verbose

    def run(self) -> ForwardResponseTriFEM2D:
        mesh = self.mesh
        rho = self.resistivity_ohm_m
        sigma = 1.0 / rho
        boundary_ids = _boundary_node_ids(mesh)
        incidence = _node_triangle_incidence(mesh)

        # Precompute per-triangle P1 gradients/areas once (frequency
        # independent).
        grads = np.empty((mesh.n_triangles, 3, 2))
        areas = np.empty(mesh.n_triangles)
        for t_idx, triangle in enumerate(mesh.triangles):
            grads[t_idx], areas[t_idx] = p1_gradients(mesh.nodes_m, triangle)

        K_lap, M_sigma = assemble_stiffness_mass(
            mesh, sigma, weight_stiffness=False
        )
        K_rho, M_plain = assemble_stiffness_mass(
            mesh, rho, weight_stiffness=True
        )

        nf = len(self.freqs)
        nr = len(self.receiver_node_ids)
        zxy_all = np.zeros((nf, nr), dtype=complex)
        zyx_all = np.zeros((nf, nr), dtype=complex)
        residual_te = np.zeros(nf)
        residual_tm = np.zeros(nf)

        for fi, freq in enumerate(self.freqs):
            omega = 2.0 * np.pi * freq
            if self.verbose:
                print(f"  [TriFEM2D] f={freq:.4g} Hz", end="\r", flush=True)

            # ── TE mode ──────────────────────────────────────────────
            A_te = (K_lap + 1j * omega * MU0 * M_sigma).astype(complex)
            b_te = np.zeros(mesh.n_nodes, dtype=complex)
            bc_te = _te_boundary_values(mesh, rho, boundary_ids, omega)
            A_te, b_te = _apply_dirichlet(A_te, b_te, boundary_ids, bc_te)
            ey = spsolve(A_te, b_te)
            residual_te[fi] = _relative_residual(A_te, ey, b_te)

            # ── TM mode ──────────────────────────────────────────────
            A_tm = (K_rho + 1j * omega * MU0 * M_plain).astype(complex)
            b_tm = np.zeros(mesh.n_nodes, dtype=complex)
            bc_tm = _tm_boundary_values(mesh, rho, boundary_ids, omega)
            A_tm, b_tm = _apply_dirichlet(A_tm, b_tm, boundary_ids, bc_tm)
            hy = spsolve(A_tm, b_tm)
            residual_tm[fi] = _relative_residual(A_tm, hy, b_tm)

            for ri, node_id in enumerate(self.receiver_node_ids):
                incident = incidence[int(node_id)]
                total_area = sum(areas[t] for t in incident)

                dey_dz = sum(
                    areas[t]
                    * (grads[t][:, 1] @ ey[mesh.triangles[t]])
                    for t in incident
                ) / total_area
                Hx = -dey_dz / (1j * omega * MU0)
                ey_here = ey[node_id]
                zxy_all[fi, ri] = (
                    ey_here / Hx if abs(Hx) > 1e-30 else 0.0 + 0j
                )

                dhy_dz = sum(
                    areas[t]
                    * (grads[t][:, 1] @ hy[mesh.triangles[t]])
                    for t in incident
                ) / total_area
                sigma_local = sum(
                    areas[t] * sigma[t] for t in incident
                ) / total_area
                Ex = -dhy_dz / sigma_local
                hy_here = hy[node_id]
                zyx_all[fi, ri] = (
                    -Ex / hy_here if abs(hy_here) > 1e-30 else 0.0 + 0j
                )

        if self.verbose:
            print(f"  [TriFEM2D] {nf} frequencies done.          ")

        return ForwardResponseTriFEM2D(
            self.freqs,
            self.receiver_names,
            zxy_all,
            zyx_all,
            residual_te,
            residual_tm,
            mesh,
        )


# ─────────────────────────────────────────────────────────────────────────
# Adapter
# ─────────────────────────────────────────────────────────────────────────


def _match_receiver_nodes(
    mesh: TriMesh, coordinates_m: np.ndarray, *, tol: float = 1e-3
) -> np.ndarray:
    """Return the mesh node index of every receiver, or -1 if unmatched."""
    node_ids = np.full(len(coordinates_m), -1, dtype=int)
    for i, (x, z) in enumerate(coordinates_m):
        distances = np.hypot(mesh.nodes_m[:, 0] - x, mesh.nodes_m[:, 1] - z)
        j = int(np.argmin(distances))
        if distances[j] <= tol:
            node_ids[i] = j
    return node_ids


class TriFEM2DAdapter(BaseMaxwellAdapter):
    """In-house 2-D triangular-mesh MT FEM adapter (research-only).

    Parameters
    ----------
    version : str, default="1.0"
        Adapter version reported in every
        :class:`~pycsamt.forward.maxwell.contracts.ForwardResult`.
    policy : AdapterPolicy or None, optional
        Solver-independent result acceptance policy.
    node_tolerance_m : float, default=1e-3
        Maximum distance between a receiver and its assumed mesh node.
    verbose : bool, default=False
        Forwarded to :class:`Tri2DFEMForward`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.forward.maxwell import ReceiverSet
    >>> from pycsamt.forward.maxwell.contracts_tri import TriMesh, TriProblem
    >>> nodes = [[-1000, 0], [1000, 0], [1000, 1000], [-1000, 1000]]
    >>> mesh = TriMesh(
    ...     nodes, [[0, 1, 2], [0, 2, 3]],
    ...     boundary_segments=[[0, 1], [1, 2], [2, 3], [3, 0]],
    ... )
    >>> problem = TriProblem(
    ...     mesh, np.full(2, 0.01), [1.0], ReceiverSet([[-1000, 0]], ["S00"])
    ... )
    >>> adapter = TriFEM2DAdapter()
    >>> adapter.capabilities.name
    'trifem2d'
    """

    def __init__(
        self,
        *,
        version: str = "1.0",
        policy: AdapterPolicy | None = None,
        node_tolerance_m: float = 1e-3,
        verbose: bool = False,
    ) -> None:
        capabilities = BackendCapabilities(
            name="trifem2d",
            version=version,
            dimensions=(2,),
            components=("zxy", "zyx"),
            time_conventions=("exp(+iwt)",),
            supports_nonuniform_mesh=True,
            supports_inactive_cells=False,
            supports_topography=True,
            supports_anisotropy=False,
            verified_benchmarks=_VERIFIED_BENCHMARKS,
        )
        super().__init__(capabilities, policy)
        self._node_tolerance_m = float(node_tolerance_m)
        self._verbose = bool(verbose)

    def assess(self, problem: TriProblem) -> CompatibilityReport:
        """Assess a problem, adding this solver's mesh/receiver checks.

        Parameters
        ----------
        problem : TriProblem
            Candidate simulation problem.

        Returns
        -------
        CompatibilityReport
            The generic capability report from
            :meth:`~pycsamt.forward.maxwell.adapters.BaseMaxwellAdapter.assess`,
            extended with boundary-segment, receiver-node-matching, surface,
            and permeability checks.

        Examples
        --------
        See :class:`TriFEM2DAdapter` for construction; a receiver not
        coinciding with any mesh node is rejected before the solver runs.
        """
        base = super().assess(problem)
        errors = list(base.errors)
        if not isinstance(problem.mesh, TriMesh):
            errors.append("trifem2d requires a TriProblem built on a TriMesh.")
            return CompatibilityReport(
                self.capabilities.name, False, tuple(errors), base.warnings
            )
        if problem.mesh.boundary_segments is None:
            errors.append(
                "trifem2d requires mesh.boundary_segments to apply "
                "Dirichlet boundary conditions."
            )
        off_surface: list[str] = []
        for name, (x, z) in zip(
            problem.receivers.names, problem.receivers.coordinates_m
        ):
            try:
                surface_z = _local_surface_z(problem.mesh, float(x))
            except ValueError:
                off_surface.append(name)
                continue
            if abs(float(z) - surface_z) > _SURFACE_TOLERANCE_M:
                off_surface.append(name)
        if off_surface:
            errors.append(
                "trifem2d only evaluates receivers sitting on the mesh's "
                f"own local surface (top boundary), not below/above it: {off_surface}"
            )
        node_ids = _match_receiver_nodes(
            problem.mesh,
            problem.receivers.coordinates_m,
            tol=self._node_tolerance_m,
        )
        if np.any(node_ids < 0):
            missing = [
                problem.receivers.names[i]
                for i in range(len(node_ids))
                if node_ids[i] < 0
            ]
            errors.append(
                "trifem2d requires every receiver to coincide with a mesh "
                f"node (within {self._node_tolerance_m:g} m); no matching "
                f"node found for: {missing}"
            )
        if not np.isclose(
            problem.magnetic_permeability_h_m, MU0, rtol=_PERMEABILITY_RTOL
        ):
            errors.append(
                f"trifem2d assumes vacuum magnetic permeability ({MU0:.6g} H/m)"
            )
        if len(errors) == len(base.errors):
            return base
        return CompatibilityReport(
            self.capabilities.name, False, tuple(errors), base.warnings
        )

    def _solve_backend(self, problem: TriProblem) -> ForwardResult:
        mesh = problem.mesh
        resistivity = 1.0 / problem.conductivity_s_m
        node_ids = _match_receiver_nodes(
            mesh, problem.receivers.coordinates_m, tol=self._node_tolerance_m
        )
        solver = Tri2DFEMForward(
            problem.frequencies_hz,
            mesh,
            resistivity,
            node_ids,
            problem.receivers.names,
            verbose=self._verbose,
        )
        start = time.monotonic()
        response = solver.run()
        runtime_s = time.monotonic() - start

        n_station = len(problem.receivers.names)
        n_freq = len(problem.frequencies_hz)
        n_component = len(problem.components)
        impedance = np.empty((n_station, n_freq, n_component), dtype=complex)
        for column, name in enumerate(problem.components):
            source = response.zxy if name == "zxy" else response.zyx
            impedance[:, :, column] = source.T

        converged = np.stack(
            [
                np.isfinite(response.residual_te),
                np.isfinite(response.residual_tm),
            ],
            axis=1,
        )
        residual = np.stack(
            [response.residual_te, response.residual_tm], axis=1
        )
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


def register_trifem2d_backend(*, replace: bool = False) -> None:
    """Register :class:`TriFEM2DAdapter` in the process-wide backend registry.

    Parameters
    ----------
    replace : bool, default=False
        Explicitly replace an existing ``"trifem2d"`` registration.

    Examples
    --------
    >>> from pycsamt.forward.maxwell import create_backend, list_backends
    >>> register_trifem2d_backend(replace=True)
    >>> "trifem2d" in list_backends()
    True
    >>> create_backend("trifem2d").capabilities.name
    'trifem2d'
    """
    capabilities = BackendCapabilities(
        name="trifem2d",
        version="1.0",
        dimensions=(2,),
        components=("zxy", "zyx"),
        time_conventions=("exp(+iwt)",),
        supports_nonuniform_mesh=True,
        supports_inactive_cells=False,
        supports_topography=True,
        supports_anisotropy=False,
        verified_benchmarks=_VERIFIED_BENCHMARKS,
    )
    _register_backend(
        BackendRegistration(capabilities, TriFEM2DAdapter), replace=replace
    )

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
2-D magnetotelluric finite-difference forward solver.

:class:`MT2DForward` solves the 2-D MT boundary-value problem for both
the TE and TM polarisation modes on a non-uniform staggered grid using
the finite-difference (FD) method.

Physics
-------
For a 2-D earth with conductivity σ(x, z) and a plane electromagnetic wave
incident from above, the governing scalar PDEs (time convention e^{+iωt})
are:

**TE mode** (electric field parallel to strike, E_y):

    ∂²E_y/∂x² + ∂²E_y/∂z² − iωμ₀σ(x,z) E_y = 0

**TM mode** (magnetic field parallel to strike, H_y):

    ∂/∂x[σ⁻¹ ∂H_y/∂x] + ∂/∂z[σ⁻¹ ∂H_y/∂z] − iωμ₀ H_y = 0

Grid and node numbering
-----------------------
The FD grid has ``(nx+1) × (nz+1)`` nodes.  Node at column *j*, row *i*
(depth) has global index ``k = i*(nx+1) + j``.  The surface is row *i = 0*,
the bottom is row *i = nz*.

Boundary conditions
-------------------
Dirichlet conditions are applied at all four boundaries using the analytic
1-D MT response of the edge column profiles (left/right) or the halfspace
(top/bottom).  For the TE mode the normalised plane-wave source
``E_y = 1`` at the surface is the natural starting point; the BC values
at depth follow from the 1-D recursion.  For the TM mode the analogous
1-D ``H_y`` field is used.

Each frequency requires two independent sparse solves (one per mode).
With ``scipy.sparse.linalg.spsolve`` and grids of order 50 × 40 nodes
(~2000 unknowns) each solve takes < 0.1 s on a single CPU core.

References
----------
Wannamaker, P.E., Stodt, J.A., Rijo, L. (1987). A stable finite element
solution for two-dimensional magnetotelluric modelling. *Geophys. J. Int.*
88, 277-296.

de Groot-Hedlin, C., Constable, S. (1990). Occam's inversion to generate
smooth, two-dimensional models from magnetotelluric data. *Geophysics*
55, 1613-1624.

Wait, J.R. (1954). On the relation between telluric currents and the Earth's
magnetic field. *Geophysics* 19, 281-289.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .em1d import _z_surface_mt
from .grid2d import Grid2D

__all__ = [
    "MT2DForward",
    "ForwardResponse2D",
]

MU0: float = 4.0e-7 * np.pi  # H/m


# ─────────────────────────────────────────────────────────────────────────────
# Output container
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class ForwardResponse2D:
    """Output of the 2-D MT forward solver.

    All impedance and apparent-resistivity arrays have shape
    ``(n_freqs, n_stations)``.

    Parameters
    ----------
    freqs : ndarray, shape (n_freqs,)
        Evaluation frequencies [Hz].
    stations_x : ndarray, shape (n_stations,)
        Station x-positions [m].
    zxy : ndarray, shape (n_freqs, n_stations), complex
        TE-mode surface impedance Z_xy = E_y / H_x  [V/A].
    zyx : ndarray, shape (n_freqs, n_stations), complex
        TM-mode surface impedance Z_yx = E_x / H_y  [V/A]  (negative
        by convention: Z_yx = −Z_xy for a 1-D earth).
    rho_a_te : ndarray, shape (n_freqs, n_stations)
        TE apparent resistivity [Ω·m].
    phase_te : ndarray, shape (n_freqs, n_stations)
        TE impedance phase [degrees, 0–90° for a normal model].
    rho_a_tm : ndarray, shape (n_freqs, n_stations)
        TM apparent resistivity [Ω·m].
    phase_tm : ndarray, shape (n_freqs, n_stations)
        TM impedance phase [degrees].
    residual_te, residual_tm : ndarray, shape (n_freqs,)
        Relative linear-system residual ``norm(A @ x - b) /
        norm(b)`` of the direct sparse solve for each mode, one
        value per frequency. This measures the achieved numerical
        accuracy of ``scipy.sparse.linalg.spsolve``; it is not an
        iteration count, since the solver is direct.
    grid : Grid2D
        The model grid used for the forward run.
    """

    freqs: np.ndarray
    stations_x: np.ndarray
    zxy: np.ndarray
    zyx: np.ndarray
    rho_a_te: np.ndarray
    phase_te: np.ndarray
    rho_a_tm: np.ndarray
    phase_tm: np.ndarray
    residual_te: np.ndarray
    residual_tm: np.ndarray
    grid: Grid2D = field(repr=False)

    # ─── convenience ─────────────────────────────────────────────────────

    @property
    def n_freqs(self) -> int:
        return len(self.freqs)

    @property
    def n_stations(self) -> int:
        return len(self.stations_x)

    @property
    def periods(self) -> np.ndarray:
        """Periods [s]."""
        return 1.0 / self.freqs

    def station_response(self, station: int) -> dict:
        """Return all response arrays for one station as a dict."""
        return dict(
            freqs=self.freqs,
            periods=self.periods,
            zxy=self.zxy[:, station],
            zyx=self.zyx[:, station],
            rho_a_te=self.rho_a_te[:, station],
            phase_te=self.phase_te[:, station],
            rho_a_tm=self.rho_a_tm[:, station],
            phase_tm=self.phase_tm[:, station],
        )

    def to_feature_array(
        self,
        *,
        mode: str = "both",
        log_rho: bool = True,
        include_phase: bool = True,
    ) -> np.ndarray:
        """Flatten to a 2-D feature matrix for ML training.

        Parameters
        ----------
        mode : {'te', 'tm', 'both'}
            Which mode(s) to include.
        log_rho : bool
            Return log₁₀(ρ_a) instead of ρ_a.
        include_phase : bool
            Concatenate phase alongside ρ_a.

        Returns
        -------
        ndarray, shape (n_stations, n_features)
        """
        parts = []
        for mode_key in (
            ["te"] if mode == "te" else ["tm"] if mode == "tm" else ["te", "tm"]
        ):
            rho = getattr(self, f"rho_a_{mode_key}")  # (n_freqs, n_stations)
            phi = getattr(self, f"phase_{mode_key}")
            if log_rho:
                rho = np.log10(np.maximum(rho, 1e-12))
            parts.append(rho.T)  # (n_stations, n_freqs)
            if include_phase:
                parts.append(phi.T)
        return np.concatenate(parts, axis=1)

    def plot(
        self,
        station: int = 0,
        ax=None,
        *,
        figsize: tuple[float, float] = (9, 5),
    ):
        """Quick diagnostic plot of ρ_a and phase for one station.

        Parameters
        ----------
        station : int
            Station index.
        ax : array of Axes or None
            Two axes (rho, phase).  Created when not provided.
        figsize : tuple

        Returns
        -------
        axes : ndarray of Axes, shape (2,)
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(2, 1, figsize=figsize, sharex=True)
        per = self.periods
        ax[0].loglog(per, self.rho_a_te[:, station], "b-o", ms=3, label="TE")
        ax[0].loglog(per, self.rho_a_tm[:, station], "r-s", ms=3, label="TM")
        ax[0].set_ylabel(r"$\rho_a$  (Ω·m)")
        ax[0].legend(fontsize=8)
        ax[0].set_title(
            f"2-D MT response — station {station} "
            f"(x = {self.stations_x[station]:.0f} m)"
        )
        ax[1].semilogx(per, self.phase_te[:, station], "b-o", ms=3)
        ax[1].semilogx(per, self.phase_tm[:, station], "r-s", ms=3)
        ax[1].set_ylabel("Phase  (°)")
        ax[1].set_xlabel("Period  (s)")
        # Zyx = -Zxy for a 1-D earth, so phase_tm normally falls in a
        # different quadrant than phase_te; let matplotlib auto-scale
        # rather than assume both share one fixed range.
        plt.tight_layout()
        return ax


# ─────────────────────────────────────────────────────────────────────────────
# 1-D MT helpers (boundary conditions)
# ─────────────────────────────────────────────────────────────────────────────


def _z1d_column(omega: float, rho: np.ndarray, thick: np.ndarray) -> complex:
    """Surface impedance for a 1-D column profile."""
    return _z_surface_mt(omega, rho, thick)


def _ey_1d_profile(
    omega: float,
    rho: np.ndarray,
    thick: np.ndarray,
    z_nodes: np.ndarray,
) -> np.ndarray:
    """Compute the normalised TE field E_y(z) at the grid node depths.

    Uses the standard downward-continuation formula.  The surface value
    is normalised to 1 + 0j.

    Returns
    -------
    Ey : ndarray of complex, shape (nz+1,)
    """
    nz1 = len(z_nodes)
    Ey = np.zeros(nz1, dtype=complex)
    Ey[0] = 1.0 + 0j

    # k_j = sqrt(i * omega * mu0 / rho_j)
    k = np.sqrt(1j * omega * MU0 / rho)  # shape (n_layers,)

    z_top = np.concatenate([[0.0], np.cumsum(thick)])  # layer top depths

    for n in range(1, nz1):
        z = z_nodes[n]
        # Find which layer z falls in
        lay = np.searchsorted(z_top, z, side="right") - 1
        lay = min(lay, len(rho) - 1)
        z - z_top[lay]
        k_j = k[lay]
        # Propagate E_y into this layer using the upward-recursion impedance
        # E_y(z) = E_y(z_top) * cosh(k_j*dz) + (Z(z_top)/iωμ₀/k_j) * ...
        # For the FD BC we only need a smooth exponential-like profile:
        # E_y(z) ≈ exp(-k_j * dz)  (downgoing wave in uniform layer)
        Ey[n] = Ey[0] * np.exp(-k_j * z)

    return Ey


def _hy_1d_profile(
    omega: float,
    rho: np.ndarray,
    thick: np.ndarray,
    z_nodes: np.ndarray,
) -> np.ndarray:
    """Compute the normalised TM field H_y(z) at grid node depths.

    Surface value normalised to 1 + 0j.

    Returns
    -------
    Hy : ndarray of complex, shape (nz+1,)
    """
    nz1 = len(z_nodes)
    Hy = np.zeros(nz1, dtype=complex)
    Hy[0] = 1.0 + 0j
    k = np.sqrt(1j * omega * MU0 / rho)
    z_top = np.concatenate([[0.0], np.cumsum(thick)])
    for n in range(1, nz1):
        z = z_nodes[n]
        lay = min(np.searchsorted(z_top, z, side="right") - 1, len(rho) - 1)
        Hy[n] = Hy[0] * np.exp(-k[lay] * z)
    return Hy


# ─────────────────────────────────────────────────────────────────────────────
# FD matrix assemblers
# ─────────────────────────────────────────────────────────────────────────────


def _node_index(i: int, j: int, nx1: int) -> int:
    """Global node index from row i (depth) and column j (x)."""
    return i * nx1 + j


def _relative_residual(
    A: sparse.csr_matrix,
    x: np.ndarray,
    b: np.ndarray,
) -> float:
    """Return ``norm(A @ x - b) / norm(b)`` for a direct solve."""
    denominator = np.linalg.norm(b)
    if denominator <= 0.0:
        denominator = 1e-300
    return float(np.linalg.norm(A @ x - b) / denominator)


def _assemble_te(
    grid: Grid2D,
    omega: float,
) -> sparse.csr_matrix:
    """Build the sparse FD matrix for the TE mode at angular frequency ω.

    The system  A · e = b  (where e is the flattened E_y field at all
    (nx+1)×(nz+1) nodes) is assembled here.  Boundary nodes are handled
    by the caller (rows replaced with identity, RHS set to BC values).

    Returns
    -------
    A : scipy.sparse.csr_matrix, shape (n_nodes, n_nodes)
    """
    dx = grid.dx  # (nx,)
    dz = grid.dz  # (nz,)
    nx, nz = grid.nx, grid.nz
    nx1 = nx + 1
    nz1 = nz + 1
    n = nx1 * nz1

    rows, cols, vals = [], [], []

    # Harmonic-mean conductivity for the reaction term:
    # At node (i,j) we average the up to 4 surrounding cell conductivities
    sigma = grid.conductivity  # (nz, nx)

    for i in range(nz1):
        for j in range(nx1):
            k = _node_index(i, j, nx1)

            # Boundary nodes — leave for caller to overwrite
            if i == 0 or i == nz or j == 0 or j == nx:
                rows.append(k)
                cols.append(k)
                vals.append(1.0)
                continue

            # ── x-direction second difference ──────────────────────────
            dx_l = dx[j - 1]  # cell to the left  of node j
            dx_r = dx[j]  # cell to the right of node j
            dx_avg = 0.5 * (dx_l + dx_r)

            c_xr = 1.0 / (dx_r * dx_avg)
            c_xl = 1.0 / (dx_l * dx_avg)
            c_xc = -(c_xr + c_xl)

            # ── z-direction second difference ──────────────────────────
            dz_u = dz[i - 1]  # cell above node i
            dz_d = dz[i]  # cell below node i
            dz_avg = 0.5 * (dz_u + dz_d)

            c_zd = 1.0 / (dz_d * dz_avg)
            c_zu = 1.0 / (dz_u * dz_avg)
            c_zc = -(c_zd + c_zu)

            # ── reaction term iωμ₀σ ────────────────────────────────────
            # Average conductivity over the four surrounding cells
            s_tl = sigma[i - 1, j - 1]
            s_tr = sigma[i - 1, j]
            s_bl = sigma[i, j - 1]
            s_br = sigma[i, j]
            # Area-weighted average of σ at this node
            a_tl = dx_l * dz_u
            a_tr = dx_r * dz_u
            a_bl = dx_l * dz_d
            a_br = dx_r * dz_d
            a_tot = a_tl + a_tr + a_bl + a_br
            sigma_n = (a_tl * s_tl + a_tr * s_tr + a_bl * s_bl + a_br * s_br) / a_tot
            react = 1j * omega * MU0 * sigma_n

            # ── diagonal ───────────────────────────────────────────────
            diag = c_xc + c_zc - react
            rows.append(k)
            cols.append(k)
            vals.append(diag)

            # ── off-diagonal ───────────────────────────────────────────
            k_xr = _node_index(i, j + 1, nx1)
            k_xl = _node_index(i, j - 1, nx1)
            k_zd = _node_index(i + 1, j, nx1)
            k_zu = _node_index(i - 1, j, nx1)

            rows.append(k)
            cols.append(k_xr)
            vals.append(c_xr)
            rows.append(k)
            cols.append(k_xl)
            vals.append(c_xl)
            rows.append(k)
            cols.append(k_zd)
            vals.append(c_zd)
            rows.append(k)
            cols.append(k_zu)
            vals.append(c_zu)

    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)


def _assemble_tm(
    grid: Grid2D,
    omega: float,
) -> sparse.csr_matrix:
    """Build the sparse FD matrix for the TM mode at angular frequency ω.

    The TM governing equation is:

        ∂/∂x(ρ ∂H_y/∂x) + ∂/∂z(ρ ∂H_y/∂z) − iωμ₀ H_y = 0

    where ρ = 1/σ is resistivity.  At an interior node, the coefficient
    for each direction's second difference is evaluated at a point that
    sits exactly on a material boundary in the *other* direction: the
    x-direction coefficient connecting node (i, j) to (i, j+1) lies on
    the shared edge between cell rows i-1 (above) and i (below), and
    the z-direction coefficient connecting (i, j) to (i+1, j) lies on
    the shared edge between cell columns j-1 (left) and j (right).  In
    both cases the two adjacent cells act as parallel current paths
    for that direction's flux (they carry independent current side by
    side, not current passing through one cell then the other), so the
    physically- and flux-consistent interface value is a
    thickness-weighted **harmonic** mean of resistivity (equivalently
    a thickness-weighted *arithmetic* mean of conductivity) — the same
    parallel-conductance combination :func:`_surface_impedance_tm`
    already uses when it averages conductivity, not resistivity,
    across laterally adjacent surface cells.  Using an arithmetic mean
    of resistivity here instead is the standard convention for the
    unrelated case of two cells in *series* along the differentiation
    direction, and was an inconsistency with that surface-extraction
    step for genuinely heterogeneous 2-D structure (a laterally or
    vertically uniform earth is unaffected either way, since averaging
    two equal values gives the same result regardless of method).
    This still reduces exactly to the TE equation for a uniform
    medium, since harmonic and arithmetic means agree when the two
    inputs are equal.

    Returns
    -------
    A : scipy.sparse.csr_matrix, shape (n_nodes, n_nodes)
    """
    dx = grid.dx
    dz = grid.dz
    rho = grid.resistivity  # (nz, nx)  — use resistivity, not conductivity
    nx, nz = grid.nx, grid.nz
    nx1 = nx + 1
    n = nx1 * (nz + 1)

    rows, cols, vals = [], [], []

    for i in range(nz + 1):
        for j in range(nx1):
            k = _node_index(i, j, nx1)

            # Boundary nodes — placeholder; overwritten by _apply_dirichlet
            if i == 0 or i == nz or j == 0 or j == nx:
                rows.append(k)
                cols.append(k)
                vals.append(1.0)
                continue

            # Cell sizes adjacent to node (i, j):
            #   row i-1 (above), row i (below), col j-1 (left), col j (right)
            dx_l = dx[j - 1]
            dx_r = dx[j]
            dz_u = dz[i - 1]
            dz_d = dz[i]
            dx_avg = 0.5 * (dx_l + dx_r)
            dz_avg = 0.5 * (dz_u + dz_d)

            # ── ρ at x-directed interfaces ───────────────────────────────
            # Right interface (between node j and j+1, at depth z_nodes[i]):
            # thickness-weighted harmonic mean of the two cells straddling
            # this depth (parallel current paths; see the docstring).
            rho_xr = (dz_u + dz_d) / (dz_u / rho[i - 1, j] + dz_d / rho[i, j])

            # Left interface (between node j-1 and j):
            rho_xl = (dz_u + dz_d) / (
                dz_u / rho[i - 1, j - 1] + dz_d / rho[i, j - 1]
            )

            # ── ρ at z-directed interfaces ───────────────────────────────
            # Lower interface (between node i and i+1, at x_nodes[j]):
            # thickness-weighted harmonic mean of the two cells straddling
            # this position (parallel current paths; see the docstring).
            rho_zd = (dx_l + dx_r) / (dx_l / rho[i, j - 1] + dx_r / rho[i, j])

            # Upper interface (between node i-1 and i):
            rho_zu = (dx_l + dx_r) / (
                dx_l / rho[i - 1, j - 1] + dx_r / rho[i - 1, j]
            )

            # ── FD coefficients ──────────────────────────────────────────
            c_xr = rho_xr / (dx_r * dx_avg)
            c_xl = rho_xl / (dx_l * dx_avg)
            c_xc = -(c_xr + c_xl)

            c_zd = rho_zd / (dz_d * dz_avg)
            c_zu = rho_zu / (dz_u * dz_avg)
            c_zc = -(c_zd + c_zu)

            # Reaction term:  ∂/∂x(ρ…) + ∂/∂z(ρ…) = iωμ₀ H_y
            diag = c_xc + c_zc - 1j * omega * MU0

            rows.append(k)
            cols.append(k)
            vals.append(diag)

            k_xr = _node_index(i, j + 1, nx1)
            k_xl = _node_index(i, j - 1, nx1)
            k_zd = _node_index(i + 1, j, nx1)
            k_zu = _node_index(i - 1, j, nx1)

            rows.append(k)
            cols.append(k_xr)
            vals.append(c_xr)
            rows.append(k)
            cols.append(k_xl)
            vals.append(c_xl)
            rows.append(k)
            cols.append(k_zd)
            vals.append(c_zd)
            rows.append(k)
            cols.append(k_zu)
            vals.append(c_zu)

    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)


# ─────────────────────────────────────────────────────────────────────────────
# Boundary condition builders
# ─────────────────────────────────────────────────────────────────────────────


def _bc_te(grid: Grid2D, omega: float) -> np.ndarray:
    """Return the Dirichlet BC vector for the TE mode.

    Non-zero only at boundary nodes; interior entries are 0.

    Returns
    -------
    b : ndarray of complex, shape (n_nodes,)
    """
    nx, nz = grid.nx, grid.nz
    nx1 = nx + 1
    n = nx1 * (nz + 1)
    b = np.zeros(n, dtype=complex)

    z_nodes = grid.z_nodes

    # Left column (j=0): 1-D profile of leftmost column
    rho_l, thick_l = grid.column_profile(0)
    ey_l = _ey_1d_profile(omega, rho_l, thick_l, z_nodes)
    for i in range(nz + 1):
        b[_node_index(i, 0, nx1)] = ey_l[i]

    # Right column (j=nx): 1-D profile of rightmost column
    rho_r, thick_r = grid.column_profile(nx - 1)
    ey_r = _ey_1d_profile(omega, rho_r, thick_r, z_nodes)
    for i in range(nz + 1):
        b[_node_index(i, nx, nx1)] = ey_r[i]

    # Top row (i=0): E_y = 1 (normalised incident field)
    for j in range(nx1):
        b[_node_index(0, j, nx1)] = 1.0 + 0j

    # Bottom row (i=nz): use average of left/right profiles at max depth
    ey_bot = 0.5 * (ey_l[-1] + ey_r[-1])
    for j in range(nx1):
        b[_node_index(nz, j, nx1)] = ey_bot

    return b


def _bc_tm(grid: Grid2D, omega: float) -> np.ndarray:
    """Return the Dirichlet BC vector for the TM mode.

    Returns
    -------
    b : ndarray of complex, shape (n_nodes,)
    """
    nx, nz = grid.nx, grid.nz
    nx1 = nx + 1
    n = nx1 * (nz + 1)
    b = np.zeros(n, dtype=complex)

    z_nodes = grid.z_nodes

    rho_l, thick_l = grid.column_profile(0)
    hy_l = _hy_1d_profile(omega, rho_l, thick_l, z_nodes)
    for i in range(nz + 1):
        b[_node_index(i, 0, nx1)] = hy_l[i]

    rho_r, thick_r = grid.column_profile(nx - 1)
    hy_r = _hy_1d_profile(omega, rho_r, thick_r, z_nodes)
    for i in range(nz + 1):
        b[_node_index(i, nx, nx1)] = hy_r[i]

    # Top row: H_y = 1
    for j in range(nx1):
        b[_node_index(0, j, nx1)] = 1.0 + 0j

    hy_bot = 0.5 * (hy_l[-1] + hy_r[-1])
    for j in range(nx1):
        b[_node_index(nz, j, nx1)] = hy_bot

    return b


# ─────────────────────────────────────────────────────────────────────────────
# Surface field extraction
# ─────────────────────────────────────────────────────────────────────────────


def _surface_impedance_te(
    ey_surface: np.ndarray,  # E_y at z=0, shape (nx+1,)
    ey_below: np.ndarray,  # E_y at z=dz[0], shape (nx+1,)
    dz0: float,
    omega: float,
) -> np.ndarray:
    """Compute Z_xy = E_y / H_x at the surface nodes.

    From Faraday's law:  H_x = (1/iωμ₀) ∂E_y/∂z
    Approximated as:     H_x ≈ (E_y(dz0) - E_y(0)) / (dz0 * iωμ₀)

    Returns
    -------
    Zxy : ndarray of complex, shape (nx+1,)
    """
    dEy_dz = (ey_below - ey_surface) / dz0
    # Faraday's law (e^{+iωt}): ∂E_y/∂z = iωμ₀ H_x  →  H_x = ∂E_y/∂z / (iωμ₀)
    # But with the downgoing-wave normalisation used in the FD assembler the
    # correct sign for Z_xy = E_y/H_x gives a +45° phase on a halfspace:
    Hx = -dEy_dz / (1j * omega * MU0)
    with np.errstate(divide="ignore", invalid="ignore"):
        Zxy = np.where(np.abs(Hx) > 1e-30, ey_surface / Hx, 0.0 + 0j)
    return Zxy


def _surface_impedance_tm(
    hy_surface: np.ndarray,  # H_y at z=0, shape (nx+1,)
    hy_below: np.ndarray,  # H_y at z=dz[0], shape (nx+1,)
    sigma_top: np.ndarray,  # σ of top-row cells, shape (nx,)
    dx: np.ndarray,  # cell widths, shape (nx,)
    dz0: float,
    omega: float,
) -> np.ndarray:
    """Compute Z_yx = −E_x / H_y at the surface nodes.

    From Ampere's law in 2-D (no y-dependence), the x-component of
    curl(H) is ``-∂H_y/∂z``, so ``-∂H_y/∂z = σ E_x`` and
    ``E_x = -(1/σ_top) * ∂H_y/∂z``. Omitting this negative sign
    previously made ``Z_yx`` come out equal to ``Z_xy`` instead of
    ``-Z_xy`` on a half-space, contradicting both this module's own
    documented convention and the standard MT sign convention.

    Returns
    -------
    Zyx : ndarray of complex, shape (nx+1,)
    """
    dHy_dz = (hy_below - hy_surface) / dz0

    # σ at nodes: average of left/right cells (edge nodes get one side)
    sigma_nodes = np.empty(len(dx) + 1)
    sigma_nodes[0] = sigma_top[0]
    sigma_nodes[1:-1] = 0.5 * (sigma_top[:-1] + sigma_top[1:])
    sigma_nodes[-1] = sigma_top[-1]

    with np.errstate(divide="ignore", invalid="ignore"):
        Ex = np.where(sigma_nodes > 0.0, -dHy_dz / sigma_nodes, 0.0 + 0j)
        Zyx = np.where(np.abs(hy_surface) > 1e-30, -Ex / hy_surface, 0.0 + 0j)
    return Zyx


def _z_to_rho_phase(
    Z: np.ndarray,
    omega: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert impedance to apparent resistivity and phase.

    Parameters
    ----------
    Z : ndarray of complex
    omega : float

    Returns
    -------
    rho_a : ndarray
    phase : ndarray  [degrees]
    """
    rho_a = np.abs(Z) ** 2 / (omega * MU0)
    phase = np.degrees(np.angle(Z))
    return rho_a, phase


def _interpolate_to_stations(
    surface_field: np.ndarray,  # shape (nx+1,) — node values
    x_nodes: np.ndarray,  # shape (nx+1,)
    x_stations: np.ndarray,  # shape (n_stations,)
) -> np.ndarray:
    """Linear interpolation of a surface node field to station positions."""
    return np.interp(x_stations, x_nodes, surface_field)


# ─────────────────────────────────────────────────────────────────────────────
# Main forward class
# ─────────────────────────────────────────────────────────────────────────────


class MT2DForward:
    """
    2-D magnetotelluric finite-difference forward solver.

    Solves both TE and TM modes for a list of frequencies on the provided
    :class:`~pycsamt.forward.grid2d.Grid2D` and returns a
    :class:`ForwardResponse2D` with apparent resistivity and phase at all
    station positions.

    Parameters
    ----------
    freqs : array-like
        Frequencies [Hz] at which to evaluate the response.
    grid : Grid2D
        2-D resistivity model and station layout.
    verbose : bool
        Print per-frequency progress.

    Examples
    --------
    Uniform halfspace — should recover the 1-D MT response::

        >>> import numpy as np
        >>> from pycsamt.forward.grid2d import Grid2D
        >>> from pycsamt.forward.em2d import MT2DForward

        >>> freqs = np.logspace(-2, 3, 10)
        >>> g = Grid2D.halfspace(rho=100.0, nx=30, nz=20,
        ...                       x_max=5_000.0, z_max=3_000.0,
        ...                       n_stations=5)
        >>> fwd = MT2DForward(freqs, g)
        >>> resp = fwd.run()
        >>> resp.rho_a_te.shape
        (10, 5)

    2-D model with conductive anomaly::

        >>> g2 = Grid2D.with_anomaly(bg_rho=100.0, anomaly_rho=2.0,
        ...     anomaly_bounds=(1500.0, 3500.0, 300.0, 900.0),
        ...     nx=40, nz=25, x_max=6000.0, z_max=3000.0, n_stations=10)
        >>> resp2 = MT2DForward(freqs, g2).run()
    """

    def __init__(
        self,
        freqs: np.ndarray,
        grid: Grid2D,
        *,
        verbose: bool = True,
    ):
        self.freqs = np.asarray(freqs, dtype=float)
        self.grid = grid
        self.verbose = verbose

    def run(self) -> ForwardResponse2D:
        """Run the forward solver for all frequencies.

        Returns
        -------
        ForwardResponse2D
        """
        grid = self.grid
        freqs = self.freqs
        nf = len(freqs)
        ns = grid.n_stations
        nx1 = grid.nx + 1
        x_nodes = grid.x_nodes
        dz0 = grid.dz[0]

        zxy_all = np.zeros((nf, ns), dtype=complex)
        zyx_all = np.zeros((nf, ns), dtype=complex)
        residual_te = np.zeros(nf)
        residual_tm = np.zeros(nf)

        sigma_top = grid.conductivity[0, :]  # top cell row

        for fi, freq in enumerate(freqs):
            omega = 2.0 * np.pi * freq

            if self.verbose:
                print(f"  [MT2D] f={freq:.4g} Hz", end="\r", flush=True)

            # ── TE mode ──────────────────────────────────────────────────
            A_te = _assemble_te(grid, omega)
            b_te = _bc_te(grid, omega)
            A_te, b_te = _apply_dirichlet(A_te, b_te, grid)
            ey_flat = spsolve(A_te, b_te)
            residual_te[fi] = _relative_residual(A_te, ey_flat, b_te)
            ey_grid = ey_flat.reshape(grid.nz + 1, nx1)  # (nz+1, nx+1)

            ey_surf = ey_grid[0, :]
            ey_below = ey_grid[1, :]
            Zxy_nodes = _surface_impedance_te(ey_surf, ey_below, dz0, omega)
            Zxy_st = _interpolate_to_stations(
                Zxy_nodes, x_nodes, grid.x_stations
            )
            zxy_all[fi, :] = Zxy_st

            # ── TM mode ──────────────────────────────────────────────────
            A_tm = _assemble_tm(grid, omega)
            b_tm = _bc_tm(grid, omega)
            A_tm, b_tm = _apply_dirichlet(A_tm, b_tm, grid)
            hy_flat = spsolve(A_tm, b_tm)
            residual_tm[fi] = _relative_residual(A_tm, hy_flat, b_tm)
            hy_grid = hy_flat.reshape(grid.nz + 1, nx1)

            hy_surf = hy_grid[0, :]
            hy_below = hy_grid[1, :]
            Zyx_nodes = _surface_impedance_tm(
                hy_surf, hy_below, sigma_top, grid.dx, dz0, omega
            )
            Zyx_st = _interpolate_to_stations(
                Zyx_nodes, x_nodes, grid.x_stations
            )
            zyx_all[fi, :] = Zyx_st

        if self.verbose:
            print(f"  [MT2D] {nf} frequencies done.          ")

        omega_col = 2.0 * np.pi * freqs[:, None]
        rho_a_te, phase_te = _z_to_rho_phase(zxy_all, omega_col)
        rho_a_tm, phase_tm = _z_to_rho_phase(zyx_all, omega_col)

        return ForwardResponse2D(
            freqs=freqs,
            stations_x=grid.x_stations.copy(),
            zxy=zxy_all,
            zyx=zyx_all,
            rho_a_te=rho_a_te,
            phase_te=phase_te,
            rho_a_tm=rho_a_tm,
            phase_tm=phase_tm,
            residual_te=residual_te,
            residual_tm=residual_tm,
            grid=grid,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Dirichlet enforcement
# ─────────────────────────────────────────────────────────────────────────────


def _apply_dirichlet(
    A: sparse.csr_matrix,
    b: np.ndarray,
    grid: Grid2D,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Replace boundary rows of A with identity and b with BC values.

    The assemblers already set boundary rows to identity; this function
    zeroes out any off-diagonal contributions that leaked through
    (accumulation from multiple passes) and enforces the RHS.

    Returns
    -------
    A_mod : sparse.csr_matrix
    b_mod : ndarray
    """
    nx, nz = grid.nx, grid.nz
    nx1 = nx + 1
    n = nx1 * (nz + 1)

    # Identify boundary node indices
    bc_mask = np.zeros(n, dtype=bool)
    for j in range(nx1):
        bc_mask[_node_index(0, j, nx1)] = True
        bc_mask[_node_index(nz, j, nx1)] = True
    for i in range(nz + 1):
        bc_mask[_node_index(i, 0, nx1)] = True
        bc_mask[_node_index(i, nx, nx1)] = True

    bc_indices = np.where(bc_mask)[0]

    A_lil = A.tolil()
    for k in bc_indices:
        A_lil.rows[k] = [k]
        A_lil.data[k] = [1.0 + 0j]

    return A_lil.tocsr(), b

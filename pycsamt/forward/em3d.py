# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Quasi-3D magnetotelluric forward solver.

:class:`MT3DForward` approximates the full 3-D MT impedance tensor by
running :class:`~pycsamt.forward.em2d.MT2DForward` on pairs of orthogonal
2-D cross-sections (XZ and YZ planes) extracted from a
:class:`~pycsamt.forward.grid3d.Grid3D` model.

Physics of the quasi-3D approximation
--------------------------------------
For a 3-D resistivity model σ(x, y, z), the full impedance tensor is:

    **Z** = [ Z_xx  Z_xy ]
            [ Z_yx  Z_yy ]

The quasi-3D method splits the tensor estimation into two independent
2-D FD solves:

1. **XZ profiles** — one per unique y-row of stations.  Each solve
   treats the model as 2-D with structure in the xz-plane.  From the
   2-D TE mode (E_y / H_x) we get *Z_xy*; from TM (E_x / H_y) we
   get *Z_yx*.

2. **YZ profiles** — one per unique x-column of stations.  The solver
   uses the yz-slice resistivity with y as the horizontal axis.  TE in
   this plane gives *Z_xy* from the NS structure; TM gives *Z_yx* from
   the NS structure.

The final quasi-3D tensor is the arithmetic average of contributions
from both profile directions:

    Z_xy ← ½ (Z_xy_xz + Z_xy_yz)
    Z_yx ← ½ (Z_yx_xz + Z_yx_yz)
    Z_xx = Z_yy = 0   (valid for structures without azimuthal asymmetry)

For a 1-D earth both profiles give the same response and the average
reduces to the exact 1-D solution.  For 2-D structures (variation in
one horizontal direction only), one profile captures the full signal
and the other gives the background, so the average is intermediate.
For genuine 3-D structures, both profiles contribute and the result is
an approximation that correctly captures the leading-order lateral
variation — sufficient for generating AI training data.

Accuracy and limitations
------------------------
* Off-diagonal components (Z_xx, Z_yy) are set to zero; the true 3-D
  Z_xx and Z_yy are typically 5–20 % of Z_xy in magnitude but are
  important for strike estimation.
* The approximation breaks down for strongly 3-D structures where
  galvanic coupling between adjacent columns is significant (e.g. thin
  resistive dykes).
* For the AI training use case, the quasi-3D data vastly improves on
  the fully independent-1D pseudo-3D baseline by injecting realistic
  lateral coupling physics.

References
----------
Wannamaker, P.E. (1999). Affordable magnetotellurics: Interpretation in
natural environments. *SEG Distinguished Instructors Short Course*, 7.

Mackie, R.L. et al. (1993). Three-dimensional electromagnetic modeling
using finite differences. *Geophysics*, 58(2), 215-226.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import (
    LinearOperator,
    bicgstab,
    spilu,
    spsolve,
)

from .batch import SurveyDataset3D
from .em2d import MT2DForward, _z_to_rho_phase
from .grid3d import Grid3D

__all__ = [
    "MT3DForward",
    "ForwardResponse3D",
]

# Warn when edge count exceeds this threshold
_FD3D_EDGE_WARN = 80_000

# Use direct solver (spsolve) for systems no larger than this.
# Below ~20 k edges spsolve is faster than ILU-BiCGSTAB and, crucially,
# converges at all frequencies (including low ω where skin depth >> grid).
_FD3D_DIRECT_THRESHOLD = 20_000

MU0: float = 4.0e-7 * np.pi


# ─────────────────────────────────────────────────────────────────────────────
# Output container
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class ForwardResponse3D:
    """Full (approximate) impedance tensor from the 3-D MT forward solver.

    All arrays have shape ``(n_freqs, n_stations)``.

    Parameters
    ----------
    freqs : ndarray, shape (n_freqs,)
    stations_xy : ndarray, shape (n_stations, 2)
        Station (x, y) positions [m].
    zxy, zyx, zxx, zyy : ndarray of complex
        Impedance tensor components [V/A].
    rho_a_xy, rho_a_yx, rho_a_xx, rho_a_yy : ndarray
        Apparent resistivities [Ω·m].
    phase_xy, phase_yx, phase_xx, phase_yy : ndarray
        Impedance phases [degrees].
    method : str
        Solver method (``'quasi3d'``).
    grid : Grid3D
        Source model.
    """

    freqs: np.ndarray
    stations_xy: np.ndarray

    zxy: np.ndarray
    zyx: np.ndarray
    zxx: np.ndarray
    zyy: np.ndarray

    rho_a_xy: np.ndarray
    phase_xy: np.ndarray
    rho_a_yx: np.ndarray
    phase_yx: np.ndarray
    rho_a_xx: np.ndarray
    phase_xx: np.ndarray
    rho_a_yy: np.ndarray
    phase_yy: np.ndarray

    method: str = "quasi3d"
    grid: Grid3D = field(repr=False, default=None)

    # ── convenience ──────────────────────────────────────────────────────────

    @property
    def n_freqs(self) -> int:
        return len(self.freqs)

    @property
    def n_stations(self) -> int:
        return len(self.stations_xy)

    @property
    def periods(self) -> np.ndarray:
        return 1.0 / self.freqs

    # ── feature extraction for ML training ───────────────────────────────────

    def to_feature_array(
        self,
        *,
        components: str = "xy_yx",
        log_rho: bool = True,
        include_phase: bool = True,
    ) -> np.ndarray:
        """Flatten to a 2-D feature matrix for ML input.

        Parameters
        ----------
        components : str
            Comma/underscore-separated component names.
            E.g. ``"xy_yx"`` (default) uses Z_xy and Z_yx;
            ``"xy"`` uses only Z_xy; ``"all"`` uses all four.
        log_rho : bool
            Return log₁₀(ρ_a) instead of ρ_a.
        include_phase : bool
            Concatenate phase alongside ρ_a.

        Returns
        -------
        X : ndarray, shape (n_stations, n_features)
        """
        if components == "all":
            comp_list = ["xy", "yx", "xx", "yy"]
        else:
            comp_list = [
                c.strip()
                for c in components.replace(",", "_").split("_")
                if c.strip()
            ]

        parts = []
        for c in comp_list:
            rho = getattr(self, f"rho_a_{c}")  # (n_freqs, n_stations)
            phi = getattr(self, f"phase_{c}")
            if log_rho:
                rho = np.log10(np.maximum(rho, 1e-12))
            parts.append(rho.T)  # (n_stations, n_freqs)
            if include_phase:
                parts.append(phi.T)
        return np.concatenate(parts, axis=1)

    def to_survey_dataset(
        self,
        y_models: np.ndarray | None = None,
        *,
        components: str = "xy_yx",
    ) -> SurveyDataset3D:
        """Wrap this response as a single-survey :class:`SurveyDataset3D`.

        The result can be concatenated with other surveys to build a
        training dataset for :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`.

        Parameters
        ----------
        y_models : ndarray or None
            Per-station model target vectors, shape ``(n_stations, n_params)``.
            When ``None`` a zero placeholder is used.
        components : str
            Passed to :meth:`to_feature_array`.

        Returns
        -------
        SurveyDataset3D
            Single-survey dataset (``n_surveys = 1``).
        """
        X = self.to_feature_array(components=components)  # (n_st, n_feat)
        if y_models is None:
            y_models = np.zeros((self.n_stations, 1), dtype=np.float32)
        X = X[None].astype(np.float32)  # (1, n_st, n_feat)
        y = y_models[None].astype(np.float32)  # (1, n_st, n_params)
        return SurveyDataset3D(
            X=X,
            y=y,
            coords=self.stations_xy.astype(np.float32),
            freqs=self.freqs.astype(np.float32),
            solver="quasi3d",
        )

    def station_response(self, station: int) -> dict:
        """Return all tensor components for one station as a dict."""
        return {
            "freqs": self.freqs,
            "periods": self.periods,
            "zxy": self.zxy[:, station],
            "zyx": self.zyx[:, station],
            "zxx": self.zxx[:, station],
            "zyy": self.zyy[:, station],
            "rho_a_xy": self.rho_a_xy[:, station],
            "phase_xy": self.phase_xy[:, station],
            "rho_a_yx": self.rho_a_yx[:, station],
            "phase_yx": self.phase_yx[:, station],
        }

    def __repr__(self) -> str:
        return (
            f"ForwardResponse3D(n_freqs={self.n_freqs}, "
            f"n_stations={self.n_stations}, method={self.method!r})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Forward solver
# ─────────────────────────────────────────────────────────────────────────────


class MT3DForward:
    """Quasi-3D magnetotelluric forward solver.

    Runs :class:`~pycsamt.forward.em2d.MT2DForward` on orthogonal XZ and
    YZ cross-sections extracted from a :class:`~pycsamt.forward.grid3d.Grid3D`
    model and assembles the result into an approximate full impedance tensor.

    Parameters
    ----------
    freqs : array-like
        Frequencies [Hz].
    grid : Grid3D
        3-D resistivity model with surface station positions.
    method : {'quasi3d'}
        Forward solver method.  Currently only ``'quasi3d'`` is
        implemented.  A true 3-D FD Yee-grid solver (``'fd3d'``) is
        planned for a future phase.
    verbose : bool
        Print per-profile progress.

    Examples
    --------
    Halfspace — should recover the 1-D MT response::

        >>> import numpy as np
        >>> from pycsamt.forward.grid3d import Grid3D
        >>> from pycsamt.forward.em3d import MT3DForward

        >>> freqs = np.logspace(-1, 2, 8)
        >>> g = Grid3D.halfspace(rho=100.0, nx=16, ny=16, nz=12,
        ...     x_max=5_000.0, y_max=5_000.0, z_max=3_000.0,
        ...     nx_stations=3, ny_stations=3)
        >>> resp = MT3DForward(freqs, g).run()
        >>> resp.rho_a_xy.shape
        (8, 9)

    3-D conductive block::

        >>> g = Grid3D.block_anomaly(bg_rho=500.0, anomaly_rho=5.0, nx=20, ny=20, nz=15,
        ...     x_max=6_000.0, y_max=6_000.0, z_max=4_000.0,
        ...     nx_stations=4, ny_stations=4)
        >>> resp = MT3DForward(freqs, g, verbose=False).run()
        >>> resp.rho_a_xy.std(axis=1).max() > 1.0   # lateral variation present
        True
    """

    def __init__(
        self,
        freqs: np.ndarray,
        grid: Grid3D,
        *,
        method: str = "quasi3d",
        verbose: bool = True,
    ):
        self.freqs = np.asarray(freqs, dtype=float)
        self.grid = grid
        self.method = method.lower()
        self.verbose = verbose

        _valid = {"quasi3d", "fd3d"}
        if self.method not in _valid:
            raise ValueError(
                f"method must be one of {_valid!r}, got {method!r}."
            )

    def run(self) -> ForwardResponse3D:
        """Run the forward solver and return the full impedance tensor.

        Returns
        -------
        ForwardResponse3D
        """
        if self.method == "fd3d":
            return self._run_fd3d()
        return self._run_quasi3d()

    # ─────────────────────────────────────────────────────────────────────────

    def _run_quasi3d(self) -> ForwardResponse3D:
        """Quasi-3D solver: orthogonal 2-D profile stacking."""
        grid = self.grid
        freqs = self.freqs
        ns = grid.n_stations
        nf = len(freqs)

        # Accumulation buffers for xz and yz profile contributions
        zxy_xz = np.zeros((nf, ns), dtype=complex)
        zyx_xz = np.zeros((nf, ns), dtype=complex)
        zxy_yz = np.zeros((nf, ns), dtype=complex)
        zyx_yz = np.zeros((nf, ns), dtype=complex)

        # ── XZ profiles  (one per unique y-row of stations) ──────────────────
        y_cells = grid._station_y_cells()
        unique_y_cells = np.unique(y_cells)
        n_xz_profiles = len(unique_y_cells)

        if self.verbose:
            print(
                f"  [MT3D quasi3d] {n_xz_profiles} XZ profiles + "
                f"{len(np.unique(grid._station_x_cells()))} YZ profiles …"
            )

        for ki, yi in enumerate(unique_y_cells):
            if self.verbose:
                print(
                    f"    XZ profile {ki + 1}/{n_xz_profiles}  (y-cell {yi})",
                    end="\r",
                    flush=True,
                )

            g2d_xz, st_idx = grid.xz_slice(yi)
            resp_xz = MT2DForward(freqs, g2d_xz, verbose=False).run()

            for local_i, global_i in enumerate(st_idx):
                zxy_xz[:, global_i] = resp_xz.zxy[:, local_i]
                zyx_xz[:, global_i] = resp_xz.zyx[:, local_i]

        # ── YZ profiles  (one per unique x-column of stations) ───────────────
        x_cells = grid._station_x_cells()
        unique_x_cells = np.unique(x_cells)
        n_yz_profiles = len(unique_x_cells)

        for ki, xi in enumerate(unique_x_cells):
            if self.verbose:
                print(
                    f"    YZ profile {ki + 1}/{n_yz_profiles}  (x-cell {xi})",
                    end="\r",
                    flush=True,
                )

            g2d_yz, st_idx = grid.yz_slice(xi)
            resp_yz = MT2DForward(freqs, g2d_yz, verbose=False).run()

            # In the YZ solver, horizontal axis = y_3D.
            # TE  → E-field perpendicular to yz-plane = E_x_3D
            #     → zxy_yz_solver = E_x_3D / H_y_3D = Z_xy_3D contribution
            # TM  → H-field perpendicular to yz-plane = H_x_3D
            #     → zyx_yz_solver = E_y_3D / H_x_3D = Z_yx_3D contribution
            for local_i, global_i in enumerate(st_idx):
                zxy_yz[:, global_i] = resp_yz.zxy[:, local_i]
                zyx_yz[:, global_i] = resp_yz.zyx[:, local_i]

        if self.verbose:
            print("  [MT3D quasi3d] done.                          ")

        # ── Assemble quasi-3D tensor ─────────────────────────────────────────
        # Average contributions from both profile directions
        zxy = 0.5 * (zxy_xz + zxy_yz)
        zyx = 0.5 * (zyx_xz + zyx_yz)

        # Off-diagonal: zero in the 2-D approximation
        zxx = np.zeros((nf, ns), dtype=complex)
        zyy = np.zeros((nf, ns), dtype=complex)

        omega = 2.0 * np.pi * freqs[:, None]  # (nf, 1) for broadcasting

        rho_a_xy, phase_xy = _z_to_rho_phase(zxy, omega)
        rho_a_yx, phase_yx = _z_to_rho_phase(zyx, omega)
        rho_a_xx, phase_xx = _z_to_rho_phase(zxx, omega)
        rho_a_yy, phase_yy = _z_to_rho_phase(zyy, omega)

        return ForwardResponse3D(
            freqs=freqs,
            stations_xy=grid.stations_xy.copy(),
            zxy=zxy,
            zyx=zyx,
            zxx=zxx,
            zyy=zyy,
            rho_a_xy=rho_a_xy,
            phase_xy=phase_xy,
            rho_a_yx=rho_a_yx,
            phase_yx=phase_yx,
            rho_a_xx=rho_a_xx,
            phase_xx=phase_xx,
            rho_a_yy=rho_a_yy,
            phase_yy=phase_yy,
            method=self.method,
            grid=grid,
        )

    # ─────────────────────────────────────────────────────────────────────────
    # True 3-D FD Yee-grid solver
    # ─────────────────────────────────────────────────────────────────────────

    def _run_fd3d(self) -> ForwardResponse3D:
        """True 3-D Yee-grid FD solver.

        Builds the curl-curl operator A = C^T C + iωμ₀ M_σ on the Yee
        staggered grid, enforces Dirichlet boundary conditions from the
        1-D MT solution of the corner column profiles, and solves for both
        x- and y-polarised E-fields using BiCGSTAB + ILU preconditioner.
        The full impedance tensor Z is assembled from the two polarisation
        solutions at each station.

        Accuracy note
        -------------
        Grid-resolution error is O(Δh) — the same as the quasi-3D solver
        but now capturing true 3-D galvanic coupling and non-zero Z_xx /
        Z_yy components.  For grids where the skin depth spans < 4 cells
        the error increases; use at least 6 cells per skin depth for
        quantitative work.
        """
        grid = self.grid
        freqs = self.freqs
        nx, ny, nz = grid.nx, grid.ny, grid.nz
        ns = grid.n_stations
        nf = len(freqs)

        n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
        n_edges = n_ex + n_ey + n_ez

        if n_edges > _FD3D_EDGE_WARN:
            warnings.warn(
                f"FD3D: grid has {n_edges} edges (nx={nx}, ny={ny}, nz={nz}). "
                f"Grids with > {_FD3D_EDGE_WARN} edges are slow in pure Python. "
                "Consider using method='quasi3d' or reducing the grid size.",
                RuntimeWarning,
                stacklevel=3,
            )

        if self.verbose:
            print(
                f"  [MT3D fd3d] {n_edges} unknowns  "
                f"(nx={nx}, ny={ny}, nz={nz})  {nf} freq …"
            )

        # ── Frequency-independent assembly ────────────────────────────────────
        C = _build_curl_matrix(grid)
        sigma = _edge_conductivity(grid)  # shape (n_edges,)
        M_sig = sparse.diags(sigma, format="csr")

        bc_mask = _boundary_edge_mask(nx, ny, nz)  # boolean (n_edges,)
        bc_idx = np.where(bc_mask)[0]

        # ── Storage ──────────────────────────────────────────────────────────
        zxy_all = np.zeros((nf, ns), dtype=complex)
        zyx_all = np.zeros((nf, ns), dtype=complex)
        zxx_all = np.zeros((nf, ns), dtype=complex)
        zyy_all = np.zeros((nf, ns), dtype=complex)

        for fi, freq in enumerate(freqs):
            omega = 2.0 * np.pi * freq
            if self.verbose:
                print(
                    f"  [MT3D fd3d] freq {fi + 1}/{nf}  {freq:.4g} Hz",
                    end="\r",
                    flush=True,
                )

            # System matrix A = C^T C + iωμ₀ σ
            A_base = C.T @ C + 1j * omega * MU0 * M_sig  # (n_edges, n_edges)

            # ── Solve two polarisations ───────────────────────────────────────
            E_xpol = self._fd3d_solve(
                A_base, bc_mask, bc_idx, grid, omega, "x"
            )
            E_ypol = self._fd3d_solve(
                A_base, bc_mask, bc_idx, grid, omega, "y"
            )

            # ── Extract impedance tensor at surface stations ──────────────────
            zxx, zxy, zyx, zyy = _extract_z_tensor_fd3d(
                E_xpol, E_ypol, grid, omega, n_ex, n_ey, n_ez
            )
            zxx_all[fi] = zxx
            zxy_all[fi] = zxy
            zyx_all[fi] = zyx
            zyy_all[fi] = zyy

        if self.verbose:
            print("  [MT3D fd3d] done.                          ")

        omega2d = 2.0 * np.pi * freqs[:, None]
        rho_a_xy, phase_xy = _z_to_rho_phase(zxy_all, omega2d)
        rho_a_yx, phase_yx = _z_to_rho_phase(zyx_all, omega2d)
        rho_a_xx, phase_xx = _z_to_rho_phase(zxx_all, omega2d)
        rho_a_yy, phase_yy = _z_to_rho_phase(zyy_all, omega2d)

        return ForwardResponse3D(
            freqs=freqs,
            stations_xy=grid.stations_xy.copy(),
            zxy=zxy_all,
            zyx=zyx_all,
            zxx=zxx_all,
            zyy=zyy_all,
            rho_a_xy=rho_a_xy,
            phase_xy=phase_xy,
            rho_a_yx=rho_a_yx,
            phase_yx=phase_yx,
            rho_a_xx=rho_a_xx,
            phase_xx=phase_xx,
            rho_a_yy=rho_a_yy,
            phase_yy=phase_yy,
            method="fd3d",
            grid=grid,
        )

    @staticmethod
    def _fd3d_solve(
        A: sparse.csr_matrix,
        bc_mask: np.ndarray,
        bc_idx: np.ndarray,
        grid: Grid3D,
        omega: float,
        pol: str,
    ) -> np.ndarray:
        """Total-field solve — returns E at all Yee-grid edges.

        Enforces Dirichlet BCs from the 1-D MT background plane-wave field
        (computed by :func:`_bc_vector_3d`) and solves the interior with
        BiCGSTAB + ILU(A_mod).  The preconditioner is built from A_mod
        (after Dirichlet enforcement) rather than the original A so that
        the boundary-modified structure is correctly approximated.
        """
        b = _bc_vector_3d(grid, omega, pol)
        A_mod, b_mod = _apply_dirichlet_3d(A, b, bc_mask, bc_idx)

        # Direct solve for small grids — robust at all frequencies because the
        # skin-depth can greatly exceed the grid at low ω, leaving C^T C
        # nearly singular and iterative solvers ill-conditioned.
        if A_mod.shape[0] <= _FD3D_DIRECT_THRESHOLD:
            return spsolve(A_mod, b_mod)

        prec = None
        try:
            ilu = spilu(A_mod.tocsc(), drop_tol=1e-6, fill_factor=20)
            prec = LinearOperator(A_mod.shape, ilu.solve)
        except Exception:
            pass

        if prec is not None:
            E, info = bicgstab(A_mod, b_mod, M=prec, tol=1e-7, maxiter=800)
        else:
            E, info = bicgstab(A_mod, b_mod, tol=1e-7, maxiter=800)
        if info != 0:
            warnings.warn(
                f"FD3D BiCGSTAB did not converge (pol={pol!r}, info={info}). "
                "Result may be inaccurate.  Try coarser grid or fewer freqs.",
                RuntimeWarning,
                stacklevel=4,
            )
        return E


# ═══════════════════════════════════════════════════════════════════════════════
# FD3D internal helpers  (module-level, used by MT3DForward._run_fd3d)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Edge / face dimension helpers ─────────────────────────────────────────────


def _fd3d_edge_counts(nx: int, ny: int, nz: int) -> tuple[int, int, int]:
    """Return (n_ex, n_ey, n_ez) edge counts for a nx×ny×nz cell grid."""
    n_ex = nx * (ny + 1) * (nz + 1)
    n_ey = (nx + 1) * ny * (nz + 1)
    n_ez = (nx + 1) * (ny + 1) * nz
    return n_ex, n_ey, n_ez


def _fd3d_face_counts(nx: int, ny: int, nz: int) -> tuple[int, int, int]:
    n_hx = (nx + 1) * ny * nz
    n_hy = nx * (ny + 1) * nz
    n_hz = nx * ny * (nz + 1)
    return n_hx, n_hy, n_hz


# ── Vectorised edge global-index helpers ─────────────────────────────────────


def _ex_idx(I, J, K, ny, nz) -> np.ndarray:
    """Global index of Ex edge at (I,J,K).  I=0..nx-1, J=0..ny, K=0..nz."""
    return I * (ny + 1) * (nz + 1) + J * (nz + 1) + K


def _ey_idx(I, J, K, nx, ny, nz) -> np.ndarray:
    """Global index of Ey edge at (I,J,K).  I=0..nx, J=0..ny-1, K=0..nz."""
    n_ex = nx * (ny + 1) * (nz + 1)
    return n_ex + I * ny * (nz + 1) + J * (nz + 1) + K


def _ez_idx(I, J, K, nx, ny, nz) -> np.ndarray:
    """Global index of Ez edge at (I,J,K).  I=0..nx, J=0..ny, K=0..nz-1."""
    n_ex = nx * (ny + 1) * (nz + 1)
    n_ey = (nx + 1) * ny * (nz + 1)
    return n_ex + n_ey + I * (ny + 1) * nz + J * nz + K


# ── Curl matrix ───────────────────────────────────────────────────────────────


def _build_curl_matrix(grid: Grid3D) -> sparse.csr_matrix:
    """Build the discrete curl operator C (n_faces × n_edges).

    Each row corresponds to one magnetic face (Hx, Hy, or Hz) and has
    exactly 4 non-zero entries corresponding to the four E-field edges
    that surround the face.  Entries include the 1/Δ spacing factors so
    that C represents the true continuous curl in SI units.

    Returns
    -------
    C : scipy.sparse.csr_matrix, shape (n_faces, n_edges)
    """
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz

    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    n_hx, n_hy, n_hz = _fd3d_face_counts(nx, ny, nz)
    n_faces = n_hx + n_hy + n_hz
    n_edges = n_ex + n_ey + n_ez

    # Each face has 4 non-zero entries → max 4 × n_faces entries total
    row_list, col_list, val_list = [], [], []

    def _add(rows, cols, vals):
        row_list.append(rows)
        col_list.append(cols)
        val_list.append(vals)

    # ── Hx faces  (curl E)_x = ∂Ez/∂y − ∂Ey/∂z ─────────────────────────────
    # i=0..nx,  j=0..ny-1,  k=0..nz-1
    ii, jj, kk = np.meshgrid(
        np.arange(nx + 1), np.arange(ny), np.arange(nz), indexing="ij"
    )
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()
    row_hx = ii * ny * nz + jj * nz + kk  # row in C

    # +Ez[i, j+1, k] / dy[j]
    _add(row_hx, _ez_idx(ii, jj + 1, kk, nx, ny, nz), 1.0 / dy[jj])
    # -Ez[i, j,   k] / dy[j]
    _add(row_hx, _ez_idx(ii, jj, kk, nx, ny, nz), -1.0 / dy[jj])
    # -Ey[i, j,   k+1] / dz[k]
    _add(row_hx, _ey_idx(ii, jj, kk + 1, nx, ny, nz), -1.0 / dz[kk])
    # +Ey[i, j,   k]   / dz[k]
    _add(row_hx, _ey_idx(ii, jj, kk, nx, ny, nz), 1.0 / dz[kk])

    # ── Hy faces  (curl E)_y = ∂Ex/∂z − ∂Ez/∂x ─────────────────────────────
    # i=0..nx-1,  j=0..ny,  k=0..nz-1
    ii, jj, kk = np.meshgrid(
        np.arange(nx), np.arange(ny + 1), np.arange(nz), indexing="ij"
    )
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()
    row_hy = n_hx + ii * (ny + 1) * nz + jj * nz + kk

    # +Ex[i, j, k+1] / dz[k]
    _add(row_hy, _ex_idx(ii, jj, kk + 1, ny, nz), 1.0 / dz[kk])
    # -Ex[i, j, k]   / dz[k]
    _add(row_hy, _ex_idx(ii, jj, kk, ny, nz), -1.0 / dz[kk])
    # -Ez[i+1, j, k] / dx[i]
    _add(row_hy, _ez_idx(ii + 1, jj, kk, nx, ny, nz), -1.0 / dx[ii])
    # +Ez[i,   j, k] / dx[i]
    _add(row_hy, _ez_idx(ii, jj, kk, nx, ny, nz), 1.0 / dx[ii])

    # ── Hz faces  (curl E)_z = ∂Ey/∂x − ∂Ex/∂y ─────────────────────────────
    # i=0..nx-1,  j=0..ny-1,  k=0..nz
    ii, jj, kk = np.meshgrid(
        np.arange(nx), np.arange(ny), np.arange(nz + 1), indexing="ij"
    )
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()
    row_hz = n_hx + n_hy + ii * ny * (nz + 1) + jj * (nz + 1) + kk

    # +Ey[i+1, j, k] / dx[i]
    _add(row_hz, _ey_idx(ii + 1, jj, kk, nx, ny, nz), 1.0 / dx[ii])
    # -Ey[i,   j, k] / dx[i]
    _add(row_hz, _ey_idx(ii, jj, kk, nx, ny, nz), -1.0 / dx[ii])
    # -Ex[i, j+1, k] / dy[j]
    _add(row_hz, _ex_idx(ii, jj + 1, kk, ny, nz), -1.0 / dy[jj])
    # +Ex[i, j,   k] / dy[j]
    _add(row_hz, _ex_idx(ii, jj, kk, ny, nz), 1.0 / dy[jj])

    rows = np.concatenate(row_list)
    cols = np.concatenate(col_list)
    vals = np.concatenate(val_list)

    return sparse.coo_matrix(
        (vals, (rows, cols)), shape=(n_faces, n_edges)
    ).tocsr()


# ── Edge conductivity (mass matrix diagonal) ─────────────────────────────────


def _edge_conductivity(grid: Grid3D) -> np.ndarray:
    """Return conductivity values at every edge midpoint.

    For each edge, the conductivity is the arithmetic mean of the
    conductivities of the cells surrounding that edge (1–4 cells
    depending on edge position).

    Returns
    -------
    sigma_edge : ndarray, shape (n_edges,)
    """
    _dx, _dy, _dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    sigma = grid.conductivity  # (nz, ny, nx)

    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    n_edges = n_ex + n_ey + n_ez
    sig_edge = np.zeros(n_edges)

    # ── Ex edges  (i=0..nx-1, j=0..ny, k=0..nz) ─────────────────────────────
    # Adjacent cells in j-k plane (up to 4)
    for j in range(ny + 1):
        for k in range(nz + 1):
            cells = []
            if j > 0 and k > 0:
                cells.append(sigma[k - 1, j - 1, :])
            if j < ny and k > 0:
                cells.append(sigma[k - 1, j, :])
            if j > 0 and k < nz:
                cells.append(sigma[k, j - 1, :])
            if j < ny and k < nz:
                cells.append(sigma[k, j, :])
            s = np.mean(cells, axis=0) if cells else np.zeros(nx)
            idx = _ex_idx(np.arange(nx), j, k, ny, nz)
            sig_edge[idx] = s

    # ── Ey edges  (i=0..nx, j=0..ny-1, k=0..nz) ─────────────────────────────
    for i in range(nx + 1):
        for k in range(nz + 1):
            cells = []
            if i > 0 and k > 0:
                cells.append(sigma[k - 1, :, i - 1])
            if i < nx and k > 0:
                cells.append(sigma[k - 1, :, i])
            if i > 0 and k < nz:
                cells.append(sigma[k, :, i - 1])
            if i < nx and k < nz:
                cells.append(sigma[k, :, i])
            s = np.mean(cells, axis=0) if cells else np.zeros(ny)
            idx = _ey_idx(i, np.arange(ny), k, nx, ny, nz)
            sig_edge[idx] = s

    # ── Ez edges  (i=0..nx, j=0..ny, k=0..nz-1) ─────────────────────────────
    for i in range(nx + 1):
        for j in range(ny + 1):
            cells = []
            if i > 0 and j > 0:
                cells.append(sigma[:, j - 1, i - 1])
            if i < nx and j > 0:
                cells.append(sigma[:, j - 1, i])
            if i > 0 and j < ny:
                cells.append(sigma[:, j, i - 1])
            if i < nx and j < ny:
                cells.append(sigma[:, j, i])
            s = np.mean(cells, axis=0) if cells else np.zeros(nz)
            idx = _ez_idx(i, j, np.arange(nz), nx, ny, nz)
            sig_edge[idx] = s

    return sig_edge


# ── Boundary edge mask ────────────────────────────────────────────────────────


def _boundary_edge_mask(nx: int, ny: int, nz: int) -> np.ndarray:
    """Return a boolean mask of length n_edges; True for boundary edges.

    Boundary edges are those lying on any of the 6 grid faces:
    * Top/bottom (z=0, z=z_max):  Ex and Ey edges with k=0 or k=nz
    * Left/right (x=0, x=x_max):  Ey and Ez edges with i=0 or i=nx
    * Front/back (y=0, y=y_max):  Ex and Ez edges with j=0 or j=ny
    """
    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    n_edges = n_ex + n_ey + n_ez
    mask = np.zeros(n_edges, dtype=bool)

    # ── Ex boundary: j=0, j=ny, k=0, k=nz ───────────────────────────────────
    ii = np.arange(nx)
    for j in (0, ny):
        for k in range(nz + 1):
            mask[_ex_idx(ii, j, k, ny, nz)] = True
    for j in range(ny + 1):
        for k in (0, nz):
            mask[_ex_idx(ii, j, k, ny, nz)] = True

    # ── Ey boundary: i=0, i=nx, k=0, k=nz ───────────────────────────────────
    jj = np.arange(ny)
    for i in (0, nx):
        for k in range(nz + 1):
            mask[_ey_idx(i, jj, k, nx, ny, nz)] = True
    for i in range(nx + 1):
        for k in (0, nz):
            mask[_ey_idx(i, jj, k, nx, ny, nz)] = True

    # ── Ez boundary: i=0, i=nx, j=0, j=ny ───────────────────────────────────
    kk = np.arange(nz)
    for i in (0, nx):
        for j in range(ny + 1):
            mask[_ez_idx(i, j, kk, nx, ny, nz)] = True
    for i in range(nx + 1):
        for j in (0, ny):
            mask[_ez_idx(i, j, kk, nx, ny, nz)] = True

    return mask


# ── Boundary condition vector ─────────────────────────────────────────────────


def _bc_vector_3d(grid: Grid3D, omega: float, pol: str) -> np.ndarray:
    """Return the Dirichlet BC vector for one polarisation.

    For **x-polarisation** (``pol='x'``):
    * Top (k=0): Ex = 1 (incident normalised plane wave), Ey = 0.
    * All other boundary Ex/Ey: Ex = exp(−k_1d z), Ey = 0.
    * All Ez boundaries: Ez = 0.

    For **y-polarisation** (``pol='y'``): same but Ex↔Ey roles swapped.

    The 1-D decaying factor uses the average conductivity of the four
    corner columns as a representative background.

    Returns
    -------
    b : ndarray of complex, shape (n_edges,)
    """
    _dx, _dy, _dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    sigma = grid.conductivity  # (nz, ny, nx)
    z_nodes = grid.z_nodes  # (nz+1,)

    n_ex, n_ey, n_ez = _fd3d_edge_counts(nx, ny, nz)
    b = np.zeros(n_ex + n_ey + n_ez, dtype=complex)

    # Average conductivity of 4 corner columns for the 1-D decay
    sigma_corners = np.mean(
        [
            sigma[:, 0, 0],
            sigma[:, 0, -1],
            sigma[:, -1, 0],
            sigma[:, -1, -1],
        ],
        axis=0,
    )  # shape (nz,)
    sigma_avg = float(np.mean(sigma_corners))
    k_1d = np.sqrt(1j * omega * MU0 * sigma_avg)  # complex wavenumber

    def E_decay(z):
        """Normalised 1-D plane-wave E-field at depth z [m]."""
        return np.exp(-k_1d * z)

    pol = pol.lower()
    ii = np.arange(nx)
    jj = np.arange(ny)

    if pol == "x":
        # Ex edges: boundary j=0/ny, k=0/nz
        # Top surface k=0: Ex = 1 (all j)
        for j in range(ny + 1):
            b[_ex_idx(ii, j, 0, ny, nz)] = 1.0
        # Bottom k=nz and side j=0/ny
        for j in (0, ny):
            for k in range(1, nz + 1):
                b[_ex_idx(ii, j, k, ny, nz)] = E_decay(z_nodes[k])
        for j in range(1, ny):
            for k in (nz,):
                b[_ex_idx(ii, j, k, ny, nz)] = E_decay(z_nodes[k])
        # Ey edges: all boundary → 0 (already zero)
        # Ez edges: all boundary → 0 (already zero)

    else:  # y-polarisation
        # Ey edges: boundary k=0 top → Ey = 1
        for i in range(nx + 1):
            b[_ey_idx(i, jj, 0, nx, ny, nz)] = 1.0
        for i in (0, nx):
            for k in range(1, nz + 1):
                b[_ey_idx(i, jj, k, nx, ny, nz)] = E_decay(z_nodes[k])
        for i in range(1, nx):
            for k in (nz,):
                b[_ey_idx(i, jj, k, nx, ny, nz)] = E_decay(z_nodes[k])
        # Ex, Ez: all boundary → 0

    return b


# ── Dirichlet enforcement ─────────────────────────────────────────────────────


def _apply_dirichlet_3d(
    A: sparse.csr_matrix,
    b: np.ndarray,
    bc_mask: np.ndarray,
    bc_idx: np.ndarray,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Enforce Dirichlet BCs symmetrically.

    Uses the projection ``A_mod = P A P + I_bc`` where *P* is a diagonal
    restriction that zeros out boundary rows **and** columns, and *I_bc*
    restores the identity diagonal at boundary positions.  This preserves
    the Hermitian structure of the curl-curl operator, keeping the ILU
    preconditioner well-conditioned.

    The RHS is updated to:
        b_mod[interior] = b[interior] - A[interior, bc] @ E_bc
        b_mod[boundary] = E_bc            (boundary values)

    For a uniform halfspace E_bc is the exact 1-D MT solution, so the
    source b_mod is smooth and BiCGSTAB converges rapidly.
    """
    # ── Subtract known boundary contributions from interior RHS ─────────────
    b_bc = np.zeros_like(b)
    b_bc[bc_idx] = b[bc_idx]
    b_mod = b.copy() - A @ b_bc
    b_mod[bc_idx] = b[bc_idx]

    # ── Projection: zero all boundary rows AND columns, add I at bc ─────────
    int_w = (~bc_mask).astype(np.float64)
    P = sparse.diags(int_w, format="csr")
    I_bc = sparse.diags((bc_mask).astype(np.float64), format="csr")
    A_mod = (P @ A @ P + I_bc).tocsr()
    return A_mod, b_mod


# ── Surface impedance extraction ──────────────────────────────────────────────


def _extract_z_tensor_fd3d(
    E_xpol: np.ndarray,
    E_ypol: np.ndarray,
    grid: Grid3D,
    omega: float,
    n_ex: int,
    n_ey: int,
    n_ez: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract the full impedance tensor at surface stations.

    Uses the two polarisation solutions to assemble:

        Z = [[Z_xx, Z_xy], [Z_yx, Z_yy]]

    via the linear system:

        [E_x^xpol  E_x^ypol]   [Z_xx  Z_xy] [H_x^xpol  H_x^ypol]
        [E_y^xpol  E_y^ypol] = [Z_yx  Z_yy] [H_y^xpol  H_y^ypol]

    Surface H is approximated from the finite-difference curl of E at
    the top face (k=0).

    Parameters
    ----------
    E_xpol, E_ypol : ndarray, shape (n_edges,)
    grid : Grid3D
    omega : float
    n_ex, n_ey, n_ez : int   Edge count offsets.

    Returns
    -------
    zxx, zxy, zyx, zyy : ndarray of complex, shape (n_stations,)
    """
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    nx, ny, nz = grid.nx, grid.ny, grid.nz
    ns = grid.n_stations

    # ── Reshape E into component arrays ──────────────────────────────────────
    def _ex_arr(E):
        return E[:n_ex].reshape(nx, ny + 1, nz + 1)

    def _ey_arr(E):
        return E[n_ex : n_ex + n_ey].reshape(nx + 1, ny, nz + 1)

    def _ez_arr(E):
        return E[n_ex + n_ey :].reshape(nx + 1, ny + 1, nz)

    # ── Surface E-field at nodes (i,j,k=0) ───────────────────────────────────
    # Ex[i,j,0]: i=0..nx-1, j=0..ny (top-face Ex edges)
    # Ey[i,j,0]: i=0..nx, j=0..ny-1 (top-face Ey edges)

    Ex_xp = _ex_arr(E_xpol)[:, :, 0]  # (nx, ny+1)
    Ey_xp = _ey_arr(E_xpol)[:, :, 0]  # (nx+1, ny)
    Ez_xp = _ez_arr(E_xpol)[:, :, 0]  # (nx+1, ny+1) ... first depth layer

    Ex_yp = _ex_arr(E_ypol)[:, :, 0]
    Ey_yp = _ey_arr(E_ypol)[:, :, 0]
    Ez_yp = _ez_arr(E_ypol)[:, :, 0]

    Ex_xp_next = _ex_arr(E_xpol)[:, :, 1]  # k=1 layer for ∂/∂z
    Ey_xp_next = _ey_arr(E_xpol)[:, :, 1]
    Ex_yp_next = _ex_arr(E_ypol)[:, :, 1]
    Ey_yp_next = _ey_arr(E_ypol)[:, :, 1]

    # ── Surface H at face (i,j,k=0) via curl of E ────────────────────────────
    # From curl E = −iωμ₀ H:
    #   H_x =  (∂Ey/∂z − ∂Ez/∂y) / (iωμ₀)
    #   H_y =  (∂Ez/∂x − ∂Ex/∂z) / (iωμ₀)

    # Hx: shape (nx+1, ny) at face (i,j,k=0) — i=0..nx, j=0..ny-1
    def _Hx(Ex_a, Ey_a, Ez_a, Ey_next):
        dEy_dz = (Ey_next - Ey_a) / dz[0]  # (nx+1, ny)
        dEz_dy = (Ez_a[:, 1:] - Ez_a[:, :-1]) / dy[None, :]  # (nx+1, ny)
        return (dEy_dz - dEz_dy) / (1j * omega * MU0)

    # Hy face (i,j,0): i=0..nx-1, j=0..ny
    def _Hy(Ex_a, Ez_a, Ex_next):
        dEz_dx = (Ez_a[1:, :] - Ez_a[:-1, :]) / dx[:, None]  # (nx, ny+1)
        dEx_dz = (Ex_next - Ex_a) / dz[0]  # (nx, ny+1)
        return (dEz_dx - dEx_dz) / (1j * omega * MU0)

    Hx_xp = _Hx(Ex_xp, Ey_xp, Ez_xp, Ey_xp_next)  # (nx+1, ny)
    Hy_xp = _Hy(Ex_xp, Ez_xp, Ex_xp_next)  # (nx, ny+1)
    Hx_yp = _Hx(Ex_yp, Ey_yp, Ez_yp, Ey_yp_next)
    Hy_yp = _Hy(Ex_yp, Ez_yp, Ex_yp_next)

    # ── Interpolate to station positions ─────────────────────────────────────
    x_nodes = grid.x_nodes  # (nx+1,)
    y_nodes = grid.y_nodes  # (ny+1,)
    x_st, y_st = grid.stations_xy[:, 0], grid.stations_xy[:, 1]

    def _interp2(field, x_vals, y_vals):
        """Bilinear interpolation of 2-D field to (x_st, y_st)."""
        xi = np.clip(
            np.searchsorted(x_vals, x_st, side="right") - 1,
            0,
            len(x_vals) - 2,
        )
        yi = np.clip(
            np.searchsorted(y_vals, y_st, side="right") - 1,
            0,
            len(y_vals) - 2,
        )
        tx = (x_st - x_vals[xi]) / (x_vals[xi + 1] - x_vals[xi])
        ty = (y_st - y_vals[yi]) / (y_vals[yi + 1] - y_vals[yi])
        return (
            (1 - tx) * (1 - ty) * field[xi, yi]
            + tx * (1 - ty) * field[xi + 1, yi]
            + (1 - tx) * ty * field[xi, yi + 1]
            + tx * ty * field[xi + 1, yi + 1]
        )

    x_ex = x_nodes[:-1]  # Ex defined at x_centers[i=0..nx-1]
    x_ey = x_nodes  # Ey at x_nodes[i=0..nx]
    y_ex = y_nodes  # Ex at y_nodes[j=0..ny]
    y_ey = y_nodes[:-1]  # Ey at y_centers[j=0..ny-1]

    x_hx = x_nodes  # Hx face at x_nodes[i=0..nx]
    y_hx = y_nodes[:-1]  # Hx face at y_centers
    x_hy = x_nodes[:-1]  # Hy face at x_centers
    y_hy = y_nodes  # Hy face at y_nodes

    Ex_xp_st = _interp2(Ex_xp, x_ex, y_ex)
    Ey_xp_st = _interp2(Ey_xp, x_ey, y_ey)
    Hx_xp_st = _interp2(Hx_xp, x_hx, y_hx)
    Hy_xp_st = _interp2(Hy_xp, x_hy, y_hy)

    Ex_yp_st = _interp2(Ex_yp, x_ex, y_ex)
    Ey_yp_st = _interp2(Ey_yp, x_ey, y_ey)
    Hx_yp_st = _interp2(Hx_yp, x_hx, y_hx)
    Hy_yp_st = _interp2(Hy_yp, x_hy, y_hy)

    # ── Solve for Z at each station: E_mat = Z × H_mat ───────────────────────
    zxx = np.zeros(ns, dtype=complex)
    zxy = np.zeros(ns, dtype=complex)
    zyx = np.zeros(ns, dtype=complex)
    zyy = np.zeros(ns, dtype=complex)

    for s in range(ns):
        E_mat = np.array(
            [[Ex_xp_st[s], Ex_yp_st[s]], [Ey_xp_st[s], Ey_yp_st[s]]],
            dtype=complex,
        )
        H_mat = np.array(
            [[Hx_xp_st[s], Hx_yp_st[s]], [Hy_xp_st[s], Hy_yp_st[s]]],
            dtype=complex,
        )
        det = H_mat[0, 0] * H_mat[1, 1] - H_mat[0, 1] * H_mat[1, 0]
        if abs(det) > 1e-30:
            Z = E_mat @ np.linalg.inv(H_mat)
            zxx[s], zxy[s] = Z[0, 0], Z[0, 1]
            zyx[s], zyy[s] = Z[1, 0], Z[1, 1]

    return zxx, zxy, zyx, zyy

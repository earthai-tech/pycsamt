# -*- coding: utf-8 -*-
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

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from .em2d import MT2DForward, _z_to_rho_phase
from .grid3d import Grid3D
from ..forward.batch import SurveyDataset3D

__all__ = [
    "MT3DForward",
    "ForwardResponse3D",
]

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

    freqs:       np.ndarray
    stations_xy: np.ndarray

    zxy:  np.ndarray
    zyx:  np.ndarray
    zxx:  np.ndarray
    zyy:  np.ndarray

    rho_a_xy:  np.ndarray
    phase_xy:  np.ndarray
    rho_a_yx:  np.ndarray
    phase_yx:  np.ndarray
    rho_a_xx:  np.ndarray
    phase_xx:  np.ndarray
    rho_a_yy:  np.ndarray
    phase_yy:  np.ndarray

    method: str   = "quasi3d"
    grid:   Grid3D = field(repr=False, default=None)

    # ── convenience ──────────────────────────────────────────────────────────

    @property
    def n_freqs(self) -> int:    return len(self.freqs)
    @property
    def n_stations(self) -> int: return len(self.stations_xy)
    @property
    def periods(self) -> np.ndarray: return 1.0 / self.freqs

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
            comp_list = [c.strip() for c in components.replace(",", "_").split("_") if c.strip()]

        parts = []
        for c in comp_list:
            rho  = getattr(self, f"rho_a_{c}")   # (n_freqs, n_stations)
            phi  = getattr(self, f"phase_{c}")
            if log_rho:
                rho = np.log10(np.maximum(rho, 1e-12))
            parts.append(rho.T)          # (n_stations, n_freqs)
            if include_phase:
                parts.append(phi.T)
        return np.concatenate(parts, axis=1)

    def to_survey_dataset(
        self,
        y_models: Optional[np.ndarray] = None,
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
        X = self.to_feature_array(components=components)       # (n_st, n_feat)
        if y_models is None:
            y_models = np.zeros((self.n_stations, 1), dtype=np.float32)
        X = X[None].astype(np.float32)       # (1, n_st, n_feat)
        y = y_models[None].astype(np.float32)  # (1, n_st, n_params)
        return SurveyDataset3D(
            X=X, y=y,
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
        self.freqs   = np.asarray(freqs, dtype=float)
        self.grid    = grid
        self.method  = method.lower()
        self.verbose = verbose

        if self.method != "quasi3d":
            raise NotImplementedError(
                f"method={method!r} is not yet implemented. "
                "Use method='quasi3d'."
            )

    def run(self) -> ForwardResponse3D:
        """Run the forward solver and return the full impedance tensor.

        Returns
        -------
        ForwardResponse3D
        """
        return self._run_quasi3d()

    # ─────────────────────────────────────────────────────────────────────────

    def _run_quasi3d(self) -> ForwardResponse3D:
        """Quasi-3D solver: orthogonal 2-D profile stacking."""
        grid   = self.grid
        freqs  = self.freqs
        ns     = grid.n_stations
        nf     = len(freqs)

        # Accumulation buffers for xz and yz profile contributions
        zxy_xz = np.zeros((nf, ns), dtype=complex)
        zyx_xz = np.zeros((nf, ns), dtype=complex)
        zxy_yz = np.zeros((nf, ns), dtype=complex)
        zyx_yz = np.zeros((nf, ns), dtype=complex)

        # ── XZ profiles  (one per unique y-row of stations) ──────────────────
        y_cells         = grid._station_y_cells()
        unique_y_cells  = np.unique(y_cells)
        n_xz_profiles   = len(unique_y_cells)

        if self.verbose:
            print(f"  [MT3D quasi3d] {n_xz_profiles} XZ profiles + "
                  f"{len(np.unique(grid._station_x_cells()))} YZ profiles …")

        for ki, yi in enumerate(unique_y_cells):
            if self.verbose:
                print(f"    XZ profile {ki+1}/{n_xz_profiles}  (y-cell {yi})",
                      end="\r", flush=True)

            g2d_xz, st_idx = grid.xz_slice(yi)
            resp_xz        = MT2DForward(freqs, g2d_xz, verbose=False).run()

            for local_i, global_i in enumerate(st_idx):
                zxy_xz[:, global_i] = resp_xz.zxy[:, local_i]
                zyx_xz[:, global_i] = resp_xz.zyx[:, local_i]

        # ── YZ profiles  (one per unique x-column of stations) ───────────────
        x_cells         = grid._station_x_cells()
        unique_x_cells  = np.unique(x_cells)
        n_yz_profiles   = len(unique_x_cells)

        for ki, xi in enumerate(unique_x_cells):
            if self.verbose:
                print(f"    YZ profile {ki+1}/{n_yz_profiles}  (x-cell {xi})",
                      end="\r", flush=True)

            g2d_yz, st_idx = grid.yz_slice(xi)
            resp_yz        = MT2DForward(freqs, g2d_yz, verbose=False).run()

            # In the YZ solver, horizontal axis = y_3D.
            # TE  → E-field perpendicular to yz-plane = E_x_3D
            #     → zxy_yz_solver = E_x_3D / H_y_3D = Z_xy_3D contribution
            # TM  → H-field perpendicular to yz-plane = H_x_3D
            #     → zyx_yz_solver = E_y_3D / H_x_3D = Z_yx_3D contribution
            for local_i, global_i in enumerate(st_idx):
                zxy_yz[:, global_i] = resp_yz.zxy[:, local_i]
                zyx_yz[:, global_i] = resp_yz.zyx[:, local_i]

        if self.verbose:
            print(f"  [MT3D quasi3d] done.                          ")

        # ── Assemble quasi-3D tensor ─────────────────────────────────────────
        # Average contributions from both profile directions
        zxy = 0.5 * (zxy_xz + zxy_yz)
        zyx = 0.5 * (zyx_xz + zyx_yz)

        # Off-diagonal: zero in the 2-D approximation
        zxx = np.zeros((nf, ns), dtype=complex)
        zyy = np.zeros((nf, ns), dtype=complex)

        omega = 2.0 * np.pi * freqs[:, None]   # (nf, 1) for broadcasting

        rho_a_xy, phase_xy = _z_to_rho_phase(zxy, omega)
        rho_a_yx, phase_yx = _z_to_rho_phase(zyx, omega)
        rho_a_xx, phase_xx = _z_to_rho_phase(zxx, omega)
        rho_a_yy, phase_yy = _z_to_rho_phase(zyy, omega)

        return ForwardResponse3D(
            freqs=freqs,
            stations_xy=grid.stations_xy.copy(),
            zxy=zxy, zyx=zyx, zxx=zxx, zyy=zyy,
            rho_a_xy=rho_a_xy, phase_xy=phase_xy,
            rho_a_yx=rho_a_yx, phase_yx=phase_yx,
            rho_a_xx=rho_a_xx, phase_xx=phase_xx,
            rho_a_yy=rho_a_yy, phase_yy=phase_yy,
            method=self.method,
            grid=grid,
        )

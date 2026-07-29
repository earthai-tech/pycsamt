# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Physics-informed AI 3-D MT inversion.

:class:`PINNInverter3D` optimises a quasi-3D
resistivity volume — one 1-D column per station
arranged in a 2-D horizontal network — by
minimising a physics-informed loss that combines
per-station 1-D data misfits with a graph-Laplacian
spatial smoothness penalty.

The station network is represented as a graph whose
adjacency ``A`` encodes inter-station proximity.
Spatial coupling is imposed via the graph Laplacian
:math:`L = D - A`:

.. math::

    \mathcal{L}_{\rm spatial} =
    \frac{1}{\|m\|}
    \sum_{k=1}^{L}
    m_{:,k}^\top L\, m_{:,k}

where :math:`m_{s,k} = \log_{10}\rho_{s,k}`,
:math:`L` is the graph Laplacian built from the
normalized adjacency, and the sum is over layers.

Full loss
---------

.. math::

    \mathcal{L} =
    \mathcal{L}_{\rm data}
    + \lambda_z \mathcal{L}_{\rm vert}
    + \lambda_g \mathcal{L}_{\rm spatial}

where :math:`\mathcal{L}_{\rm data}` and
:math:`\mathcal{L}_{\rm vert}` are the same
data-misfit and vertical-smoothness terms used in
:class:`~pycsamt.ai.inversion.pinn2d.PINNInverter2D`.

Example
-------
>>> from pycsamt.ai.inversion import PINNInverter3D
>>> inv = PINNInverter3D(  # doctest: +SKIP
...     "edi/survey/",
...     n_layers=10,
...     depth_max=3000.0,
...     radius=3000.0,
...     epochs=300,
... )
>>> inv.fit()  # doctest: +SKIP
PINNInverter3D(n_stations=25, fitted)
>>> vol = inv.resistivity_volume()  # doctest: +SKIP
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .._base import BasePINNInverter
from ._pinn_ops import fit_3d_joint
from ._sites_bridge import (
    SiteObs2D,
    _interp_to_grid,
    _make_common_grid,
    sites_to_coords_3d,
    sites_to_obs_2d,
)

__all__ = ["PINNInverter3D"]


# ── backward-compat re-export (hybrid3d imports this)


class PINNInverter3D(BasePINNInverter):
    r"""
    Physics-informed quasi-3D MT inversion.

    Optimises a per-station 1-D column model for all
    stations simultaneously using the 1-D MT forward
    with graph-Laplacian spatial smoothness coupling
    neighbour stations through their shared adjacency.

    Parameters
    ----------
    sites : Any
        Path, ``EDIFile``, ``EDICollection``, ``Site``,
        ``Sites``, ``APISurvey``, or iterable.
    n_layers : int, default 10
        Number of model layers per station.
    depth_max : float, default 2000.0
        Target maximum depth in metres.
    n_freqs : int, default 32
        Frequency-grid points for the common grid.
    mode : {'te', 'tm', 'both'}, default 'te'
        Which polarisation to fit.
    smoothness_weight : float, default 0.01
        Vertical smoothness weight.
    graph_weight : float, default 0.005
        Graph-Laplacian spatial smoothness weight.
    radius : float, default 5000.0
        Edge radius in metres used to build the
        station adjacency when *adjacency* is None.
    adjacency : ndarray (S, S) or None
        Pre-computed normalised adjacency matrix.
        When ``None`` it is built from station
        positions with :func:`build_adjacency`.
    station_coords : ndarray (S, 2) or None
        Explicit ``(x, y)`` station positions [m].
        Auto-extracted from site metadata if None.
    station_spacing : float, default 500.0
        Uniform grid spacing [m] used when no
        geographic coordinates are available.
    epochs : int, default 300
        Adam iterations.
    lr : float, default 1e-2
        Adam learning rate.
    comp_te : str, default ``'xy'``
    comp_tm : str, default ``'yx'``
    device : str or None
    recursive : bool, default True
    on_dup : str, default ``'replace'``
    verbose : int, default 0
    """

    def __init__(
        self,
        sites: Any,
        *,
        n_layers: int = 10,
        depth_max: float = 2000.0,
        n_freqs: int = 32,
        mode: str = "te",
        smoothness_weight: float = 0.01,
        graph_weight: float = 0.005,
        radius: float = 5000.0,
        adjacency: np.ndarray | None = None,
        station_coords: np.ndarray | None = None,
        station_spacing: float = 500.0,
        epochs: int = 300,
        lr: float = 1e-2,
        comp_te: str = "xy",
        comp_tm: str = "yx",
        device: str | None = None,
        recursive: bool = True,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> None:
        if mode not in ("te", "tm", "both"):
            raise ValueError(f"mode must be 'te', 'tm', or 'both'; got {mode!r}.")
        super().__init__(
            n_layers=n_layers,
            depth_max=depth_max,
            device=device,
        )
        self.n_freqs = int(n_freqs)
        self.mode = mode
        self.smoothness_weight = float(smoothness_weight)
        self.graph_weight = float(graph_weight)
        self.radius = float(radius)
        self.epochs = int(epochs)
        self.lr = float(lr)
        self.comp_te = comp_te
        self.comp_tm = comp_tm
        self.verbose = verbose

        self._obs: list[SiteObs2D] = sites_to_obs_2d(
            sites,
            comp_te=comp_te,
            comp_tm=comp_tm,
            recursive=recursive,
            on_dup=on_dup,
            verbose=verbose,
        )
        self._freqs_grid = _make_common_grid(
            [o.freq for o in self._obs],
            n_freqs=n_freqs,
        )

        # Station coordinates
        if station_coords is not None:
            self._coords = np.asarray(station_coords, dtype=float)
        else:
            self._coords = sites_to_coords_3d(
                sites,
                station_spacing=station_spacing,
                recursive=recursive,
                on_dup=on_dup,
                verbose=verbose,
            )

        # Adjacency matrix
        if adjacency is not None:
            self._adjacency = np.asarray(adjacency, dtype=np.float64)
        else:
            from pycsamt.ai.nets.gcn import (
                build_adjacency,
            )

            self._adjacency = build_adjacency(
                self._coords,
                radius=self.radius,
            ).astype(np.float64)

        self._result: dict | None = None

    # ── fit ──

    def fit(
        self,
        *,
        verbose: bool = True,
        log_every: int = 50,
    ) -> PINNInverter3D:
        """
        Run the joint 3-D physics-informed inversion.

        Returns
        -------
        self
        """
        self._require_backend()
        dev = self._resolve_device()
        lr_obs, ph_obs = self._build_obs_arrays()

        if verbose:
            S = len(self._obs)
            print(
                f"PINNInverter3D: optimising "
                f"{S} stations x {self.n_layers} "
                f"layers ({self.epochs} epochs) ..."
            )

        _every = log_every if verbose else 0
        self._result = fit_3d_joint(
            lr_obs,
            ph_obs,
            self._freqs_grid,
            self._adjacency,
            n_layers=self.n_layers,
            depth_max=self.depth_max,
            lam_z=self.smoothness_weight,
            lam_g=self.graph_weight,
            lr=self.lr,
            epochs=self.epochs,
            device=dev,
            log_every=_every,
        )
        self._history = self._result["history"]
        self._is_fitted = True
        return self

    # ── outputs ──

    def resistivity_volume(self, *, as_log10: bool = True) -> np.ndarray:
        """
        Return the quasi-3D resistivity volume.

        Returns
        -------
        ndarray (n_layers, n_stations)
            Column ``i`` is the 1-D model at
            station ``i`` in the same order as
            :attr:`stations`.
        """
        self._check_fitted()
        lr = self._result["log_rho"]  # (S, L)
        section = lr.T  # (L, S)
        if as_log10:
            return section
        return 10.0**section

    def thickness_volume(self) -> np.ndarray:
        """
        Return layer thicknesses in metres.

        Returns
        -------
        ndarray (n_layers-1, n_stations)
        """
        self._check_fitted()
        lt = self._result["log_thick"]  # (S, L-1)
        return (10.0**lt).T  # (L-1, S)

    def station_coords(self) -> np.ndarray:
        """Return station (x, y) positions [m]."""
        return self._coords.copy()

    def adjacency(self) -> np.ndarray:
        """Return the station adjacency matrix."""
        return self._adjacency.copy()

    def convergence_curve(self):
        """
        Return Adam loss history.

        Returns
        -------
        pandas.DataFrame
            Columns: epoch, loss.
        """
        self._check_fitted()
        import pandas as pd

        return pd.DataFrame(
            {
                "epoch": range(
                    1,
                    len(self._result["history"]) + 1,
                ),
                "loss": self._result["history"],
            }
        )

    def residuals(self):
        """
        Observed vs predicted data fit.

        Returns
        -------
        pandas.DataFrame
            Columns: station, freq, rho_obs,
            rho_pred, phase_obs, phase_pred.
        """
        self._check_fitted()
        import pandas as pd

        from pycsamt.forward.em1d import MT1DForward
        from pycsamt.forward.synthetic import LayeredModel

        lr_2d = self._result["log_rho"]
        lt_2d = self._result["log_thick"]
        rows = []
        for i, obs in enumerate(self._obs):
            try:
                m = LayeredModel(
                    resistivity=np.maximum(10.0 ** lr_2d[i], 1e-3),
                    thickness=np.maximum(10.0 ** lt_2d[i], 1.0),
                )
                resp = MT1DForward(obs.freq).run(m)
                rp, pp = resp.rho_a, resp.phase
            except Exception:
                rp = np.full_like(obs.freq, np.nan)
                pp = np.full_like(obs.freq, np.nan)

            rho_obs = obs.rho_te if self.mode in ("te", "both") else obs.rho_tm
            ph_obs = obs.phase_te if self.mode in ("te", "both") else obs.phase_tm
            for k in range(len(obs.freq)):
                rows.append(
                    {
                        "station": obs.name,
                        "freq": obs.freq[k],
                        "rho_obs": rho_obs[k],
                        "rho_pred": rp[k],
                        "phase_obs": ph_obs[k],
                        "phase_pred": pp[k],
                    }
                )
        return pd.DataFrame(rows)

    # ── read-only properties ──

    @property
    def stations(self) -> list[str]:
        """Station names in order."""
        return [o.name for o in self._obs]

    @property
    def n_sites(self) -> int:
        """Number of loaded stations."""
        return len(self._obs)

    # ── internals ──

    def _build_obs_arrays(
        self,
    ) -> tuple[np.ndarray, np.ndarray]:
        S = len(self._obs)
        F = len(self._freqs_grid)
        lr_obs = np.full((S, F), np.nan)
        ph_obs = np.full((S, F), np.nan)

        for i, o in enumerate(self._obs):
            if self.mode == "te":
                rho_src = o.rho_te
                ph_src = o.phase_te
            elif self.mode == "tm":
                rho_src = o.rho_tm
                ph_src = o.phase_tm
            else:
                rho_src = 0.5 * (o.rho_te + o.rho_tm)
                ph_src = 0.5 * (o.phase_te + o.phase_tm)

            lg, pg = _interp_to_grid(
                o.freq,
                rho_src,
                ph_src,
                self._freqs_grid,
            )
            lr_obs[i] = lg
            ph_obs[i] = pg

        return lr_obs, ph_obs

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"PINNInverter3D("
            f"n_stations={self.n_sites}, "
            f"n_layers={self.n_layers}, "
            f"{status})"
        )

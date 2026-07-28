# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Hybrid AI + physics quasi-3D MT inversion.

:class:`HybridInverter3D` combines the spatially
coherent output of a pre-trained
:class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`
with the physics-informed refinement of
:class:`~pycsamt.ai.inversion.pinn3d.PINNInverter3D`.

**Stage 1 — GCN initial model**

The graph-convolutional network produces per-station
1-D models by message-passing over the station
network::

    X, A --> GCNInverter3D --> params_0

where ``params_0`` has shape
``(n_stations, 2*n_layers-1)`` in log10 scale.

**Stage 2 — Physics refinement**

Starting from ``params_0``, Adam minimises the same
loss used by
:class:`~pycsamt.ai.inversion.pinn3d.PINNInverter3D`:

.. math::

    \mathcal{L} =
    \mathcal{L}_{\rm data}
    + \lambda_z \mathcal{L}_{\rm vert}
    + \lambda_g \mathcal{L}_{\rm spatial}

Example
-------
>>> from pycsamt.ai.inversion import (
...     GCNInverter3D,
...     HybridInverter3D,
... )
>>> gcn = GCNInverter3D.load(  # doctest: +SKIP
...     "checkpoints/gcn3d.npz"
... )
>>> inv = HybridInverter3D(  # doctest: +SKIP
...     "edi/survey/",
...     ai_inverter=gcn,
...     epochs=150,
...     graph_weight=0.003,
... )
>>> inv.fit()  # doctest: +SKIP
HybridInverter3D(n_stations=25, fitted)
>>> vol = inv.resistivity_volume()  # doctest: +SKIP
>>> s1 = inv.stage1_volume()  # doctest: +SKIP
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .._base import BaseHybridInverter
from ._pinn_ops import fit_3d_joint
from ._sites_bridge import (
    SiteObs2D,
    _interp_to_grid,
    _make_common_grid,
    obs_to_features_1d,
    sites_to_coords_3d,
    sites_to_obs_2d,
)

__all__ = ["HybridInverter3D"]


# ── HybridInverter3D ──


class HybridInverter3D(BaseHybridInverter):
    r"""
    Two-stage hybrid AI + physics quasi-3D inversion.

    Parameters
    ----------
    sites : Any
        Path, ``EDIFile``, ``EDICollection``, ``Site``,
        ``Sites``, ``APISurvey``, or iterable.
    ai_inverter : GCNInverter3D or str or Path
        Pre-trained (fitted)
        :class:`~pycsamt.ai.inversion.inv3d.GCNInverter3D`
        or path to a ``.npz`` checkpoint.
    n_layers : int or None
        Layers per station for Stage 2.  Defaults to
        ``ai_inverter.n_layers``.
    depth_max : float, default 2000.0
        Target depth for uniform thickness init.
    n_freqs : int, default 32
        Frequency-grid size for Stage 2 optimisation.
    mode : {'te', 'tm', 'both'}, default 'te'
        Polarisation used in Stage 2.
    smoothness_weight : float, default 0.005
        Vertical smoothness weight for Stage 2.
    graph_weight : float, default 0.003
        Graph spatial smoothness weight.
    radius : float, default 5000.0
        Edge radius [m] for adjacency construction.
    adjacency : ndarray (S, S) or None
        Pre-computed adjacency.  Built from station
        positions if ``None``.
    station_coords : ndarray (S, 2) or None
        Explicit ``(x, y)`` positions [m].
    station_spacing : float, default 500.0
        Fallback uniform grid spacing [m].
    epochs : int, default 150
        Adam iterations for Stage 2.
    lr : float, default 5e-3
        Adam learning rate for Stage 2.
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
        ai_inverter: Any | str | Path,
        *,
        n_layers: int | None = None,
        depth_max: float = 2000.0,
        n_freqs: int = 32,
        mode: str = "te",
        smoothness_weight: float = 0.005,
        graph_weight: float = 0.003,
        radius: float = 5000.0,
        adjacency: np.ndarray | None = None,
        station_coords: np.ndarray | None = None,
        station_spacing: float = 500.0,
        epochs: int = 150,
        lr: float = 5e-3,
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

        self._ai_inv = self._load_ai_inverter(ai_inverter)
        self.n_layers = int(n_layers if n_layers is not None else self._ai_inv.n_layers)

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

        # Station coordinates and adjacency
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

        self._stage1_params: np.ndarray | None = None
        self._result: dict | None = None

    # ── fit ──

    def fit(
        self,
        *,
        verbose: bool = True,
        log_every: int = 50,
    ) -> HybridInverter3D:
        """
        Run both inversion stages.

        Returns
        -------
        self
        """
        self._require_backend()
        dev = self._resolve_device()

        # Stage 1: GCN prediction
        if verbose:
            print("Stage 1: GCNInverter3D forward pass ...")
        init_params = self._run_stage1()

        # init_params: (S, 2*n_layers-1)
        S = len(self._obs)
        n_l = self.n_layers
        init_log_rho = init_params[:, :n_l]
        init_log_thick = init_params[:, n_l:]

        # Uniform thickness fallback if shape mismatch
        if init_log_thick.shape[1] != n_l - 1:
            dz = max(self.depth_max / n_l, 1.0)
            init_log_thick = np.full((S, n_l - 1), np.log10(dz))

        # Build observed arrays for Stage 2
        lr_obs, ph_obs = self._build_obs_arrays()

        if verbose:
            print(
                f"Stage 2: joint physics refinement "
                f"({S} stations, {n_l} "
                f"layers, {self.epochs} epochs) ..."
            )

        _every = log_every if verbose else 0
        self._result = fit_3d_joint(
            lr_obs,
            ph_obs,
            self._freqs_grid,
            self._adjacency,
            n_layers=n_l,
            depth_max=self.depth_max,
            lam_z=self.smoothness_weight,
            lam_g=self.graph_weight,
            lr=self.lr,
            epochs=self.epochs,
            device=dev,
            log_every=_every,
            init_log_rho=init_log_rho,
            init_log_thick=init_log_thick,
        )
        self._history = self._result["history"]
        self._is_fitted = True
        return self

    # ── outputs ──

    def resistivity_volume(self, *, as_log10: bool = True) -> np.ndarray:
        """
        Return the Stage-2 quasi-3D volume.

        Returns
        -------
        ndarray (n_layers, n_stations)
        """
        self._check_fitted()
        lr = self._result["log_rho"]  # (S, L)
        s = lr.T  # (L, S)
        return s if as_log10 else 10.0**s

    def thickness_volume(self) -> np.ndarray:
        """
        Return Stage-2 layer thicknesses [m].

        Returns
        -------
        ndarray (n_layers-1, n_stations)
        """
        self._check_fitted()
        lt = self._result["log_thick"]
        return (10.0**lt).T

    def stage1_volume(self, *, as_log10: bool = True) -> np.ndarray:
        """
        Return the Stage-1 GCN quasi-3D volume.

        Returns
        -------
        ndarray (n_layers, n_stations)
        """
        if self._stage1_params is None:
            raise RuntimeError("Call fit() to populate Stage-1 results.")
        n_l = self.n_layers
        lr = self._stage1_params[:, :n_l].T  # (L,S)
        return lr if as_log10 else 10.0**lr

    def convergence_curve(self):
        """
        Return Stage-2 Adam loss history.

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

    def residuals(self, stage: int = 2):
        """
        Observed vs predicted data fit.

        Parameters
        ----------
        stage : {1, 2}, default 2

        Returns
        -------
        pandas.DataFrame
            Columns: station, stage, freq,
            rho_obs, rho_pred, phase_obs,
            phase_pred.
        """
        self._check_fitted()
        import pandas as pd

        from pycsamt.forward.em1d import MT1DForward
        from pycsamt.forward.synthetic import LayeredModel

        n_l = self.n_layers
        if stage == 1:
            if self._stage1_params is None:
                raise RuntimeError("Call fit() first.")
            lr_2d = self._stage1_params[:, :n_l]
            dz = max(self.depth_max / n_l, 1.0)
            lt_2d = np.full(
                (len(self._obs), n_l - 1),
                np.log10(dz),
            )
        else:
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
                        "stage": stage,
                        "freq": obs.freq[k],
                        "rho_obs": rho_obs[k],
                        "rho_pred": rp[k],
                        "phase_obs": ph_obs[k],
                        "phase_pred": pp[k],
                    }
                )
        return pd.DataFrame(rows)

    def station_coords(self) -> np.ndarray:
        """Return station (x, y) positions [m]."""
        return self._coords.copy()

    def adjacency(self) -> np.ndarray:
        """Return the station adjacency matrix."""
        return self._adjacency.copy()

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

    def _load_ai_inverter(self, ai_inverter):
        from .inv3d import GCNInverter3D

        if isinstance(ai_inverter, GCNInverter3D):
            if not ai_inverter._is_fitted:
                raise ValueError("ai_inverter must be a fitted GCNInverter3D instance.")
            return ai_inverter
        if isinstance(ai_inverter, (str, Path)):
            return GCNInverter3D.load(Path(ai_inverter))
        raise TypeError(
            "ai_inverter must be a fitted "
            "GCNInverter3D or a path to a checkpoint;"
            f" got {type(ai_inverter)!r}."
        )

    def _run_stage1(self) -> np.ndarray:
        """
        Apply GCNInverter3D to extract initial models.

        Returns
        -------
        params : ndarray (S, 2*n_layers-1)
            Per-station model parameters in log10
            scale.
        """
        ai = self._ai_inv
        # n_features = 2 * n_freqs for 1D features
        gcn_nfreqs = max(ai.n_features // 2, 1)

        # Use obs_to_features_1d — safe with
        # SiteObs2D lists (avoids ensure_sites).
        X, _, _ = obs_to_features_1d(self._obs, n_freqs=gcn_nfreqs)
        # X: (S, n_features)

        # Build adjacency for the GCN call
        A = self._adjacency.astype(np.float32)
        # Use training adjacency if S mismatch;
        # otherwise pass the computed one.
        params = ai.predict(
            X,
            adjacency=A,
            as_log_rho=True,
        )
        # params: (S, n_out) where n_out=2*n_layers-1
        self._stage1_params = params
        return params

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
            f"HybridInverter3D("
            f"n_stations={self.n_sites}, "
            f"n_layers={self.n_layers}, "
            f"{status})"
        )

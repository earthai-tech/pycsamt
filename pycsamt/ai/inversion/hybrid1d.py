# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Hybrid AI + physics 1-D MT/CSAMT inversion.

:class:`HybridInverter1D` combines a pre-trained
:class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`
(Option 1, supervised AI) with a physics-informed
refinement step (Option 2 gradient descent).

**Stage 1 — AI initial model**

The supervised inverter maps the observed EM response
directly to a layered Earth model::

    d_obs --> EMInverter1D --> m_0  (fast, per-station)

**Stage 2 — Physics refinement**

Starting from ``m_0``, Adam minimises the same
physics-informed loss as
:class:`~pycsamt.ai.inversion.pinn1d.PINNInverter1D`:

.. math::

    \mathcal{L}(\theta) =
    \frac{1}{F}\sum_{f}
    \left[
      \left(
        \log_{10}\frac{\rho_a^{\rm pred}}
                      {\rho_a^{\rm obs}}
      \right)^2
      + \left(
          \frac{\phi^{\rm pred} - \phi^{\rm obs}}{90}
        \right)^2
    \right]
    + \lambda\sum_k(\rho_k - \rho_{k-1})^2

Because ``m_0`` is physically plausible, Stage 2
converges faster and more reliably than a PINN started
from a naive initialisation.

Example
-------
>>> from pycsamt.ai.inversion import (
...     EMInverter1D,
...     HybridInverter1D,
... )
>>> ai = EMInverter1D.load(  # doctest: +SKIP
...     "checkpoints/mt1d_resnet.npz"
... )
>>> inv = HybridInverter1D(  # doctest: +SKIP
...     "edi/",
...     ai_inverter=ai,
...     max_iter=200,
...     smoothness_weight=0.005,
... )
>>> inv.fit()  # doctest: +SKIP
HybridInverter1D(n_stations=5, fitted)
>>> models = inv.predict()  # doctest: +SKIP
>>> df = inv.convergence_curves()  # doctest: +SKIP
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .._base import BaseHybridInverter
from ._pinn_ops import fit_station
from ._sites_bridge import (
    SiteObs1D,
    obs_to_features_1d,
    sites_to_obs_1d,
)

__all__ = ["HybridInverter1D"]


# ── HybridInverter1D ──


class HybridInverter1D(BaseHybridInverter):
    r"""
    Two-stage hybrid AI + physics 1-D inversion.

    Parameters
    ----------
    sites : Any
        Path, ``EDIFile``, ``EDICollection``, ``Site``,
        ``Sites``, ``APISurvey``, or iterable.
    ai_inverter : EMInverter1D or str or Path
        Pre-trained supervised inverter (fitted
        :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`)
        or path to a saved ``.npz`` checkpoint.
    solver : {'mt1d', 'csamt1d'}, default 'mt1d'
        EM physics used in the refinement step.
    max_iter : int, default 200
        Adam iterations for the physics refinement.
    smoothness_weight : float, default 0.005
        Regularisation weight on log-resistivity
        first differences.
    lr : float, default 5e-3
        Adam learning rate for Stage 2.
    device : str or None
        Torch device.  Auto-detects CUDA/CPU if None.
    comp : {'xy', 'yx', 'xx', 'yy'}, default 'xy'
        Impedance component for observed data.
    n_freqs : int, default 32
        Frequency-grid size for EMInverter1D input.
    recursive : bool, default True
    on_dup : str, default 'replace'
    verbose : int, default 0
    """

    def __init__(
        self,
        sites: Any,
        ai_inverter: Any | str | Path,
        *,
        solver: str = "mt1d",
        max_iter: int = 200,
        smoothness_weight: float = 0.005,
        lr: float = 5e-3,
        device: str | None = None,
        comp: str = "xy",
        n_freqs: int = 32,
        recursive: bool = True,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> None:
        if solver not in ("mt1d", "csamt1d"):
            raise ValueError(
                "solver must be 'mt1d' or 'csamt1d'; "
                f"got {solver!r}.  TEM not supported."
            )
        super().__init__(device=device)
        self.solver = solver
        self.max_iter = int(max_iter)
        self.smoothness_weight = float(smoothness_weight)
        self.lr = float(lr)
        self.comp = comp
        self.n_freqs = int(n_freqs)
        self.recursive = recursive
        self.on_dup = on_dup
        self.verbose = verbose

        self._ai_inv = self._load_ai_inverter(ai_inverter)
        self._obs: list[SiteObs1D] = sites_to_obs_1d(
            sites,
            comp=comp,
            recursive=recursive,
            on_dup=on_dup,
            strict=True,
            verbose=verbose,
        )
        self._stage1: list[dict] = []
        self._stage2: list[dict] = []

    # ── fit ──

    def fit(
        self,
        *,
        verbose: bool = True,
        log_every: int = 50,
    ) -> HybridInverter1D:
        """
        Run both inversion stages.

        Parameters
        ----------
        verbose : bool, default True
            Print per-station progress.
        log_every : int, default 50
            Epoch-detail print frequency for Stage 2.

        Returns
        -------
        self
        """
        from pycsamt.api.view.progress import get_progress_bar

        self._require_backend()
        dev = self._resolve_device()

        # Stage 1: AI initial model
        init_params = self._run_stage1(verbose=verbose)

        # Stage 2: physics refinement
        _every = log_every if verbose else 0
        self._stage2 = []
        n_layers = self._ai_inv.n_layers

        with get_progress_bar(
            total=len(self._obs),
            desc="Stage 2 refinement",
            unit="station",
            verbose=verbose,
        ) as bar:
            for obs, params in zip(self._obs, init_params):
                # Unpack AI output: [log_rho..., log_thick...]
                init_lr = params[:n_layers]
                init_lt = params[n_layers:]

                res = fit_station(
                    obs,
                    n_layers=n_layers,
                    depth_max=0.0,
                    lam=self.smoothness_weight,
                    lr=self.lr,
                    epochs=self.max_iter,
                    device=dev,
                    log_every=_every,
                    init_log_rho=init_lr,
                    init_log_thick=init_lt,
                )
                self._stage2.append(res)
                bar.set_description(f"Stage 2 refinement: {obs.name}")
                bar.update(1)

        self._is_fitted = True
        return self

    # ── predict / diagnostics ──

    def predict(self) -> list:
        """
        Return Stage-2 refined layered models.

        Returns
        -------
        list of LayeredModel
        """
        self._check_fitted()
        return self._results_to_models(self._stage2)

    def stage1_models(self) -> list:
        """
        Return Stage-1 (AI-only) layered models.

        These are the AI starting points before
        physics refinement.

        Returns
        -------
        list of LayeredModel
        """
        if not self._stage1:
            raise RuntimeError("Call fit() to populate stage-1 models.")
        return self._results_to_models(self._stage1)

    def convergence_curves(self):
        """
        Return Stage-2 Adam loss history.

        Returns
        -------
        pandas.DataFrame
            Columns: station, epoch, loss.
        """
        self._check_fitted()
        import pandas as pd

        rows = []
        for obs, res in zip(self._obs, self._stage2):
            for ep, val in enumerate(res["history"], start=1):
                rows.append(
                    {
                        "station": obs.name,
                        "epoch": ep,
                        "loss": val,
                    }
                )
        return pd.DataFrame(rows)

    def residuals(self, stage: int = 2):
        """
        Observed vs predicted data fit.

        Parameters
        ----------
        stage : {1, 2}, default 2
            Which stage's models to evaluate.

        Returns
        -------
        pandas.DataFrame
            Columns: station, freq, rho_obs, rho_pred,
            phase_obs, phase_pred.
        """
        self._check_fitted()
        import pandas as pd

        from pycsamt.forward.em1d import (
            CSAMT1DForward,
            MT1DForward,
        )
        from pycsamt.forward.synthetic import LayeredModel

        results = self._stage2 if stage == 2 else self._stage1
        Fwd = MT1DForward if self.solver == "mt1d" else CSAMT1DForward
        rows = []
        for obs, res in zip(self._obs, results):
            rho = 10.0 ** res["log_rho"]
            thick = 10.0 ** res["log_thick"]
            try:
                m = LayeredModel(
                    resistivity=np.maximum(rho, 1e-3),
                    thickness=np.maximum(thick, 1.0),
                )
                resp = Fwd(obs.freq).run(m)
                rho_pred = resp.rho_a
                ph_pred = resp.phase
            except Exception:
                rho_pred = np.full_like(obs.rho_obs, np.nan)
                ph_pred = np.full_like(obs.phase_obs, np.nan)
            for k in range(len(obs.freq)):
                rows.append(
                    {
                        "station": obs.name,
                        "stage": stage,
                        "freq": obs.freq[k],
                        "rho_obs": obs.rho_obs[k],
                        "rho_pred": rho_pred[k],
                        "phase_obs": obs.phase_obs[k],
                        "phase_pred": ph_pred[k],
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

    def _load_ai_inverter(self, ai_inverter):
        """Load or validate the AI inverter."""
        from .inv1d import EMInverter1D

        if isinstance(ai_inverter, EMInverter1D):
            if not ai_inverter._is_fitted:
                raise ValueError(
                    "ai_inverter must be a fitted EMInverter1D instance."
                )
            return ai_inverter
        if isinstance(ai_inverter, (str, Path)):
            return EMInverter1D.load(Path(ai_inverter))
        raise TypeError(
            "ai_inverter must be a fitted "
            "EMInverter1D or a path to a checkpoint; "
            f"got {type(ai_inverter)!r}."
        )

    def _run_stage1(self, *, verbose: bool) -> list[np.ndarray]:
        """
        Apply EMInverter1D to all stations.

        Returns a list of parameter arrays with shape
        ``(2*n_layers-1,)`` per station (same format as
        :meth:`~EMInverter1D.predict` output rows).
        """
        n_layers = self._ai_inv.n_layers
        n_features = self._ai_inv._n_features
        if n_features is None:
            raise RuntimeError(
                "ai_inverter._n_features is None. "
                "Re-fit or reload the checkpoint."
            )
        # Derive n_freqs for the common grid
        # n_features = 2 * n_freqs (rho block + phase)
        _n_freqs = max(n_features // 2, 1)

        if verbose:
            print(
                "Stage 1: applying EMInverter1D to "
                f"{len(self._obs)} station(s) ..."
            )
        # Use obs_to_features_1d to avoid running
        # ensure_sites on SiteObs1D dataclasses.
        X, _, _ = obs_to_features_1d(self._obs, n_freqs=_n_freqs)
        # predict returns (n_stations, 2*n_layers-1)
        params = self._ai_inv.predict(X, as_log_rho=True)
        # Store as stage1 results (compatible format)
        self._stage1 = []
        for _i, row in enumerate(params):
            # row layout: [log_rho * n_layers,
            #              log_thick * (n_layers-1)]
            log_rho_i = row[:n_layers]
            log_thick_i = row[n_layers:]
            self._stage1.append(
                {
                    "log_rho": log_rho_i,
                    "log_thick": log_thick_i,
                    "history": [],
                }
            )
        return list(params)

    def _results_to_models(self, results: list[dict]) -> list:
        """Convert result dicts to LayeredModel list."""
        from pycsamt.forward.synthetic import LayeredModel

        models = []
        for res in results:
            rho = 10.0 ** res["log_rho"]
            thick = 10.0 ** res["log_thick"]
            try:
                m = LayeredModel(
                    resistivity=np.maximum(rho, 1e-3),
                    thickness=np.maximum(thick, 1.0),
                )
            except Exception:
                m = None
            models.append(m)
        return models

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"HybridInverter1D("
            f"n_stations={self.n_sites}, "
            f"n_layers={self._ai_inv.n_layers}, "
            f"{status})"
        )

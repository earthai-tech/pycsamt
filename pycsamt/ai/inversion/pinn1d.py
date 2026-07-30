# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Physics-informed 1-D MT/CSAMT inversion.

:class:`PINNInverter1D` fits a layered Earth model
to observed EM data without labelled training examples.
Model parameters (log-resistivity and log-thickness)
are optimised with Adam by minimising:

.. math::

    \mathcal{L}(\theta) =
    \frac{1}{F}\sum_{f=1}^{F}
    \left[
      \left(
        \log_{10}\rho_a^{\rm pred}(f)
        - \log_{10}\rho_a^{\rm obs}(f)
      \right)^2
      + \left(
          \frac{\phi^{\rm pred}(f)
                - \phi^{\rm obs}(f)}{90}
        \right)^2
    \right]
    + \lambda \sum_{k=1}^{L-1}
      \!\bigl(\rho_k - \rho_{k-1}\bigr)^2

where :math:`\rho_a^{\rm pred}` comes from the Wait
(1954) 1-D MT recursion.  The recursion is
implemented in both PyTorch and TensorFlow; the
active backend is selected automatically via
:mod:`pycsamt.backends`.

Example
-------
>>> from pycsamt.ai.inversion import PINNInverter1D
>>> inv = PINNInverter1D(
...     "edi/",
...     n_layers=10,
...     depth_max=2000.0,
...     smoothness_weight=0.01,
... )  # doctest: +SKIP
>>> inv.fit(epochs=500)  # doctest: +SKIP
PINNInverter1D(n_stations=5, fitted)
>>> models = inv.predict()  # doctest: +SKIP
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .._base import BasePINNInverter
from ._pinn_ops import fit_station
from ._sites_bridge import SiteObs1D, sites_to_obs_1d

__all__ = ["PINNInverter1D"]


# backward-compat re-exports consumed by hybrid1d — torch only
try:
    from ._pinn_ops_torch import (
        _fit_station_torch as _fit_station,
    )
    from ._pinn_ops_torch import (
        _mt1d_torch,
    )
except ImportError:  # TF-only environment
    _mt1d_torch = None  # type: ignore[assignment]
    _fit_station = None  # type: ignore[assignment]


# ── PINNInverter1D ──


class PINNInverter1D(BasePINNInverter):
    r"""
    Physics-informed 1-D EM inversion.

    Fits a layered Earth model to observed MT/CSAMT
    apparent resistivity and phase by gradient descent.
    No labelled training data is required.

    Parameters
    ----------
    sites : Any
        Path, ``EDIFile``, ``EDICollection``, ``Site``,
        ``Sites``, ``APISurvey``, or iterable.
        Station data are extracted at construction time.
    solver : {'mt1d', 'csamt1d'}
        EM physics solver (TEM not supported).
    n_layers : int, default 10
        Number of earth layers including the halfspace.
    depth_max : float, default 2000.0
        Approximate investigation depth in metres, used
        to set equal initial layer thicknesses.
    smoothness_weight : float, default 0.01
        Regularisation weight :math:`\lambda` on the
        first-difference of log-resistivity.
    lr : float, default 1e-2
        Adam learning rate.
    device : str or None
        Torch device.  Auto-detects CUDA/CPU if None.
    comp : {'xy', 'yx', 'xx', 'yy'}, default 'xy'
        Impedance tensor component to use.
    recursive : bool, default True
        Passed to ``ensure_sites``.
    on_dup : str, default 'replace'
        Passed to ``ensure_sites``.
    verbose : int, default 0
        Verbosity level for site loading.
    """

    def __init__(
        self,
        sites: Any,
        *,
        solver: str = "mt1d",
        n_layers: int = 10,
        depth_max: float = 2000.0,
        smoothness_weight: float = 0.01,
        lr: float = 1e-2,
        device: str | None = None,
        comp: str = "xy",
        recursive: bool = True,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> None:
        if solver not in ("mt1d", "csamt1d"):
            raise ValueError(
                f"solver must be 'mt1d' or 'csamt1d'; got {solver!r}."
            )
        super().__init__(
            n_layers=n_layers,
            depth_max=depth_max,
            device=device,
        )
        self.solver = solver
        self.smoothness_weight = float(smoothness_weight)
        self.lr = float(lr)
        self.comp = comp
        self.recursive = recursive
        self.on_dup = on_dup
        self.verbose = verbose

        self._obs: list[SiteObs1D] = sites_to_obs_1d(
            sites,
            comp=comp,
            recursive=recursive,
            on_dup=on_dup,
            strict=True,
            verbose=verbose,
        )
        self._results: list[dict] = []

    # ── fit ──

    def fit(
        self,
        epochs: int = 500,
        *,
        verbose: bool = True,
        log_every: int = 100,
    ) -> PINNInverter1D:
        """
        Run the physics-informed optimisation.

        Parameters
        ----------
        epochs : int, default 500
            Number of Adam iterations per station.
        verbose : bool, default True
            Print progress per station.
        log_every : int, default 100
            Print epoch detail every this many steps.

        Returns
        -------
        self
        """
        self._require_backend()
        dev = self._resolve_device()
        _every = log_every if verbose else 0

        self._results = []
        for i, obs in enumerate(self._obs):
            if verbose:
                print(f"[{i + 1}/{len(self._obs)}] Inverting {obs.name} ...")
            res = fit_station(
                obs,
                n_layers=self.n_layers,
                depth_max=self.depth_max,
                lam=self.smoothness_weight,
                lr=self.lr,
                epochs=epochs,
                device=dev,
                log_every=_every,
            )
            self._results.append(res)

        self._is_fitted = True
        return self

    # ── predict ──

    def predict(self) -> list:
        """
        Return fitted layered models for all stations.

        Returns
        -------
        models : list of LayeredModel
            One per station, same order as
            :attr:`stations`.
        """
        self._check_fitted()
        from pycsamt.forward.synthetic import LayeredModel

        models = []
        for res in self._results:
            log_rho = res["log_rho"]
            log_thick = res["log_thick"]
            rho = 10.0**log_rho
            thick = 10.0**log_thick
            try:
                m = LayeredModel(
                    resistivity=np.maximum(rho, 1e-3),
                    thickness=np.maximum(thick, 1.0),
                )
            except Exception:
                m = None
            models.append(m)
        return models

    # ── diagnostics ──

    def residuals(self):
        """
        Compute observed vs predicted data for all sites.

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

        Fwd = MT1DForward if self.solver == "mt1d" else CSAMT1DForward
        rows = []
        for obs, res in zip(self._obs, self._results):
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
                        "freq": obs.freq[k],
                        "rho_obs": obs.rho_obs[k],
                        "rho_pred": rho_pred[k],
                        "phase_obs": obs.phase_obs[k],
                        "phase_pred": ph_pred[k],
                    }
                )
        return pd.DataFrame(rows)

    def loss_curves(self):
        """
        Return the Adam loss history for all stations.

        Returns
        -------
        pandas.DataFrame
            Columns: station, epoch, loss.
        """
        self._check_fitted()
        import pandas as pd

        rows = []
        for obs, res in zip(self._obs, self._results):
            for ep, val in enumerate(res["history"], start=1):
                rows.append(
                    {
                        "station": obs.name,
                        "epoch": ep,
                        "loss": val,
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

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else ("unfitted")
        return (
            f"PINNInverter1D("
            f"n_stations={self.n_sites}, "
            f"n_layers={self.n_layers}, "
            f"{status})"
        )

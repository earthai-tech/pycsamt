# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Physics-informed AI 2-D MT inversion.

:class:`PINNInverter2D` optimises a joint 2-D
resistivity section by minimising a physics-informed
loss that combines per-station 1-D MT data misfits
with lateral and vertical smoothness penalties.

The earth model is parameterised as
:math:`\log_{10}(\rho)` tensors with shape
``(n_stations, n_layers)`` for resistivity and
``(n_stations, n_layers-1)`` for layer thicknesses.
All parameters are optimised simultaneously via Adam
so the lateral-smoothness term couples adjacent
stations.

Forward model
-------------
The differentiable Wait (1954) MT-1D recursion
(implemented in :mod:`~pycsamt.ai.inversion.pinn1d`)
is applied in batch over all stations on a shared
frequency grid.

Loss function
-------------

.. math::

    \mathcal{L} =
    \underbrace{
      \frac{1}{N_v}
      \sum_{s,f}^{\text{valid}}
      \left[
        \left(
          \log_{10}\rho_a^{\rm pred} -
          \log_{10}\rho_a^{\rm obs}
        \right)^2
        + \left(
            \frac{\phi^{\rm pred} - \phi^{\rm obs}}{90}
          \right)^2
      \right]
    }_{\text{data}}
    + \lambda_z
      \underbrace{
        \frac{1}{S(L-1)}
        \sum_{s,k}
        (m_{s,k+1} - m_{s,k})^2
      }_{\text{vertical smooth}}
    + \lambda_x
      \underbrace{
        \frac{1}{(S-1)L}
        \sum_{s,k}
        (m_{s+1,k} - m_{s,k})^2
      }_{\text{lateral smooth}}

where :math:`m_{s,k} = \log_{10}\rho_{s,k}`,
:math:`N_v` is the number of valid (non-NaN)
station-frequency pairs, :math:`S` is the number
of stations, and :math:`L` is the number of layers.

Example
-------
>>> from pycsamt.ai.inversion import PINNInverter2D
>>> inv = PINNInverter2D(  # doctest: +SKIP
...     "edi/profile1/",
...     n_layers=12,
...     depth_max=3000.0,
...     epochs=300,
... )
>>> inv.fit()  # doctest: +SKIP
PINNInverter2D(n_stations=20, fitted)
>>> section = inv.resistivity_section()  # doctest: +SKIP
>>> df = inv.residuals()  # doctest: +SKIP
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .._base import BasePINNInverter
from ._pinn_ops import fit_2d_joint
from ._sites_bridge import (
    SiteObs2D,
    _interp_to_grid,
    _make_common_grid,
    sites_to_obs_2d,
)

__all__ = ["PINNInverter2D"]


# ── backward-compat re-export (hybrid2d imports this)


class PINNInverter2D(BasePINNInverter):
    r"""
    Physics-informed 2-D MT inversion.

    Optimises a pseudo-2D resistivity section by
    minimising the data misfit with the 1-D MT forward
    plus lateral and vertical smoothness penalties.

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
        Which observed polarisation to use.
        ``'both'`` averages TE and TM data misfits.
    smoothness_weight : float, default 0.01
        Vertical smoothness weight
        :math:`\lambda_z`.
    lateral_weight : float, default 0.005
        Lateral smoothness weight
        :math:`\lambda_x`.
    epochs : int, default 300
        Adam iterations.
    lr : float, default 1e-2
        Adam learning rate.
    comp_te : str, default ``'xy'``
        Impedance tensor component for TE mode.
    comp_tm : str, default ``'yx'``
        Impedance tensor component for TM mode.
    device : str or None
        Torch compute device (auto-detects if None).
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
        lateral_weight: float = 0.005,
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
            raise ValueError(
                f"mode must be 'te', 'tm', or 'both'; got {mode!r}."
            )
        super().__init__(
            n_layers=n_layers,
            depth_max=depth_max,
            device=device,
        )
        self.n_freqs = int(n_freqs)
        self.mode = mode
        self.smoothness_weight = float(smoothness_weight)
        self.lateral_weight = float(lateral_weight)
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
        self._result: dict | None = None

    # ── fit ──

    def fit(
        self,
        *,
        verbose: bool = True,
        log_every: int = 50,
    ) -> PINNInverter2D:
        """
        Run the joint 2-D physics-informed inversion.

        Parameters
        ----------
        verbose : bool, default True
        log_every : int, default 50

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
                f"PINNInverter2D: optimising "
                f"{S} stations x {self.n_layers} "
                f"layers ({self.epochs} epochs) ..."
            )

        _every = log_every if verbose else 0
        self._result = fit_2d_joint(
            lr_obs,
            ph_obs,
            self._freqs_grid,
            n_layers=self.n_layers,
            depth_max=self.depth_max,
            lam_z=self.smoothness_weight,
            lam_x=self.lateral_weight,
            lr=self.lr,
            epochs=self.epochs,
            device=dev,
            log_every=_every,
            verbose=verbose,
        )
        self._history = self._result["history"]
        self._is_fitted = True
        return self

    # ── outputs ──

    def resistivity_section(self, *, as_log10: bool = True) -> np.ndarray:
        """
        Return the 2-D resistivity section.

        Parameters
        ----------
        as_log10 : bool, default True
            If True return log10(rho); else linear rho.

        Returns
        -------
        ndarray (n_layers, n_stations)
        """
        self._check_fitted()
        lr = self._result["log_rho"]  # (S, L)
        section = lr.T  # (L, S)
        if as_log10:
            return section
        return 10.0**section

    def thickness_section(self) -> np.ndarray:
        """
        Return layer thicknesses in metres.

        Returns
        -------
        ndarray (n_layers-1, n_stations)
        """
        self._check_fitted()
        lt = self._result["log_thick"]  # (S, L-1)
        return (10.0**lt).T  # (L-1, S)

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
                "epoch": range(1, len(self._result["history"]) + 1),
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

        lr_2d = self._result["log_rho"]  # (S, L)
        lt_2d = self._result["log_thick"]  # (S, L-1)
        rows = []
        for i, obs in enumerate(self._obs):
            rho_i = lr_2d[i]
            th_i = lt_2d[i]
            try:
                from pycsamt.forward.em1d import (
                    MT1DForward,
                )
                from pycsamt.forward.synthetic import (
                    LayeredModel,
                )

                m = LayeredModel(
                    resistivity=np.maximum(10.0**rho_i, 1e-3),
                    thickness=np.maximum(10.0**th_i, 1.0),
                )
                resp = MT1DForward(obs.freq).run(m)
                rp = resp.rho_a
                pp = resp.phase
            except Exception:
                rp = np.full_like(obs.freq, np.nan)
                pp = np.full_like(obs.freq, np.nan)

            rho_obs = obs.rho_te if self.mode in ("te", "both") else obs.rho_tm
            ph_obs = (
                obs.phase_te if self.mode in ("te", "both") else obs.phase_tm
            )
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
        """Station names in profile order."""
        return [o.name for o in self._obs]

    @property
    def n_sites(self) -> int:
        """Number of loaded stations."""
        return len(self._obs)

    # ── internals ──

    def _build_obs_arrays(
        self,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Build (log_rho_obs, ph_obs) arrays for opt.

        Returns
        -------
        log_rho_obs : ndarray (S, F)
        ph_obs : ndarray (S, F)
            NaN where interpolation is out of range.
        """
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
                # mode='both': average TE and TM
                rho_src = 0.5 * (o.rho_te + o.rho_tm)
                ph_src = 0.5 * (o.phase_te + o.phase_tm)

            lr_g, ph_g = _interp_to_grid(
                o.freq,
                rho_src,
                ph_src,
                self._freqs_grid,
            )
            lr_obs[i] = lr_g
            ph_obs[i] = ph_g

        return lr_obs, ph_obs

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"PINNInverter2D("
            f"n_stations={self.n_sites}, "
            f"n_layers={self.n_layers}, "
            f"{status})"
        )

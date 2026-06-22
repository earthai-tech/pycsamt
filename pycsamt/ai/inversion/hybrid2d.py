# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Hybrid AI + physics 2-D MT inversion.

:class:`HybridInverter2D` runs a two-stage workflow:

**Stage 1 — AI initial 2-D section**

A pre-trained
:class:`~pycsamt.ai.inversion.inv2d.EMInverter2D`
maps the EM profile panel directly to a 2-D
resistivity section::

    X_panel --> EMInverter2D --> rho_2d_0

**Stage 2 — Joint physics refinement**

Starting from ``rho_2d_0``, the same joint
physics-informed optimisation used by
:class:`~pycsamt.ai.inversion.pinn2d.PINNInverter2D`
refines all stations simultaneously:

.. math::

    \mathcal{L}(\theta) =
    \mathcal{L}_{\rm data}
    + \lambda_z \mathcal{L}_{\rm vert}
    + \lambda_x \mathcal{L}_{\rm lat}

Because the AI starting model is physically
plausible, Stage 2 converges significantly faster
than a randomly initialised PINN.

Example
-------
>>> from pycsamt.ai.inversion import (
...     EMInverter2D, HybridInverter2D,
... )
>>> ai2d = EMInverter2D.load(            # doctest: +SKIP
...     "checkpoints/unet2d.npz"
... )
>>> inv = HybridInverter2D(              # doctest: +SKIP
...     "edi/profile1/",
...     ai_inverter=ai2d,
...     epochs=150,
...     smoothness_weight=0.005,
... )
>>> inv.fit()                            # doctest: +SKIP
HybridInverter2D(n_stations=20, fitted)
>>> section = inv.resistivity_section()  # doctest: +SKIP
>>> s1 = inv.stage1_section()           # doctest: +SKIP
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .._base import BaseHybridInverter
from ._sites_bridge import (
    SiteObs2D,
    sites_to_obs_2d,
    sites_to_panel_2d,
    _make_common_grid,
    _interp_to_grid,
)
from ._pinn_ops import fit_2d_joint

__all__ = ["HybridInverter2D"]


# ── HybridInverter2D ──


class HybridInverter2D(BaseHybridInverter):
    r"""
    Two-stage hybrid AI + physics 2-D inversion.

    Parameters
    ----------
    sites : Any
        Path, ``EDIFile``, ``EDICollection``, ``Site``,
        ``Sites``, ``APISurvey``, or iterable.
    ai_inverter : EMInverter2D or str or Path
        Pre-trained
        :class:`~pycsamt.ai.inversion.inv2d.EMInverter2D`
        (fitted) or path to a ``.npz`` checkpoint.
    n_layers : int or None
        Layers per station in Stage 2.  Defaults to
        ``ai_inverter.n_depth``.
    depth_max : float, default 2000.0
        Total depth for uniform thickness init when
        Stage 1 does not provide thicknesses.
    n_freqs : int, default 32
        Frequency-grid size for the panel and the
        shared optimisation grid.
    mode : {'te', 'tm', 'both'}, default 'te'
        Data polarisation used in Stage 2.
    smoothness_weight : float, default 0.005
        Vertical smoothness weight.
    lateral_weight : float, default 0.003
        Lateral smoothness weight.
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
        ai_inverter: Union[Any, str, Path],
        *,
        n_layers: Optional[int] = None,
        depth_max: float = 2000.0,
        n_freqs: int = 32,
        mode: str = "te",
        smoothness_weight: float = 0.005,
        lateral_weight: float = 0.003,
        epochs: int = 150,
        lr: float = 5e-3,
        comp_te: str = "xy",
        comp_tm: str = "yx",
        device: Optional[str] = None,
        recursive: bool = True,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> None:
        if mode not in ("te", "tm", "both"):
            raise ValueError(
                "mode must be 'te', 'tm', or 'both'; "
                f"got {mode!r}."
            )
        super().__init__(
            depth_max=depth_max,
            device=device,
        )
        self.n_freqs = int(n_freqs)
        self.mode = mode
        self.smoothness_weight = float(
            smoothness_weight
        )
        self.lateral_weight = float(lateral_weight)
        self.epochs = int(epochs)
        self.lr = float(lr)
        self.comp_te = comp_te
        self.comp_tm = comp_tm
        self.verbose = verbose

        self._ai_inv = self._load_ai_inverter(
            ai_inverter
        )
        self.n_layers = int(
            n_layers
            if n_layers is not None
            else self._ai_inv.n_depth
        )

        self._obs: List[SiteObs2D] = (
            sites_to_obs_2d(
                sites,
                comp_te=comp_te,
                comp_tm=comp_tm,
                recursive=recursive,
                on_dup=on_dup,
                verbose=verbose,
            )
        )
        self._freqs_grid = _make_common_grid(
            [o.freq for o in self._obs],
            n_freqs=n_freqs,
        )
        self._stage1_log_rho: Optional[np.ndarray] = (
            None
        )
        self._result: Optional[Dict] = None

    # ── fit ──

    def fit(
        self,
        *,
        verbose: bool = True,
        log_every: int = 50,
    ) -> "HybridInverter2D":
        """
        Run both inversion stages.

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

        # Stage 1: AI initial 2-D section
        if verbose:
            print(
                "Stage 1: EMInverter2D forward pass ..."
            )
        init_log_rho = self._run_stage1(
            verbose=verbose
        )

        # Uniform initial thicknesses
        dz = max(
            self.depth_max / self.n_layers, 1.0
        )
        init_log_thick = np.full(
            (len(self._obs), self.n_layers - 1),
            np.log10(dz),
        )

        # Build observed arrays for Stage 2
        lr_obs, ph_obs = self._build_obs_arrays()

        if verbose:
            S = len(self._obs)
            print(
                f"Stage 2: joint physics refinement "
                f"({S} stations, {self.n_layers} "
                f"layers, {self.epochs} epochs) ..."
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
            init_log_rho=init_log_rho,
            init_log_thick=init_log_thick,
        )
        self._history = self._result["history"]
        self._is_fitted = True
        return self

    # ── outputs ──

    def resistivity_section(
        self, *, as_log10: bool = True
    ) -> np.ndarray:
        """
        Return the Stage-2 2-D resistivity section.

        Parameters
        ----------
        as_log10 : bool, default True

        Returns
        -------
        ndarray (n_layers, n_stations)
        """
        self._check_fitted()
        lr = self._result["log_rho"]   # (S, L)
        section = lr.T                 # (L, S)
        if as_log10:
            return section
        return 10.0 ** section

    def thickness_section(self) -> np.ndarray:
        """
        Return Stage-2 layer thicknesses in metres.

        Returns
        -------
        ndarray (n_layers-1, n_stations)
        """
        self._check_fitted()
        lt = self._result["log_thick"]  # (S, L-1)
        return (10.0 ** lt).T           # (L-1, S)

    def stage1_section(
        self, *, as_log10: bool = True
    ) -> np.ndarray:
        """
        Return the Stage-1 AI 2-D section.

        Parameters
        ----------
        as_log10 : bool, default True

        Returns
        -------
        ndarray (n_layers, n_stations)
        """
        if self._stage1_log_rho is None:
            raise RuntimeError(
                "Call fit() to populate Stage-1 "
                "results."
            )
        s = self._stage1_log_rho.T  # (L, S)
        if as_log10:
            return s
        return 10.0 ** s

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
            Which stage's models to evaluate.

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

        if stage == 1:
            if self._stage1_log_rho is None:
                raise RuntimeError(
                    "Call fit() first."
                )
            lr_2d = self._stage1_log_rho  # (S, L)
            dz = max(
                self.depth_max / self.n_layers, 1.0
            )
            lt_2d = np.full(
                (len(self._obs), self.n_layers - 1),
                np.log10(dz),
            )
        else:
            lr_2d = self._result["log_rho"]
            lt_2d = self._result["log_thick"]

        rows = []
        for i, obs in enumerate(self._obs):
            try:
                m = LayeredModel(
                    resistivity=np.maximum(
                        10.0 ** lr_2d[i], 1e-3
                    ),
                    thickness=np.maximum(
                        10.0 ** lt_2d[i], 1.0
                    ),
                )
                resp = MT1DForward(obs.freq).run(m)
                rp = resp.rho_a
                pp = resp.phase
            except Exception:
                rp = np.full_like(obs.freq, np.nan)
                pp = np.full_like(obs.freq, np.nan)

            rho_obs = (
                obs.rho_te
                if self.mode in ("te", "both")
                else obs.rho_tm
            )
            ph_obs = (
                obs.phase_te
                if self.mode in ("te", "both")
                else obs.phase_tm
            )
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

    # ── read-only properties ──

    @property
    def stations(self) -> List[str]:
        """Station names in profile order."""
        return [o.name for o in self._obs]

    @property
    def n_sites(self) -> int:
        """Number of loaded stations."""
        return len(self._obs)

    # ── internals ──

    def _load_ai_inverter(self, ai_inverter):
        from .inv2d import EMInverter2D

        if isinstance(ai_inverter, EMInverter2D):
            if not ai_inverter._is_fitted:
                raise ValueError(
                    "ai_inverter must be a fitted "
                    "EMInverter2D instance."
                )
            return ai_inverter
        if isinstance(ai_inverter, (str, Path)):
            return EMInverter2D.load(
                Path(ai_inverter)
            )
        raise TypeError(
            "ai_inverter must be a fitted "
            "EMInverter2D or a path to a checkpoint; "
            f"got {type(ai_inverter)!r}."
        )

    def _run_stage1(
        self, *, verbose: bool
    ) -> np.ndarray:
        """
        Apply EMInverter2D to extract an initial
        2-D section.

        Returns
        -------
        init_log_rho : ndarray (S, n_layers)
            Per-station log10(rho) initialisation.
        """
        ai = self._ai_inv
        S = len(self._obs)
        n_ch = ai.n_components
        n_f = ai.n_freqs
        n_sta_ai = ai.n_stations

        # Build the input panel with obs n_stations,
        # then resize to ai.n_stations if needed.
        panel, _, _ = sites_to_panel_2d(
            self._obs,
            n_freqs=n_f,
            n_components=n_ch,
            comp_te=self.comp_te,
            comp_tm=self.comp_tm,
        )
        # panel: (1, n_ch, n_f, S)

        if panel.shape[3] != n_sta_ai:
            # Spatially resample to match AI model
            panel = _resize_panel_stations(
                panel, n_sta_ai
            )

        # Fill NaN with column median before predict
        panel = _fill_nan_panel(panel)

        # rho_2d: (1, n_depth, n_sta_ai) in log10
        rho_2d = ai.predict(
            panel.astype(np.float32),
            as_log_rho=True,
        )
        rho_2d = rho_2d[0]  # (n_depth, n_sta_ai)

        # Resample back to actual station count
        if n_sta_ai != S:
            rho_2d = _resize_section_stations(
                rho_2d, S
            )

        # rho_2d: (n_depth, S)
        # Clip to n_layers depth levels
        n_depth = rho_2d.shape[0]
        take = min(n_depth, self.n_layers)
        init_lr = np.full(
            (S, self.n_layers), 2.0
        )
        init_lr[:, :take] = rho_2d[:take, :].T

        if take < self.n_layers:
            # Pad deeper layers with deepest value
            init_lr[:, take:] = (
                init_lr[:, take - 1: take]
            )

        self._stage1_log_rho = init_lr
        return init_lr

    def _build_obs_arrays(
        self,
    ) -> Tuple[np.ndarray, np.ndarray]:
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
                rho_src = 0.5 * (
                    o.rho_te + o.rho_tm
                )
                ph_src = 0.5 * (
                    o.phase_te + o.phase_tm
                )

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
        status = (
            "fitted" if self._is_fitted else "unfitted"
        )
        return (
            f"HybridInverter2D("
            f"n_stations={self.n_sites}, "
            f"n_layers={self.n_layers}, "
            f"{status})"
        )


# ── panel resize helpers ──


def _resize_panel_stations(
    panel: np.ndarray,
    n_out: int,
) -> np.ndarray:
    """
    Bilinearly resample the station axis of a panel.

    Parameters
    ----------
    panel : ndarray (1, C, F, S_in)
    n_out : int

    Returns
    -------
    ndarray (1, C, F, n_out)
    """
    from scipy.ndimage import zoom

    _, C, F, S_in = panel.shape
    if S_in == n_out:
        return panel
    factor = n_out / S_in
    out = np.empty((1, C, F, n_out), dtype=panel.dtype)
    for c in range(C):
        # zoom over (F, S) axis
        slc = panel[0, c]          # (F, S_in)
        out[0, c] = zoom(
            slc,
            (1.0, factor),
            order=1,
            mode="nearest",
        )
    return out


def _resize_section_stations(
    section: np.ndarray,
    n_out: int,
) -> np.ndarray:
    """
    Bilinearly resample the station axis of a section.

    Parameters
    ----------
    section : ndarray (n_depth, S_in)
    n_out : int

    Returns
    -------
    ndarray (n_depth, n_out)
    """
    from scipy.ndimage import zoom

    n_depth, S_in = section.shape
    if S_in == n_out:
        return section
    factor = n_out / S_in
    return zoom(
        section,
        (1.0, factor),
        order=1,
        mode="nearest",
    )


def _fill_nan_panel(
    panel: np.ndarray,
) -> np.ndarray:
    """
    Replace NaN in each ``(F,)`` freq column with the
    channel's finite mean; fall back to 0 if all NaN.

    Parameters
    ----------
    panel : ndarray (1, C, F, S)

    Returns
    -------
    ndarray same shape, no NaN.
    """
    out = panel.copy()
    _, C, F, S = out.shape
    for c in range(C):
        for s in range(S):
            col = out[0, c, :, s]
            mask = ~np.isfinite(col)
            if mask.any():
                finite = col[~mask]
                fill = (
                    float(finite.mean())
                    if finite.size > 0
                    else 0.0
                )
                col[mask] = fill
    return out

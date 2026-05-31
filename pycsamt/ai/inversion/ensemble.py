# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
EnsembleInverter — deep ensemble for uncertainty-aware EM inversion.

Trains :math:`M` independent copies of a base :class:`EMInverter1D`
estimator with different random seeds.  Epistemic (model) uncertainty
is estimated from the inter-model variance.  Aleatoric uncertainty is
not captured here; use MC-Dropout variants for that.

Uncertainty estimation
----------------------
For a test sample :math:`\\mathbf{x}`, each ensemble member
:math:`f_m` predicts :math:`\\hat{\\mathbf{y}}_m`.  The ensemble
mean and variance are:

.. math::

    \\bar{\\mathbf{y}} &= \\frac{1}{M} \\sum_{m=1}^{M} \\hat{\\mathbf{y}}_m

    \\sigma^2 &= \\frac{1}{M-1} \\sum_{m=1}^{M}
               (\\hat{\\mathbf{y}}_m - \\bar{\\mathbf{y}})^2

The ensemble RMSE (expected calibration) is typically lower than any
single member; :math:`\\sigma` provides a calibrated uncertainty
estimate for downstream probabilistic interpretation.

References
----------
Lakshminarayanan B. et al. (2017) *NeurIPS* — Simple and Scalable
Predictive Uncertainty Estimation using Deep Ensembles.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

__all__ = ["EnsembleInverter"]


class EnsembleInverter:
    """
    Deep ensemble of :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`
    models.

    Parameters
    ----------
    base_estimator : EMInverter1D
        Template estimator (unfitted).  Its hyper-parameters are copied
        for every ensemble member.
    n_estimators : int, default 5
        Number of ensemble members.
    seeds : list of int or None
        Per-member random seeds.  When ``None``, seeds
        ``[0, 1, …, n_estimators-1]`` are used.

    Examples
    --------
    >>> from pycsamt.ai.inversion import EMInverter1D, EnsembleInverter
    >>> base = EMInverter1D(arch="resnet", n_layers=5)
    >>> ens = EnsembleInverter(base, n_estimators=5)
    >>> ens.fit(ds, epochs=50)                           # doctest: +SKIP
    EnsembleInverter(n_estimators=5, fitted)
    >>> mean, std = ens.predict_with_uncertainty(X_test) # doctest: +SKIP
    """

    def __init__(
        self,
        base_estimator: Any,
        n_estimators: int = 5,
        seeds: Optional[Sequence[int]] = None,
    ) -> None:
        self.base_estimator = base_estimator
        self.n_estimators = int(n_estimators)
        self.seeds = (
            list(seeds) if seeds is not None
            else list(range(n_estimators))
        )
        if len(self.seeds) < self.n_estimators:
            extra = range(len(self.seeds), n_estimators)
            self.seeds += list(extra)

        self._members: List[Any] = []
        self._is_fitted: bool = False

    # ─── fit ──────────────────────────────────────────────────────────────

    def fit(
        self,
        X,
        y=None,
        *,
        epochs: int = 100,
        batch_size: int = 256,
        lr: float = 1e-3,
        patience: int = 20,
        val_frac: float = 0.1,
        grad_clip: Optional[float] = 1.0,
        verbose: bool = True,
    ) -> "EnsembleInverter":
        """
        Train all ensemble members.

        Parameters
        ----------
        X : ForwardDataset, path, or ndarray
            Training data (forwarded unchanged to each member's ``fit``).
        y : ndarray or None
        epochs, batch_size, lr, patience, val_frac, grad_clip, verbose
            Forwarded to each member's ``fit`` call; each member
            receives a unique ``seed``.

        Returns
        -------
        self
        """
        self._members = []
        for i, seed in enumerate(self.seeds[:self.n_estimators]):
            member = copy.deepcopy(self.base_estimator)
            if verbose:
                print(f"\n=== Ensemble member {i + 1}/{self.n_estimators} "
                      f"(seed={seed}) ===")
            member.fit(
                X, y,
                epochs=epochs,
                batch_size=batch_size,
                lr=lr,
                patience=patience,
                val_frac=val_frac,
                grad_clip=grad_clip,
                seed=seed,
                verbose=verbose,
            )
            self._members.append(member)

        self._is_fitted = True
        return self

    # ─── predict ──────────────────────────────────────────────────────────

    def predict(self, X: np.ndarray, **kwargs) -> np.ndarray:
        """
        Return the mean ensemble prediction.

        Parameters
        ----------
        X : ndarray (n_samples, n_features)
        **kwargs : forwarded to each member's ``predict``

        Returns
        -------
        y_mean : ndarray (n_samples, n_params)
        """
        preds = self._all_predictions(X, **kwargs)
        return preds.mean(axis=0)

    def predict_with_uncertainty(
        self,
        X: np.ndarray,
        **kwargs,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return (mean, std) ensemble prediction.

        Parameters
        ----------
        X : ndarray (n_samples, n_features)

        Returns
        -------
        mean : ndarray (n_samples, n_params)
        std : ndarray (n_samples, n_params)
            Inter-member standard deviation — a measure of epistemic
            uncertainty.
        """
        preds = self._all_predictions(X, **kwargs)  # (M, n, p)
        return preds.mean(axis=0), preds.std(axis=0, ddof=min(1, len(preds) - 1))

    def predict_quantiles(
        self,
        X: np.ndarray,
        q: Sequence[float] = (0.05, 0.25, 0.5, 0.75, 0.95),
        **kwargs,
    ) -> Dict[float, np.ndarray]:
        """
        Return per-quantile predictions.

        Parameters
        ----------
        X : ndarray (n_samples, n_features)
        q : sequence of float
            Quantiles in (0, 1).

        Returns
        -------
        quantiles : dict  {q_value: ndarray (n_samples, n_params)}
        """
        preds = self._all_predictions(X, **kwargs)  # (M, n, p)
        return {
            float(qi): np.quantile(preds, qi, axis=0)
            for qi in q
        }

    def score(
        self,
        X: np.ndarray,
        y: np.ndarray,
        *,
        metric: str = "rmse",
    ) -> float:
        """
        Score the ensemble mean prediction.

        Parameters
        ----------
        X : ndarray (n_samples, n_features)
        y : ndarray (n_samples, n_params)
        metric : {'rmse', 'mae', 'r2', 'relative_rmse'}

        Returns
        -------
        score : float
        """
        from .._base import _compute_metric
        y_pred = self.predict(X)
        return _compute_metric(y, y_pred, metric)

    def coverage(
        self,
        X: np.ndarray,
        y: np.ndarray,
        *,
        n_sigma: float = 1.96,
    ) -> float:
        """
        Fraction of true values within ±n_sigma * std of the ensemble mean.

        An uncertainty estimate is well-calibrated when
        ``coverage(1.96) ≈ 0.95``.

        Parameters
        ----------
        X : ndarray
        y : ndarray
        n_sigma : float

        Returns
        -------
        coverage : float ∈ [0, 1]
        """
        mean, std = self.predict_with_uncertainty(X)
        within = np.abs(y - mean) <= n_sigma * std
        mask = np.isfinite(y) & np.isfinite(mean) & np.isfinite(std)
        if mask.any():
            return float(within[mask].mean())
        return np.nan

    # ─── serialisation ────────────────────────────────────────────────────

    def save(self, path: Union[str, Path]) -> None:
        """
        Save all ensemble members.

        Creates a directory ``<path>/`` containing
        ``member_00.npz``, ``member_01.npz``, ….

        Parameters
        ----------
        path : str or Path
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        for i, member in enumerate(self._members):
            member.save(path / f"member_{i:02d}.npz")
        # Save ensemble metadata
        np.savez_compressed(
            path / "_ensemble_meta.npz",
            n_estimators=np.array([self.n_estimators]),
            seeds=np.array(self.seeds),
        )

    @classmethod
    def load(
        cls,
        path: Union[str, Path],
        base_class: Optional[type] = None,
    ) -> "EnsembleInverter":
        """
        Load a saved ensemble.

        Parameters
        ----------
        path : str or Path
            Directory created by :meth:`save`.
        base_class : type or None
            Estimator class to use for loading.  Defaults to
            :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`.

        Returns
        -------
        EnsembleInverter
        """
        path = Path(path)
        if base_class is None:
            from .inv1d import EMInverter1D
            base_class = EMInverter1D

        meta = np.load(path / "_ensemble_meta.npz", allow_pickle=True)
        n = int(meta["n_estimators"][0])
        seeds = meta["seeds"].tolist()

        members = []
        for i in range(n):
            members.append(base_class.load(path / f"member_{i:02d}.npz"))

        # Reconstruct without fitting
        dummy = object.__new__(base_class)
        obj = cls(dummy, n_estimators=n, seeds=seeds)
        obj._members = members
        obj._is_fitted = True
        return obj

    # ─── internal ─────────────────────────────────────────────────────────

    def _check_fitted(self) -> None:
        if not self._is_fitted or not self._members:
            raise RuntimeError("EnsembleInverter is not fitted.  Call fit() first.")

    def _all_predictions(self, X: np.ndarray, **kwargs) -> np.ndarray:
        self._check_fitted()
        preds = np.stack(
            [m.predict(X, **kwargs) for m in self._members], axis=0
        )  # (M, n_samples, n_params)
        return preds

    # ─── plotting ─────────────────────────────────────────────────────────

    def plot_uncertainty_profile(
        self,
        X: np.ndarray,
        sample_idx: int = 0,
        *,
        y_true: Optional[np.ndarray] = None,
        n_sigma: float = 1.96,
        **plot_kwargs,
    ):
        """
        Plot a 1-D resistivity profile with uncertainty bands.

        Parameters
        ----------
        X : ndarray (n_samples, n_features)
        sample_idx : int
            Which sample to plot.
        y_true : ndarray (n_params,) or None
            True model vector for overlay.
        n_sigma : float
            Band width in standard deviations.
        **plot_kwargs
            Forwarded to :func:`~pycsamt.ai.plot.diagnostics.plot_uncertainty_bands`.

        Returns
        -------
        fig : Figure
        """
        from ..plot.diagnostics import plot_uncertainty_bands

        mean, std = self.predict_with_uncertainty(X)
        m = mean[sample_idx]   # (n_params,)
        s = std[sample_idx]

        n_layers = self.base_estimator.n_layers
        rho_m = m[:n_layers]
        rho_hi = rho_m + n_sigma * s[:n_layers]
        rho_lo = rho_m - n_sigma * s[:n_layers]

        # Build depth axis from predicted thicknesses
        thick = m[n_layers:]
        depths = np.concatenate([[0.0], np.cumsum(np.maximum(thick, 1.0))])

        rho_true_plot = None
        if y_true is not None:
            rho_true_plot = np.asarray(y_true)[:n_layers]

        fig = plot_uncertainty_bands(
            depths[:n_layers], rho_m, rho_hi, rho_lo,
            y_true=rho_true_plot,
            xlabel=r"$\log_{10}(\rho)$ (Ω·m)",
            ylabel="Depth (m)",
            title=f"Sample {sample_idx} — ensemble uncertainty",
            **plot_kwargs,
        )
        return fig

    # ─── dunder ───────────────────────────────────────────────────────────

    def __len__(self) -> int:
        return len(self._members)

    def __getitem__(self, idx: int):
        return self._members[idx]

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"EnsembleInverter(n_estimators={self.n_estimators}, {status})"
        )

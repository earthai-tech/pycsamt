# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PyTorch Dataset wrapper and normalisation utilities for EM training data.

:class:`Normalizer` applies z-score normalisation (fit on training data,
applied to val/test without refit).  :class:`EMDataset` wraps a
:class:`~pycsamt.forward.batch.ForwardDataset` and handles:

* z-score normalisation of X and y
* log₁₀ transform of thickness values in y
* NaN masking (NaN targets are left as ``torch.nan`` for masked loss)
* Optional on-the-fly noise augmentation

Usage
-----
>>> from pycsamt.forward.batch import generate_dataset
>>> from pycsamt.ai.training.dataset import EMDataset, Normalizer
>>> ds = generate_dataset(n_samples=1000, seed=0, verbose=False)
>>> train_ds = EMDataset(ds, log_thickness=True)
>>> train_ds.x_norm.mean.shape   # (n_features,)
"""
from __future__ import annotations

import numpy as np

__all__ = ["Normalizer", "EMDataset"]


# ─────────────────────────────────────────────────────────────────────────────
# Normalizer
# ─────────────────────────────────────────────────────────────────────────────

class Normalizer:
    """
    Z-score normaliser that handles NaN values.

    Parameters
    ----------
    eps : float
        Small constant added to std to avoid division by zero.

    Attributes
    ----------
    mean, std : ndarray
        Per-feature statistics computed by :meth:`fit`.
    """

    def __init__(self, eps: float = 1e-8):
        self.eps = eps
        self.mean: np.ndarray | None = None
        self.std: np.ndarray | None = None

    def fit(self, X: np.ndarray) -> Normalizer:
        """Compute mean and std from *X*, ignoring NaN."""
        self.mean = np.nanmean(X, axis=0)
        self.std = np.nanstd(X, axis=0) + self.eps
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Apply z-score normalisation."""
        if self.mean is None:
            raise RuntimeError("Call fit() before transform().")
        return (X - self.mean) / self.std

    def inverse_transform(self, X_norm: np.ndarray) -> np.ndarray:
        """Reverse z-score normalisation."""
        if self.mean is None:
            raise RuntimeError("Call fit() before inverse_transform().")
        return X_norm * self.std + self.mean

    def fit_transform(self, X: np.ndarray) -> np.ndarray:
        return self.fit(X).transform(X)

    def to_dict(self) -> dict:
        return {"mean": self.mean.tolist(), "std": self.std.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> Normalizer:
        obj = cls()
        obj.mean = np.asarray(d["mean"])
        obj.std = np.asarray(d["std"])
        return obj


# ─────────────────────────────────────────────────────────────────────────────
# EMDataset
# ─────────────────────────────────────────────────────────────────────────────

class EMDataset:
    """
    PyTorch ``Dataset``-compatible wrapper for a
    :class:`~pycsamt.forward.batch.ForwardDataset`.

    Handles normalisation, log₁₀ thickness transform, and optional
    on-the-fly noise augmentation.

    Parameters
    ----------
    forward_ds : ForwardDataset
        Source data produced by :func:`~pycsamt.forward.batch.generate_dataset`.
    n_layers : int or None
        If given, only samples with ``meta.n_layers == n_layers`` are kept.
        ``None`` keeps all samples (variable n_layers, NaN-padded).
    log_thickness : bool
        Apply log₁₀ to thickness values in y before normalising.
        Strongly recommended when thicknesses span > 2 orders of magnitude.
    x_norm : Normalizer or None
        Pre-fitted normaliser for X.  If ``None`` a new one is fitted on
        this dataset (use on training split only).
    y_norm : Normalizer or None
        Pre-fitted normaliser for y.
    augment_noise : float
        If > 0, add Gaussian noise (this level) to X on-the-fly each epoch.

    Attributes
    ----------
    X : ndarray, shape (n_valid, n_features)
    y : ndarray, shape (n_valid, n_params)
    x_norm, y_norm : Normalizer
    n_features : int
    n_params : int
    """

    def __init__(
        self,
        forward_ds,
        *,
        n_layers: int | None = None,
        log_thickness: bool = True,
        x_norm: Normalizer | None = None,
        y_norm: Normalizer | None = None,
        augment_noise: float = 0.0,
    ):
        X_raw = forward_ds.X.astype(np.float32)
        y_raw = forward_ds.y.astype(np.float32)

        # Filter by n_layers if requested
        if n_layers is not None and forward_ds.meta is not None:
            mask = forward_ds.meta["n_layers"] == n_layers
            if mask.sum() == 0:
                raise ValueError(
                    f"No samples with n_layers={n_layers} found in dataset."
                )
            X_raw = X_raw[mask]
            y_raw = y_raw[mask]
            # Trim y to exact size: n_layers rho + (n_layers-1) thick
            n_out = 2 * n_layers - 1
            y_raw = y_raw[:, :n_out]

        # Log-transform thickness values
        if log_thickness:
            y_raw = self._log_thickness(y_raw, n_layers)

        # Replace inf with NaN for safety
        X_raw = np.where(np.isfinite(X_raw), X_raw, np.nan)
        y_raw = np.where(np.isfinite(y_raw), y_raw, np.nan)

        # Normalise X
        if x_norm is None:
            self.x_norm = Normalizer().fit(X_raw)
        else:
            self.x_norm = x_norm
        X_norm = self.x_norm.transform(X_raw)

        # Normalise y
        if y_norm is None:
            self.y_norm = Normalizer().fit(y_raw)
        else:
            self.y_norm = y_norm
        y_norm_arr = self.y_norm.transform(y_raw)

        # Store as float32, keeping NaN for masked loss
        self.X = X_norm.astype(np.float32)
        self.y = y_norm_arr.astype(np.float32)
        self.augment_noise = float(augment_noise)
        self.n_features = self.X.shape[1]
        self.n_params = self.y.shape[1]
        self._n_layers = n_layers
        self._log_thickness = log_thickness

    # ─── static helpers ───────────────────────────────────────────────────

    @staticmethod
    def _log_thickness(y: np.ndarray, n_layers: int | None) -> np.ndarray:
        """
        Apply log₁₀ to the thickness sub-vector in *y*.

        If ``n_layers`` is known, thicknesses are at indices
        ``n_layers : 2*n_layers - 1``.  Otherwise, assume the second
        half of the non-NaN columns contains thicknesses.
        """
        y = y.copy()
        if n_layers is not None:
            idx = slice(n_layers, 2 * n_layers - 1)
        else:
            n_cols = y.shape[1]
            mid = n_cols // 2
            idx = slice(mid, n_cols)
        t = y[:, idx]
        valid = t > 0
        y[:, idx] = np.where(valid, np.log10(np.maximum(t, 1e-6)), np.nan)
        return y

    # ─── PyTorch Dataset protocol ────────────────────────────────────────

    def __len__(self) -> int:
        return len(self.X)

    def __getitem__(self, idx):
        """Return ``(x_tensor, y_tensor)`` pair."""
        try:
            import torch
        except ImportError:
            raise ImportError("PyTorch required for EMDataset.__getitem__")
        x = torch.from_numpy(self.X[idx])
        if self.augment_noise > 0.0:
            x = x + torch.randn_like(x) * self.augment_noise
        y = torch.from_numpy(self.y[idx])
        return x, y

    # ─── helpers ─────────────────────────────────────────────────────────

    def inverse_y(self, y_norm: np.ndarray) -> np.ndarray:
        """
        Undo normalisation + log₁₀ thickness transform on predicted y.

        Returns raw parameter vector: ``[rho_0…rho_{n-1}, thick_0…thick_{n-2}]``
        where rho values are in Ω·m (not log) and thicknesses in metres.
        """
        y = self.y_norm.inverse_transform(y_norm).astype(float)
        if self._log_thickness and self._n_layers is not None:
            n = self._n_layers
            y[:, n:] = 10.0 ** y[:, n:]
        return y

    def inverse_x(self, x_norm: np.ndarray) -> np.ndarray:
        """Undo X normalisation."""
        return self.x_norm.inverse_transform(x_norm)

    def split(
        self,
        val_frac: float = 0.1,
        seed: int | None = None,
    ) -> tuple[EMDataset, EMDataset]:
        """
        Split into train and validation ``EMDataset`` objects.

        The val split uses the same normaliser fitted on the train split.

        Returns
        -------
        (train_ds, val_ds) : tuple of EMDataset
        """
        rng = np.random.default_rng(seed)
        n = len(self.X)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx = idx[:n_val]
        train_idx = idx[n_val:]

        def _subset(indices):
            ds = EMDataset.__new__(EMDataset)
            ds.X = self.X[indices]
            ds.y = self.y[indices]
            ds.x_norm = self.x_norm
            ds.y_norm = self.y_norm
            ds.augment_noise = self.augment_noise
            ds.n_features = self.n_features
            ds.n_params = self.n_params
            ds._n_layers = self._n_layers
            ds._log_thickness = self._log_thickness
            return ds

        train_ds = _subset(train_idx)
        train_ds.augment_noise = self.augment_noise   # keep augmentation on train
        val_ds = _subset(val_idx)
        val_ds.augment_noise = 0.0                    # no augmentation on val
        return train_ds, val_ds

    def __repr__(self) -> str:
        return (
            f"EMDataset(n={len(self.X)}, "
            f"n_features={self.n_features}, n_params={self.n_params})"
        )

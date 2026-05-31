# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AnomalyDetector — unsupervised profile-level anomaly detection.

A fully-connected autoencoder is trained on per-site feature vectors
derived from a survey profile.  Sites whose reconstruction error
significantly exceeds the training distribution are flagged as
anomalous — they may represent bad data, equipment problems, or
genuinely anomalous geology that warrants closer inspection.

Architecture
------------

.. math::

    \\hat{\\mathbf{x}} = g(f(\\mathbf{x}))

where :math:`f: \\mathbb{R}^d \\to \\mathbb{R}^k` is the encoder
(:math:`k \\ll d`) and :math:`g: \\mathbb{R}^k \\to \\mathbb{R}^d` is
the decoder.  The anomaly score for site :math:`i` is

.. math::

    s_i = \\|\\mathbf{x}_i - \\hat{\\mathbf{x}}_i\\|_2^2 / d.

A site is flagged when :math:`s_i` exceeds the
``threshold_percentile``-th percentile of training scores.

When PyTorch is not available the detector falls back to PCA-based
reconstruction using :class:`sklearn.decomposition.PCA`.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .._base import BaseEMProcessor

__all__ = ["AnomalyDetector"]


# ─────────────────────────────────────────────────────────────────────────────
# Internal network builder
# ─────────────────────────────────────────────────────────────────────────────

def _build_fc_ae(
    n_features: int,
    latent_dim: int,
    channels: Tuple[int, ...],
) -> Any:
    """Build a fully-connected autoencoder (FC-AE) module."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError("PyTorch is required to build the FC autoencoder") from exc

    def _mlp_block(in_d: int, out_d: int, activation: bool = True):
        layers: list = [nn.Linear(in_d, out_d), nn.BatchNorm1d(out_d)]
        if activation:
            layers.append(nn.ReLU())
        return nn.Sequential(*layers)

    ch = list(channels)

    class _FCAE(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            enc_dims = [n_features] + ch + [latent_dim]
            dec_dims = [latent_dim] + ch[::-1] + [n_features]

            enc_layers: list = []
            for i in range(len(enc_dims) - 1):
                enc_layers.append(
                    _mlp_block(enc_dims[i], enc_dims[i + 1],
                               activation=(i < len(enc_dims) - 2))
                )
            self.encoder = nn.Sequential(*enc_layers)

            dec_layers: list = []
            for i in range(len(dec_dims) - 1):
                dec_layers.append(
                    _mlp_block(dec_dims[i], dec_dims[i + 1],
                               activation=(i < len(dec_dims) - 2))
                )
            self.decoder = nn.Sequential(*dec_layers)

        def forward(self, x):
            return self.decoder(self.encoder(x))

    return _FCAE()


# ─────────────────────────────────────────────────────────────────────────────
# AnomalyDetector
# ─────────────────────────────────────────────────────────────────────────────

class AnomalyDetector(BaseEMProcessor):
    """
    Profile-level unsupervised anomaly detector.

    Parameters
    ----------
    n_features : int
        Feature vector length per site
        (typically ``n_freqs × n_components``).
    latent_dim : int, default 32
        Bottleneck dimension of the autoencoder.
    channels : tuple of int, default (128, 64)
        Hidden layer widths in the encoder (decoder is mirrored).
    threshold_percentile : float, default 95.0
        Sites with reconstruction error above this percentile of the
        training distribution are flagged as anomalous.
    device : str or None
        Compute device.  Ignored when using the PCA fallback.

    Notes
    -----
    When PyTorch is not installed the detector silently uses a PCA
    reconstruction model via scikit-learn.  The interface is identical;
    only the accuracy differs.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randn(100, 120).astype("float32")
    >>> det = AnomalyDetector(n_features=120, latent_dim=16)
    >>> det.fit(X, epochs=10, verbose=False)          # doctest: +SKIP
    AnomalyDetector(n_features=120, latent_dim=16)
    >>> scores = det.transform(X)                     # doctest: +SKIP
    >>> flags = det.flag_anomalies(X)                 # doctest: +SKIP
    """

    def __init__(
        self,
        n_features: int,
        latent_dim: int = 32,
        *,
        channels: Tuple[int, ...] = (128, 64),
        threshold_percentile: float = 95.0,
        device: Optional[str] = None,
    ) -> None:
        self.n_features = int(n_features)
        self.latent_dim = int(latent_dim)
        self.channels = tuple(channels)
        self.threshold_percentile = float(threshold_percentile)
        self.device = device

        self._network: Any = None
        self._pca: Any = None  # sklearn PCA fallback
        self._use_pca: bool = False
        self._threshold: Optional[float] = None
        self._x_mean: Optional[np.ndarray] = None
        self._x_std: Optional[np.ndarray] = None
        self._is_fitted: bool = False
        self._history: Dict[str, list] = {}

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(
        self,
        X: np.ndarray,
        *,
        epochs: int = 80,
        batch_size: int = 32,
        lr: float = 1e-3,
        val_frac: float = 0.1,
        seed: Optional[int] = None,
        verbose: bool = True,
    ) -> "AnomalyDetector":
        """
        Train the anomaly detector on profile data.

        Parameters
        ----------
        X : ndarray, shape (n_sites, n_features)
            Per-site feature vectors (normal / clean data).
        epochs : int
            Training epochs (ignored when using PCA fallback).
        batch_size : int
        lr : float
        val_frac : float
        seed : int or None
        verbose : bool

        Returns
        -------
        self
        """
        X = np.asarray(X, dtype=np.float32)
        # replace NaN columns with 0
        X = np.where(np.isfinite(X), X, 0.0)

        self._x_mean = X.mean(axis=0, keepdims=True)
        self._x_std = X.std(axis=0, keepdims=True) + 1e-8
        Xn = (X - self._x_mean) / self._x_std

        # Try torch; fall back to PCA
        try:
            self._fit_torch(Xn, epochs=epochs, batch_size=batch_size,
                            lr=lr, val_frac=val_frac, seed=seed, verbose=verbose)
            self._use_pca = False
        except ImportError:
            self._fit_pca(Xn, verbose=verbose)
            self._use_pca = True

        # Calibrate threshold on training reconstruction errors
        train_scores = self._reconstruction_error(Xn)
        self._threshold = float(np.nanpercentile(train_scores, self.threshold_percentile))
        self._is_fitted = True
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Compute per-site reconstruction error (anomaly score).

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)

        Returns
        -------
        scores : ndarray, shape (n_samples,)
            Per-site mean squared reconstruction error.  Higher values
            indicate more anomalous sites.
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")
        X = np.asarray(X, dtype=np.float32)
        X = np.where(np.isfinite(X), X, 0.0)
        Xn = (X - self._x_mean) / self._x_std
        return self._reconstruction_error(Xn)

    def flag_anomalies(self, X: np.ndarray) -> np.ndarray:
        """
        Return a boolean mask of anomalous sites.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)

        Returns
        -------
        flags : bool ndarray, shape (n_samples,)
            ``True`` for sites flagged as anomalous.
        """
        scores = self.transform(X)
        return scores > self._threshold

    # ─── internal ─────────────────────────────────────────────────────────

    def _fit_torch(
        self, Xn: np.ndarray, *, epochs: int, batch_size: int,
        lr: float, val_frac: float, seed: Optional[int], verbose: bool,
    ) -> None:
        import torch
        import torch.nn as nn
        from torch.utils.data import TensorDataset, DataLoader

        rng = np.random.default_rng(seed)
        dev = self._resolve_device()

        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx, trn_idx = idx[:n_val], idx[n_val:]
        Xtr, Xva = Xn[trn_idx], Xn[val_idx]

        self._network = _build_fc_ae(
            self.n_features, self.latent_dim, self.channels
        ).to(dev)

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=8, min_lr=1e-6
        )
        mse = nn.MSELoss()

        tr_ds = TensorDataset(torch.from_numpy(Xtr))
        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []

        for ep in range(1, epochs + 1):
            self._network.train()
            loader = DataLoader(tr_ds, batch_size=batch_size, shuffle=True)
            ep_loss = 0.0
            for (xb,) in loader:
                xb = xb.to(dev)
                pred = self._network(xb)
                loss = mse(pred, xb)
                opt.zero_grad()
                loss.backward()
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= len(Xtr)

            self._network.eval()
            with torch.no_grad():
                vb = torch.from_numpy(Xva).to(dev)
                v_loss = mse(self._network(vb), vb).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val:
                best_val = v_loss
                import copy
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 5) == 0 or ep == 1):
                print(f"  AnomalyDetector  ep {ep:>4d}/{epochs}  "
                      f"train={ep_loss:.5f}  val={v_loss:.5f}")

        if best_state is not None:
            self._network.load_state_dict(best_state)
        self._history = {"train_loss": train_losses, "val_loss": val_losses}

    def _fit_pca(self, Xn: np.ndarray, *, verbose: bool) -> None:
        try:
            from sklearn.decomposition import PCA
        except ImportError as exc:
            raise ImportError(
                "Either PyTorch or scikit-learn is required for AnomalyDetector"
            ) from exc

        n_comp = min(self.latent_dim, Xn.shape[0], Xn.shape[1])
        self._pca = PCA(n_components=n_comp, random_state=0)
        self._pca.fit(Xn)
        if verbose:
            ev = self._pca.explained_variance_ratio_.sum()
            print(f"  AnomalyDetector (PCA fallback)  "
                  f"k={n_comp}  explained_var={ev:.3f}")

    def _reconstruction_error(self, Xn: np.ndarray) -> np.ndarray:
        if self._use_pca or self._network is None:
            if self._pca is None:
                raise RuntimeError("Model is not fitted.")
            Xr = self._pca.inverse_transform(self._pca.transform(Xn))
        else:
            import torch
            dev = next(self._network.parameters()).device
            self._network.eval()
            with torch.no_grad():
                t = torch.from_numpy(Xn.astype(np.float32)).to(dev)
                Xr = self._network(t).cpu().numpy()
        return np.mean((Xn - Xr) ** 2, axis=1)

    def _resolve_device(self) -> str:
        if self.device is not None:
            return self.device
        try:
            import torch
            if torch.cuda.is_available():
                return "cuda"
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
        except ImportError:
            pass
        return "cpu"

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> Dict[str, Any]:
        return {
            "n_features": self.n_features,
            "latent_dim": self.latent_dim,
            "channels": list(self.channels),
            "threshold_percentile": self.threshold_percentile,
            "device": self.device,
        }

    def _get_weights(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        if self._x_mean is not None:
            out["_x_mean"] = self._x_mean
            out["_x_std"] = self._x_std
        if self._threshold is not None:
            out["_threshold"] = np.array([self._threshold])
        if self._network is not None:
            for k, v in self._network.state_dict().items():
                out[k] = v.cpu().numpy()
        elif self._pca is not None:
            try:
                import pickle, io
                buf = io.BytesIO()
                pickle.dump(self._pca, buf)
                out["_pca_pickle"] = np.frombuffer(buf.getvalue(), dtype=np.uint8)
            except Exception:
                pass
        return out

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        self._x_mean = weights.pop("_x_mean", None)
        self._x_std = weights.pop("_x_std", None)
        thr = weights.pop("_threshold", None)
        self._threshold = float(thr[0]) if thr is not None else None
        pca_blob = weights.pop("_pca_pickle", None)
        if pca_blob is not None:
            try:
                import pickle, io
                self._pca = pickle.load(io.BytesIO(bytes(pca_blob)))
                self._use_pca = True
                self._is_fitted = True
                return
            except Exception:
                pass
        if weights:
            import torch
            self._network = _build_fc_ae(
                self.n_features, self.latent_dim, self.channels
            )
            state = {k: torch.from_numpy(v) for k, v in weights.items()}
            self._network.load_state_dict(state)
            self._use_pca = False
        self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        backend = "pca" if self._use_pca else "torch"
        return (
            f"AnomalyDetector(n_features={self.n_features}, "
            f"latent_dim={self.latent_dim}, {backend}, {status})"
        )

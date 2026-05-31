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

When neither PyTorch nor TensorFlow is available the detector falls back
to PCA-based reconstruction using :class:`sklearn.decomposition.PCA`.
"""
from __future__ import annotations

import copy
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .._base import BaseEMProcessor
from .._backend_utils import resolve_device, get_weights, set_weights, active_backend

__all__ = ["AnomalyDetector"]


# ─────────────────────────────────────────────────────────────────────────────
# Network builders
# ─────────────────────────────────────────────────────────────────────────────

def _build_fc_ae_torch(
    n_features: int,
    latent_dim: int,
    channels: Tuple[int, ...],
) -> Any:
    """Fully-connected autoencoder — PyTorch."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError("PyTorch is required to build the FC autoencoder") from exc

    def _block(in_d: int, out_d: int, activation: bool = True):
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
                    _block(enc_dims[i], enc_dims[i + 1],
                           activation=(i < len(enc_dims) - 2))
                )
            self.encoder = nn.Sequential(*enc_layers)

            dec_layers: list = []
            for i in range(len(dec_dims) - 1):
                dec_layers.append(
                    _block(dec_dims[i], dec_dims[i + 1],
                           activation=(i < len(dec_dims) - 2))
                )
            self.decoder = nn.Sequential(*dec_layers)

        def forward(self, x):
            return self.decoder(self.encoder(x))

    return _FCAE()


def _build_fc_ae_tf(
    n_features: int,
    latent_dim: int,
    channels: Tuple[int, ...],
) -> Any:
    """Fully-connected autoencoder — TensorFlow/Keras."""
    try:
        import tensorflow as tf
        from tensorflow.keras import layers, Model
    except ImportError as exc:
        raise ImportError("TensorFlow is required for AnomalyDetector (TF backend)") from exc

    ch = list(channels)
    inp = tf.keras.Input(shape=(n_features,), name="input")

    # Encoder
    x = inp
    for d in ch:
        x = layers.Dense(d)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
    z = layers.Dense(latent_dim, name="latent")(x)

    # Decoder (symmetric)
    x = z
    for d in reversed(ch):
        x = layers.Dense(d)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
    out = layers.Dense(n_features, name="output")(x)

    return Model(inp, out, name="fc_ae")


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
    When neither PyTorch nor TensorFlow is installed the detector
    silently uses a PCA reconstruction model via scikit-learn.
    The interface is identical; only the accuracy differs.

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
        self._pca: Any = None
        self._use_pca: bool = False
        self._backend_name: Optional[str] = None
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
        X = np.where(np.isfinite(X), X, 0.0)

        self._x_mean = X.mean(axis=0, keepdims=True)
        self._x_std = X.std(axis=0, keepdims=True) + 1e-8
        Xn = (X - self._x_mean) / self._x_std

        try:
            self._backend_name = active_backend()
            if self._backend_name == "tensorflow":
                self._fit_tensorflow(
                    Xn, epochs=epochs, batch_size=batch_size,
                    lr=lr, val_frac=val_frac, seed=seed, verbose=verbose,
                )
            else:
                self._fit_torch(
                    Xn, epochs=epochs, batch_size=batch_size,
                    lr=lr, val_frac=val_frac, seed=seed, verbose=verbose,
                )
            self._use_pca = False
        except (RuntimeError, ImportError):
            self._backend_name = "pca"
            self._fit_pca(Xn, verbose=verbose)
            self._use_pca = True

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
            Per-site mean squared reconstruction error.
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
        """
        return self.transform(X) > self._threshold

    # ─── internal training paths ──────────────────────────────────────────

    def _fit_torch(
        self, Xn: np.ndarray, *, epochs: int, batch_size: int,
        lr: float, val_frac: float, seed: Optional[int], verbose: bool,
    ) -> None:
        import torch
        import torch.nn as nn
        from torch.utils.data import TensorDataset, DataLoader

        rng = np.random.default_rng(seed)
        dev = resolve_device(self.device)

        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx, trn_idx = idx[:n_val], idx[n_val:]
        Xtr, Xva = Xn[trn_idx], Xn[val_idx]

        self._network = _build_fc_ae_torch(
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
            ep_loss = 0.0
            for (xb,) in DataLoader(tr_ds, batch_size=batch_size, shuffle=True):
                xb = xb.to(dev)
                loss = mse(self._network(xb), xb)
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
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 5) == 0 or ep == 1):
                print(f"  AnomalyDetector  ep {ep:>4d}/{epochs}  "
                      f"train={ep_loss:.5f}  val={v_loss:.5f}")

        if best_state is not None:
            self._network.load_state_dict(best_state)
        self._history = {"train_loss": train_losses, "val_loss": val_losses}

    def _fit_tensorflow(
        self, Xn: np.ndarray, *, epochs: int, batch_size: int,
        lr: float, val_frac: float, seed: Optional[int], verbose: bool,
    ) -> None:
        import tensorflow as tf

        rng = np.random.default_rng(seed)
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx, trn_idx = idx[:n_val], idx[n_val:]
        Xtr, Xva = Xn[trn_idx], Xn[val_idx]

        dev = resolve_device(self.device)
        with tf.device(dev):
            self._network = _build_fc_ae_tf(
                self.n_features, self.latent_dim, self.channels
            )
            self._network.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
                loss="mse",
            )
            hist = self._network.fit(
                Xtr, Xtr,
                validation_data=(Xva, Xva),
                epochs=epochs,
                batch_size=batch_size,
                callbacks=[
                    tf.keras.callbacks.EarlyStopping(
                        monitor="val_loss", patience=12,
                        restore_best_weights=True, min_delta=1e-6,
                    ),
                    tf.keras.callbacks.ReduceLROnPlateau(
                        monitor="val_loss", factor=0.5, patience=8, min_lr=1e-6,
                    ),
                ],
                verbose=1 if verbose else 0,
            )
        self._history = {
            "train_loss": hist.history["loss"],
            "val_loss": hist.history["val_loss"],
        }

    def _fit_pca(self, Xn: np.ndarray, *, verbose: bool) -> None:
        try:
            from sklearn.decomposition import PCA
        except ImportError as exc:
            raise ImportError(
                "PyTorch, TensorFlow, or scikit-learn is required for AnomalyDetector"
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
        elif self._backend_name == "tensorflow":
            Xr = self._network.predict(Xn.astype(np.float32), verbose=0)
        else:
            import torch
            dev = next(self._network.parameters()).device
            self._network.eval()
            with torch.no_grad():
                t = torch.from_numpy(Xn.astype(np.float32)).to(dev)
                Xr = self._network(t).cpu().numpy()
        return np.mean((Xn - Xr) ** 2, axis=1)

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
        if self._backend_name is not None:
            out["_backend"] = np.array(self._backend_name)
        if self._network is not None:
            for k, v in get_weights(self._network).items():
                out[k] = v
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
        backend_blob = weights.pop("_backend", None)
        self._backend_name = str(backend_blob) if backend_blob is not None else "torch"

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
            if self._backend_name == "tensorflow":
                self._network = _build_fc_ae_tf(
                    self.n_features, self.latent_dim, self.channels
                )
            else:
                self._network = _build_fc_ae_torch(
                    self.n_features, self.latent_dim, self.channels
                )
            set_weights(self._network, weights)
            self._use_pca = False
        self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        backend = self._backend_name or ("pca" if self._use_pca else "torch")
        return (
            f"AnomalyDetector(n_features={self.n_features}, "
            f"latent_dim={self.latent_dim}, {backend}, {status})"
        )

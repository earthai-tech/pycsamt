# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
EMDenoiser — convolutional autoencoder for MT impedance tensor denoising.

Trains a 1-D convolutional autoencoder (CAE) on clean (or simulated)
MT data.  At inference time the trained network suppresses broadband
Gaussian noise and narrowband artefacts (powerline harmonics, dead-band
contamination) while preserving the smooth spectral shape of the true
signal.

Architecture
------------
The network is a symmetric encoder-decoder:

.. math::

    \\mathbf{z} = f_{\\text{enc}}(\\mathbf{x} + \\boldsymbol{\\varepsilon}),
    \\quad
    \\hat{\\mathbf{x}} = f_{\\text{dec}}(\\mathbf{z})

where :math:`\\boldsymbol{\\varepsilon}` is synthetic training noise and
the network minimises :math:`\\|\\mathbf{x} - \\hat{\\mathbf{x}}\\|_2^2`.

Input tensor shape
------------------
``(n_samples, n_components, n_freqs)``

* ``n_components = 4`` (default) → ``[log|Zxy|, \\phi_{xy}, log|Zyx|, \\phi_{yx}]``
* ``n_components = 8`` → all four impedance tensor components
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np

from .._backend_utils import (
    active_backend,
    get_weights,
    resolve_device,
    set_weights,
)
from .._base import BaseEMProcessor

__all__ = ["EMDenoiser", "prepare_z_features"]


# ─────────────────────────────────────────────────────────────────────────────
# Input preparation helper (no torch / TF dependency)
# ─────────────────────────────────────────────────────────────────────────────


def prepare_z_features(
    sites: Any,
    *,
    n_components: int = 4,
    log_amp: bool = True,
    freq_ref: np.ndarray | None = None,
) -> np.ndarray:
    """
    Extract an impedance feature array from a site collection.

    Parameters
    ----------
    sites : SiteCollection, dict, or list of Z objects
        MT impedance data in any format accepted by
        :func:`~pycsamt.emtools._core.ensure_sites`.
    n_components : {4, 8}
        4 uses off-diagonal components only; 8 includes diagonal.
    log_amp : bool
        If ``True``, amplitudes are transformed as
        :math:`\\log_{10}(|Z| + \\varepsilon)` before packing.
    freq_ref : ndarray or None
        Common frequency grid to interpolate onto.  When ``None``,
        the grid of the first site is used.

    Returns
    -------
    X : ndarray, shape ``(n_sites, n_components, n_freqs)``
        Feature array ready for :class:`EMDenoiser`.
    """
    try:
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
            ensure_sites,
        )
    except ImportError as exc:
        raise ImportError(
            "emtools is required for prepare_z_features"
        ) from exc

    S = ensure_sites(sites, recursive=True, on_dup="replace")
    rows: list[np.ndarray] = []
    freq_grid: np.ndarray | None = freq_ref

    for _i, ed in enumerate(_iter_items(S)):
        result = _get_z_block(ed, with_errors=False)
        if len(result) == 3:
            _, z, fr = result
        else:
            _, z, fr, _ = result
        if z is None:
            continue

        if freq_grid is None:
            freq_grid = fr

        comps: list[np.ndarray] = []

        def _amp(z_comp: np.ndarray) -> np.ndarray:
            amp = np.abs(z_comp)
            if log_amp:
                return np.log10(np.maximum(amp, 1e-24))
            return amp

        def _phase(z_comp: np.ndarray) -> np.ndarray:
            return np.degrees(np.angle(z_comp))

        comps += [_amp(z[:, 0, 1]), _phase(z[:, 0, 1])]  # Zxy
        comps += [_amp(z[:, 1, 0]), _phase(z[:, 1, 0])]  # Zyx

        if n_components == 8:
            comps += [_amp(z[:, 0, 0]), _phase(z[:, 0, 0])]  # Zxx
            comps += [_amp(z[:, 1, 1]), _phase(z[:, 1, 1])]  # Zyy

        mat = np.stack(comps, axis=0)  # (n_comp, n_freqs_site)

        if freq_grid is not fr:
            mat_interp = np.full((mat.shape[0], len(freq_grid)), np.nan)
            for ci in range(mat.shape[0]):
                mat_interp[ci] = np.interp(
                    np.log10(freq_grid + 1e-30),
                    np.log10(fr + 1e-30),
                    mat[ci],
                )
            mat = mat_interp

        rows.append(mat)

    if not rows:
        raise ValueError("No valid Z data found in the provided sites.")

    return np.stack(rows, axis=0).astype(
        np.float32
    )  # (n_sites, n_comp, n_freqs)


# ─────────────────────────────────────────────────────────────────────────────
# Network builders — PyTorch (channels-first) and TensorFlow (channels-last)
# ─────────────────────────────────────────────────────────────────────────────


def _build_cae_torch(
    n_components: int,
    n_freqs: int,
    channels: tuple[int, ...],
    dropout: float,
) -> Any:
    """1-D CAE using PyTorch Conv1d (channels-first)."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError("PyTorch is required for EMDenoiser") from exc

    mid = max(n_freqs // 4, 4)
    ch = list(channels)

    class _CAE(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Conv1d(n_components, ch[0], 5, padding=2),
                nn.BatchNorm1d(ch[0]),
                nn.LeakyReLU(0.1),
                nn.Conv1d(ch[0], ch[1], 3, padding=1),
                nn.BatchNorm1d(ch[1]),
                nn.LeakyReLU(0.1),
                nn.AdaptiveAvgPool1d(mid),
                nn.Conv1d(ch[1], ch[2], 3, padding=1),
                nn.BatchNorm1d(ch[2]),
                nn.LeakyReLU(0.1),
                nn.Dropout(dropout),
            )
            self.decoder = nn.Sequential(
                nn.Conv1d(ch[2], ch[1], 3, padding=1),
                nn.BatchNorm1d(ch[1]),
                nn.LeakyReLU(0.1),
                nn.Upsample(size=n_freqs),
                nn.Conv1d(ch[1], ch[0], 5, padding=2),
                nn.BatchNorm1d(ch[0]),
                nn.LeakyReLU(0.1),
                nn.Conv1d(ch[0], n_components, 5, padding=2),
            )

        def forward(self, x):  # x: (B, n_comp, n_freqs)
            return self.decoder(self.encoder(x))

    return _CAE()


def _build_cae_tf(
    n_components: int,
    n_freqs: int,
    channels: tuple[int, ...],
    dropout: float,
) -> Any:
    """
    1-D CAE using Keras Conv1D (channels-last).

    Input / output convention: ``(batch, n_freqs, n_components)``.
    Data is transposed by the caller from/to the canonical
    ``(batch, n_components, n_freqs)`` form.
    """
    try:
        import tensorflow as tf
        from tensorflow.keras import Model, layers
    except ImportError as exc:
        raise ImportError(
            "TensorFlow is required for EMDenoiser (TF backend)"
        ) from exc

    ch = list(channels)
    # pool_size ≈ n_freqs//(n_freqs//4) = 4; safe default for common grid sizes
    pool = max(1, n_freqs // max(n_freqs // 4, 4))

    inp = tf.keras.Input(shape=(n_freqs, n_components), name="input")

    # Encoder
    x = layers.Conv1D(ch[0], 5, padding="same")(inp)
    x = layers.BatchNormalization()(x)
    x = layers.LeakyReLU(0.1)(x)
    x = layers.Conv1D(ch[1], 3, padding="same")(x)
    x = layers.BatchNormalization()(x)
    x = layers.LeakyReLU(0.1)(x)
    x = layers.AveragePooling1D(pool_size=pool, padding="same")(x)
    x = layers.Conv1D(ch[2], 3, padding="same")(x)
    x = layers.BatchNormalization()(x)
    x = layers.LeakyReLU(0.1)(x)
    x = layers.SpatialDropout1D(dropout)(x)

    # Decoder
    x = layers.Conv1D(ch[1], 3, padding="same")(x)
    x = layers.BatchNormalization()(x)
    x = layers.LeakyReLU(0.1)(x)
    x = layers.UpSampling1D(pool)(x)
    x = layers.Conv1D(ch[0], 5, padding="same")(x)
    x = layers.BatchNormalization()(x)
    x = layers.LeakyReLU(0.1)(x)
    out = layers.Conv1D(n_components, 5, padding="same", name="output")(x)

    # Trim or pad to exact n_freqs in case pool arithmetic leaves ±1 off
    out = layers.Lambda(lambda t, nf=n_freqs: t[:, :nf, :], name="trim")(out)

    return Model(inp, out, name="cae")


# ─────────────────────────────────────────────────────────────────────────────
# EMDenoiser
# ─────────────────────────────────────────────────────────────────────────────


class EMDenoiser(BaseEMProcessor):
    """
    Convolutional autoencoder for MT impedance tensor denoising.

    The network is trained on clean (or synthetic) impedance data by
    adding controlled noise and learning to reconstruct the original
    signal.  Supports both PyTorch and TensorFlow backends via
    :mod:`pycsamt.backends`.

    Parameters
    ----------
    n_freqs : int or None, default None
        Number of frequency channels in the input.  When ``None``,
        inferred from ``X.shape[2]`` on the first call to :meth:`fit`.
    n_components : {4, 8}, default 4
        Feature channels per frequency: 4 = off-diagonal :math:`Z_{xy}`,
        :math:`Z_{yx}` amplitudes and phases; 8 = all four tensor
        components.
    channels : tuple of int, default (64, 128, 64)
        Channel widths for the encoder/decoder stages.
    dropout : float, default 0.1
        Dropout probability in the encoder bottleneck.
    device : str or None
        ``"cpu"``, ``"cuda"``, ``"mps"``, ``"/GPU:0"``, or ``None``
        for auto-detect via the active backend.

    Notes
    -----
    When neither PyTorch nor TensorFlow is available the denoiser falls
    back to a per-channel Gaussian smoothing filter (requires scipy).

    Call :func:`prepare_z_features` to convert a site collection to the
    ``(n_samples, n_components, n_freqs)`` array expected here.

    Examples
    --------
    >>> import numpy as np
    >>> X_clean = np.random.randn(200, 4, 32).astype("float32")
    >>> den = EMDenoiser()
    >>> den.fit(X_clean, epochs=5, verbose=False)   # doctest: +SKIP
    EMDenoiser(n_freqs=32, n_components=4)
    >>> X_denoised = den.transform(X_clean)          # doctest: +SKIP
    """

    def __init__(
        self,
        n_freqs: int | None = None,
        n_components: int = 4,
        *,
        channels: tuple[int, ...] = (64, 128, 64),
        dropout: float = 0.1,
        device: str | None = None,
    ) -> None:
        self.n_freqs = None if n_freqs is None else int(n_freqs)
        self.n_components = int(n_components)
        self.channels = tuple(channels)
        self.dropout = float(dropout)
        self.device = device

        self._network: Any = None
        self._use_numpy: bool = False
        self._is_fitted: bool = False
        self._backend_name: str | None = None
        self._x_mean: np.ndarray | None = None
        self._x_std: np.ndarray | None = None
        self._history: dict[str, list] = {}

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(
        self,
        X: np.ndarray,
        *,
        noise_level: float = 0.05,
        epochs: int = 50,
        batch_size: int = 64,
        lr: float = 1e-3,
        val_frac: float = 0.1,
        seed: int | None = None,
        verbose: bool = True,
    ) -> EMDenoiser:
        """
        Train the denoiser on clean impedance data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_components, n_freqs)
            Clean training data.
        noise_level : float
            Standard deviation of additive Gaussian noise relative to
            the per-channel standard deviation of ``X``.
        epochs : int
            Maximum number of training epochs.
        batch_size : int
        lr : float
            Initial Adam learning rate.
        val_frac : float
            Fraction of training samples held out for validation.
        seed : int or None
        verbose : bool

        Returns
        -------
        self
        """
        rng = np.random.default_rng(seed)
        X = np.asarray(X, dtype=np.float32)

        if X.ndim != 3:
            raise ValueError(
                f"X must be 3-D (n_samples, n_components, n_freqs), got shape {X.shape}"
            )

        if self.n_freqs is None:
            self.n_freqs = X.shape[2]
        elif X.shape[2] != self.n_freqs:
            raise ValueError(
                f"Expected n_freqs={self.n_freqs}, got {X.shape[2]}"
            )

        self._x_mean = X.mean(axis=(0, 2), keepdims=True)
        self._x_std = X.std(axis=(0, 2), keepdims=True) + 1e-8
        Xn = (X - self._x_mean) / self._x_std

        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx, trn_idx = idx[:n_val], idx[n_val:]
        Xtr, Xva = Xn[trn_idx], Xn[val_idx]

        try:
            self._backend_name = active_backend()
            if self._backend_name == "tensorflow":
                self._fit_tensorflow(
                    Xtr,
                    Xva,
                    noise_level,
                    epochs,
                    batch_size,
                    lr,
                    rng,
                    verbose,
                )
            elif self._backend_name != "none":
                self._fit_torch(
                    Xtr,
                    Xva,
                    noise_level,
                    epochs,
                    batch_size,
                    lr,
                    rng,
                    verbose,
                )
            else:
                raise RuntimeError("no backend")
            self._use_numpy = False
        except (RuntimeError, ImportError):
            self._backend_name = "numpy"
            self._use_numpy = True
            if verbose:
                print(
                    "  EMDenoiser (scipy Gaussian fallback) — no DL backend found."
                )

        self._is_fitted = True
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Apply the trained denoiser.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_components, n_freqs)
            Noisy impedance data.

        Returns
        -------
        X_denoised : ndarray, same shape as input
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")

        X = np.asarray(X, dtype=np.float32)
        Xn = (X - self._x_mean) / self._x_std

        if self._use_numpy:
            try:
                from scipy.ndimage import gaussian_filter1d

                out = gaussian_filter1d(Xn, sigma=1.5, axis=-1)
            except ImportError:
                out = Xn.copy()
        elif self._backend_name == "tensorflow":
            Xn_tf = Xn.transpose(0, 2, 1)  # → (batch, n_freqs, n_components)
            out_tf = self._network.predict(Xn_tf, verbose=0).astype(
                np.float32
            )
            out = out_tf.transpose(
                0, 2, 1
            )  # → (batch, n_components, n_freqs)
        else:
            import torch

            dev = next(self._network.parameters()).device
            self._network.eval()
            with torch.no_grad():
                t = torch.from_numpy(Xn).to(dev)
                out = self._network(t).cpu().numpy()

        return out * self._x_std + self._x_mean

    # ─── internal training paths ──────────────────────────────────────────

    def _fit_torch(
        self,
        Xtr: np.ndarray,
        Xva: np.ndarray,
        noise_level: float,
        epochs: int,
        batch_size: int,
        lr: float,
        rng: np.random.Generator,
        verbose: bool,
    ) -> None:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        dev = resolve_device(self.device)
        self._network = _build_cae_torch(
            self.n_components, self.n_freqs, self.channels, self.dropout
        ).to(dev)

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=7, min_lr=1e-6
        )
        mse = nn.MSELoss()
        per_ch_std = Xtr.std(axis=(0, 2), keepdims=True)

        def _make_noisy(batch: np.ndarray) -> torch.Tensor:
            noise = rng.normal(0.0, noise_level, batch.shape).astype(
                np.float32
            )
            noise *= per_ch_std
            return torch.from_numpy(batch + noise).to(dev)

        tr_ds = TensorDataset(torch.from_numpy(Xtr))
        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []

        for ep in range(1, epochs + 1):
            self._network.train()
            ep_loss = 0.0
            for (xb,) in DataLoader(
                tr_ds, batch_size=batch_size, shuffle=True
            ):
                xb = xb.to(dev)
                noisy = _make_noisy(xb.cpu().numpy())
                pred = self._network(noisy)
                loss = mse(pred, xb)
                opt.zero_grad()
                loss.backward()
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= len(Xtr)

            self._network.eval()
            with torch.no_grad():
                vb = torch.from_numpy(Xva).to(dev)
                v_pred = self._network(_make_noisy(Xva))
                v_loss = mse(v_pred, vb).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val:
                best_val = v_loss
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 10) == 0 or ep == 1):
                print(
                    f"  Epoch {ep:>4d}/{epochs}  "
                    f"train={ep_loss:.5f}  val={v_loss:.5f}"
                )

        if best_state is not None:
            self._network.load_state_dict(best_state)
        self._history = {"train_loss": train_losses, "val_loss": val_losses}

    def _fit_tensorflow(
        self,
        Xtr: np.ndarray,
        Xva: np.ndarray,
        noise_level: float,
        epochs: int,
        batch_size: int,
        lr: float,
        rng: np.random.Generator,
        verbose: bool,
    ) -> None:
        import tensorflow as tf

        # Keras Conv1D is channels-last: (batch, n_freqs, n_components)
        Xtr_tf = Xtr.transpose(0, 2, 1)
        Xva_tf = Xva.transpose(0, 2, 1)

        per_ch_std = Xtr.std(axis=(0, 2), keepdims=True)  # (1, n_comp, 1)

        def _add_noise(X_chlast: np.ndarray) -> np.ndarray:
            noise = rng.normal(0.0, noise_level, X_chlast.shape).astype(
                np.float32
            )
            # per_ch_std is (1, n_comp, 1); transpose to (1, 1, n_comp) for channels-last
            noise *= per_ch_std.transpose(0, 2, 1)
            return X_chlast + noise

        dev = resolve_device(self.device)
        with tf.device(dev):
            self._network = _build_cae_tf(
                self.n_components, self.n_freqs, self.channels, self.dropout
            )
            self._network.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
                loss="mse",
            )
            Xtr_noisy = _add_noise(Xtr_tf)
            hist = self._network.fit(
                Xtr_noisy,
                Xtr_tf,
                validation_data=(_add_noise(Xva_tf), Xva_tf),
                epochs=epochs,
                batch_size=batch_size,
                callbacks=[
                    tf.keras.callbacks.EarlyStopping(
                        monitor="val_loss",
                        patience=10,
                        restore_best_weights=True,
                        min_delta=1e-6,
                    ),
                    tf.keras.callbacks.ReduceLROnPlateau(
                        monitor="val_loss",
                        factor=0.5,
                        patience=7,
                        min_lr=1e-6,
                    ),
                ],
                verbose=1 if verbose else 0,
            )
        self._history = {
            "train_loss": hist.history["loss"],
            "val_loss": hist.history["val_loss"],
        }

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "n_freqs": self.n_freqs,  # may still be None if save() called before fit()
            "n_components": self.n_components,
            "channels": list(self.channels),
            "dropout": self.dropout,
            "device": self.device,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {}
        if self._x_mean is not None:
            out["_x_mean"] = self._x_mean
            out["_x_std"] = self._x_std
        if self._backend_name is not None:
            out["_backend"] = np.array(self._backend_name)
        if self._network is not None:
            for k, v in get_weights(self._network).items():
                out[k] = v
        return out

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        self._x_mean = weights.pop("_x_mean", None)
        self._x_std = weights.pop("_x_std", None)
        backend_blob = weights.pop("_backend", None)
        self._backend_name = (
            str(backend_blob) if backend_blob is not None else "torch"
        )

        if self._backend_name == "tensorflow":
            self._network = _build_cae_tf(
                self.n_components, self.n_freqs, self.channels, self.dropout
            )
        else:
            self._network = _build_cae_torch(
                self.n_components, self.n_freqs, self.channels, self.dropout
            )
        set_weights(self._network, weights)
        self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"EMDenoiser(n_freqs={self.n_freqs}, "
            f"n_components={self.n_components}, {status})"
        )

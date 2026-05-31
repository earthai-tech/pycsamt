# -*- coding: utf-8 -*-
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

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .._base import BaseEMProcessor, EMCheckpoint

__all__ = ["EMDenoiser", "prepare_z_features"]


# ─────────────────────────────────────────────────────────────────────────────
# Input preparation helper (no torch dependency)
# ─────────────────────────────────────────────────────────────────────────────

def prepare_z_features(
    sites: Any,
    *,
    n_components: int = 4,
    log_amp: bool = True,
    freq_ref: Optional[np.ndarray] = None,
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
        from pycsamt.emtools._core import ensure_sites, _iter_items, _get_z_block
    except ImportError as exc:
        raise ImportError("emtools is required for prepare_z_features") from exc

    S = ensure_sites(sites, recursive=True, on_dup="replace")
    rows: List[np.ndarray] = []
    freq_grid: Optional[np.ndarray] = freq_ref

    for i, ed in enumerate(_iter_items(S)):
        result = _get_z_block(ed, with_errors=False)
        if len(result) == 3:
            _, z, fr = result
        else:
            _, z, fr, _ = result
        if z is None:
            continue

        if freq_grid is None:
            freq_grid = fr

        log10e = np.log10(np.e)  # unused, use log10 directly
        comps: List[np.ndarray] = []

        def _amp(z_comp: np.ndarray) -> np.ndarray:
            amp = np.abs(z_comp)
            if log_amp:
                return np.log10(np.maximum(amp, 1e-24))
            return amp

        def _phase(z_comp: np.ndarray) -> np.ndarray:
            return np.degrees(np.angle(z_comp))

        # off-diagonal (always included)
        comps += [_amp(z[:, 0, 1]), _phase(z[:, 0, 1])]  # Zxy
        comps += [_amp(z[:, 1, 0]), _phase(z[:, 1, 0])]  # Zyx

        if n_components == 8:
            comps += [_amp(z[:, 0, 0]), _phase(z[:, 0, 0])]  # Zxx
            comps += [_amp(z[:, 1, 1]), _phase(z[:, 1, 1])]  # Zyy

        mat = np.stack(comps, axis=0)  # (n_comp, n_freqs_site)

        # Interpolate to common freq grid
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

    return np.stack(rows, axis=0).astype(np.float32)  # (n_sites, n_comp, n_freqs)


# ─────────────────────────────────────────────────────────────────────────────
# Internal network builder (lazy torch import)
# ─────────────────────────────────────────────────────────────────────────────

def _build_cae(
    n_components: int,
    n_freqs: int,
    channels: Tuple[int, ...],
    dropout: float,
) -> Any:
    """Build a 1-D convolutional autoencoder (CAE) module."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError("PyTorch is required for EMDenoiser") from exc

    mid = max(n_freqs // 4, 4)

    class _CAE(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            ch = list(channels)
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


# ─────────────────────────────────────────────────────────────────────────────
# EMDenoiser
# ─────────────────────────────────────────────────────────────────────────────

class EMDenoiser(BaseEMProcessor):
    """
    Convolutional autoencoder for MT impedance tensor denoising.

    The network is trained on clean (or synthetic) impedance data by
    adding controlled noise and learning to reconstruct the original
    signal.  It can then be applied to field data to suppress noise.

    Parameters
    ----------
    n_freqs : int
        Number of frequency channels in the input.
    n_components : {4, 8}, default 4
        Feature channels per frequency: 4 = off-diagonal :math:`Z_{xy}`,
        :math:`Z_{yx}` amplitudes and phases; 8 = all four tensor
        components.
    channels : tuple of int, default (64, 128, 64)
        Channel widths for the encoder/decoder stages.
    dropout : float, default 0.1
        Dropout probability in the encoder bottleneck.
    device : str or None
        ``"cpu"``, ``"cuda"``, ``"mps"``, or ``None`` for auto-detect.

    Notes
    -----
    Call :func:`prepare_z_features` to convert a site collection to the
    ``(n_samples, n_components, n_freqs)`` array expected here.

    Examples
    --------
    >>> import numpy as np
    >>> X_clean = np.random.randn(200, 4, 32).astype("float32")
    >>> den = EMDenoiser(n_freqs=32, n_components=4)
    >>> den.fit(X_clean, epochs=5, verbose=False)   # doctest: +SKIP
    EMDenoiser(n_freqs=32, n_components=4)
    >>> X_denoised = den.transform(X_clean)          # doctest: +SKIP
    """

    def __init__(
        self,
        n_freqs: int,
        n_components: int = 4,
        *,
        channels: Tuple[int, ...] = (64, 128, 64),
        dropout: float = 0.1,
        device: Optional[str] = None,
    ) -> None:
        self.n_freqs = int(n_freqs)
        self.n_components = int(n_components)
        self.channels = tuple(channels)
        self.dropout = float(dropout)
        self.device = device

        self._network: Any = None
        self._is_fitted: bool = False
        self._x_mean: Optional[np.ndarray] = None
        self._x_std: Optional[np.ndarray] = None
        self._history: Dict[str, list] = {}

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
        seed: Optional[int] = None,
        verbose: bool = True,
    ) -> "EMDenoiser":
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
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import TensorDataset, DataLoader
        except ImportError as exc:
            raise ImportError("PyTorch is required for EMDenoiser.fit") from exc

        rng = np.random.default_rng(seed)
        X = np.asarray(X, dtype=np.float32)

        # normalise per channel (axis 0 = samples)
        self._x_mean = X.mean(axis=(0, 2), keepdims=True)
        self._x_std = X.std(axis=(0, 2), keepdims=True) + 1e-8
        Xn = (X - self._x_mean) / self._x_std

        # train / val split
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx, trn_idx = idx[:n_val], idx[n_val:]
        Xtr, Xva = Xn[trn_idx], Xn[val_idx]

        dev = self._resolve_device()
        self._network = _build_cae(
            self.n_components, self.n_freqs, self.channels, self.dropout
        ).to(dev)

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=7, min_lr=1e-6
        )
        mse = nn.MSELoss()

        per_ch_std = Xtr.std(axis=(0, 2), keepdims=True)

        def _make_noisy(batch: np.ndarray) -> "torch.Tensor":
            noise = rng.normal(0.0, noise_level, batch.shape).astype(np.float32)
            noise *= per_ch_std  # scale noise per channel
            return torch.from_numpy(batch + noise).to(dev)

        tr_ds = TensorDataset(torch.from_numpy(Xtr))
        va_ds = TensorDataset(torch.from_numpy(Xva))

        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []

        for ep in range(1, epochs + 1):
            self._network.train()
            loader = DataLoader(tr_ds, batch_size=batch_size, shuffle=True)
            ep_loss = 0.0
            for (xb,) in loader:
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
                import copy
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 10) == 0 or ep == 1):
                print(
                    f"  Epoch {ep:>4d}/{epochs}  "
                    f"train={ep_loss:.5f}  val={v_loss:.5f}"
                )

        if best_state is not None:
            self._network.load_state_dict(best_state)

        self._history = {"train_loss": train_losses, "val_loss": val_losses}
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

        try:
            import torch
        except ImportError as exc:
            raise ImportError("PyTorch is required for EMDenoiser.transform") from exc

        X = np.asarray(X, dtype=np.float32)
        Xn = (X - self._x_mean) / self._x_std
        dev = next(self._network.parameters()).device
        self._network.eval()
        with torch.no_grad():
            t = torch.from_numpy(Xn).to(dev)
            out = self._network(t).cpu().numpy()
        return out * self._x_std + self._x_mean

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> Dict[str, Any]:
        return {
            "n_freqs": self.n_freqs,
            "n_components": self.n_components,
            "channels": list(self.channels),
            "dropout": self.dropout,
            "device": self.device,
        }

    def _get_weights(self) -> Dict[str, np.ndarray]:
        if self._network is None:
            return {}
        weights = {
            k: v.cpu().numpy()
            for k, v in self._network.state_dict().items()
        }
        if self._x_mean is not None:
            weights["_x_mean"] = self._x_mean
            weights["_x_std"] = self._x_std
        return weights

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        import torch
        self._network = _build_cae(
            self.n_components, self.n_freqs, self.channels, self.dropout
        )
        self._x_mean = weights.pop("_x_mean", None)
        self._x_std = weights.pop("_x_std", None)
        state = {k: torch.from_numpy(v) for k, v in weights.items()}
        self._network.load_state_dict(state)
        self._is_fitted = True

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

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"EMDenoiser(n_freqs={self.n_freqs}, "
            f"n_components={self.n_components}, {status})"
        )

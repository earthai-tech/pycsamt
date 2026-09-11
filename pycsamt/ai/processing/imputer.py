# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
EMImputer — masked-reconstruction gap filling for MT impedance data.

Why this is not EMDenoiser
---------------------------
:class:`~pycsamt.ai.processing.denoise.EMDenoiser` assumes every
(station, frequency) cell has *some* value and improves its quality.
Real surveys also have cells with *no* value at all -- a frequency
dropped mid-acquisition, a component that never converged, a station
skipped for one band -- and nothing else in pyCSAMT fills those gaps
with anything better than leaving them ``NaN`` through to inversion.
``EMImputer`` targets that gap specifically: reconstructing the
*missing* cells from the *observed* ones, using spectral smoothness
(this station, nearby frequencies) -- a matrix-completion /
masked-autoencoder problem, not a denoising one.

Architecture
------------
A self-supervised masked-reconstruction objective, the same family as
BERT-style masked-language-modelling or image inpainting:

.. math::

    \\hat{\\mathbf{x}} =
    f_\\theta(\\mathbf{x} \\odot \\mathbf{m},\\ \\mathbf{m}),
    \\qquad
    \\mathcal{L} = \\frac{\\big\\| (\\hat{\\mathbf{x}} - \\mathbf{x}) \\odot
    \\mathbf{m}_\\text{syn} \\big\\|_2^2}{\\sum \\mathbf{m}_\\text{syn}}

where :math:`\\mathbf{m} \\in \\{0, 1\\}` is the mask fed to the network
(1 = "the network is told this cell is observed") and
:math:`\\mathbf{m}_\\text{syn}` marks the subset of *genuinely observed*
cells that were synthetically hidden this step, which is exactly where
the loss is evaluated. On every training step, :meth:`EMImputer.fit`
hides a random ``mask_frac`` fraction of the truly-observed cells (any
real ``NaN`` gaps in the training data are always hidden and never
used as a target, since their true value is unknown) -- mirroring the
way :class:`~pycsamt.ai.processing.denoise.EMDenoiser` adds synthetic
noise, so training data only needs to be *mostly* complete rather than
needing real missing-data examples with known ground truth.

The network is the same 1-D convolutional encoder-decoder family as
``EMDenoiser`` (:func:`~pycsamt.ai.processing.denoise._build_cae_torch`
/ ``_build_cae_tf``), with one structural difference: the encoder
takes ``2 * n_components`` input channels (the zeroed-where-hidden
values, concatenated with the mask along the channel axis) so the
network can distinguish "zero because hidden" from "zero because that
is genuinely the value", while the decoder still only reconstructs
``n_components`` channels -- the mask itself is never reconstructed.

Input tensor shape
-------------------
``(n_samples, n_components, n_freqs)``, identical to
:func:`~pycsamt.ai.processing.denoise.prepare_z_features`'s output --
``NaN`` marks a missing cell.

Working directly with site collections
---------------------------------------
:meth:`EMImputer.apply` accepts and returns a site collection
directly, writing the reconstructed impedance back only into cells
that were genuinely missing on the input -- every observed measurement
is left byte-identical -- the same sites-in / sites-out convention as
:meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply`.
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
from .denoise import _reconstruct_z_block, prepare_z_features

__all__ = ["EMImputer"]


# ─────────────────────────────────────────────────────────────────────────────
# Fallback: frequency-axis interpolation (no PyTorch / TensorFlow)
# ─────────────────────────────────────────────────────────────────────────────


def _freq_axis_fill(X: np.ndarray, x_mean: np.ndarray) -> np.ndarray:
    """
    Per-sample, per-channel linear interpolation along the frequency
    axis, using only that row's finite values.

    Used when neither PyTorch nor TensorFlow is available. A channel
    with no finite value at all falls back to that channel's training
    mean (``x_mean``) -- there is nothing else to interpolate from.
    """
    out = X.copy()
    n_samples, n_comp, n_freq = X.shape
    idx = np.arange(n_freq, dtype=float)
    for s in range(n_samples):
        for c in range(n_comp):
            row = out[s, c]
            finite = np.isfinite(row)
            n_finite = int(finite.sum())
            if n_finite == 0:
                row[:] = float(x_mean[0, c, 0])
            elif n_finite < n_freq:
                row[~finite] = np.interp(
                    idx[~finite], idx[finite], row[finite]
                )
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Network builders — PyTorch (channels-first) and TensorFlow (channels-last)
# ─────────────────────────────────────────────────────────────────────────────


def _build_masked_cae_torch(
    n_components: int,
    n_freqs: int,
    channels: tuple[int, ...],
    dropout: float,
) -> Any:
    """
    1-D masked CAE using PyTorch Conv1d (channels-first).

    Encoder input is ``2 * n_components`` channels (zeroed-where-hidden
    values concatenated with the observed-cell mask); the decoder
    output is ``n_components`` channels -- the reconstruction only.
    """
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError("PyTorch is required for EMImputer") from exc

    mid = max(n_freqs // 4, 4)
    ch = list(channels)
    in_ch = 2 * n_components

    class _MaskedCAE(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Conv1d(in_ch, ch[0], 5, padding=2),
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

        def forward(self, x):  # x: (B, 2*n_comp, n_freqs)
            return self.decoder(self.encoder(x))

    return _MaskedCAE()


def _build_masked_cae_tf(
    n_components: int,
    n_freqs: int,
    channels: tuple[int, ...],
    dropout: float,
) -> Any:
    """
    1-D masked CAE using Keras Conv1D (channels-last).

    Input / output convention: ``(batch, n_freqs, 2 * n_components)``
    in, ``(batch, n_freqs, n_components)`` out. Data is transposed by
    the caller from/to the canonical ``(batch, channels, n_freqs)``
    form, the same convention as
    :func:`~pycsamt.ai.processing.denoise._build_cae_tf`.
    """
    try:
        import tensorflow as tf
        from tensorflow.keras import Model, layers
    except ImportError as exc:
        raise ImportError(
            "TensorFlow is required for EMImputer (TF backend)"
        ) from exc

    ch = list(channels)
    pool = max(1, n_freqs // max(n_freqs // 4, 4))
    in_ch = 2 * n_components

    inp = tf.keras.Input(shape=(n_freqs, in_ch), name="input")

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

    return Model(inp, out, name="masked_cae")


# ─────────────────────────────────────────────────────────────────────────────
# EMImputer
# ─────────────────────────────────────────────────────────────────────────────


class EMImputer(BaseEMProcessor):
    """
    Masked-reconstruction gap filler for MT impedance tensor data.

    Trains a 1-D convolutional autoencoder to reconstruct
    synthetically-hidden cells of otherwise-observed impedance
    spectra, then applies that network to fill genuinely missing
    (``NaN``) cells at inference time. Supports both PyTorch and
    TensorFlow backends via :mod:`pycsamt.backends`.

    Parameters
    ----------
    n_freqs : int or None, default None
        Number of frequency channels in the input. When ``None``,
        inferred from ``X.shape[2]`` on the first call to :meth:`fit`.
    n_components : {4, 8}, default 4
        Feature channels per frequency, matching
        :func:`~pycsamt.ai.processing.denoise.prepare_z_features`.
    channels : tuple of int, default (64, 128, 64)
        Channel widths for the encoder/decoder stages.
    dropout : float, default 0.1
        Dropout probability in the encoder bottleneck.
    device : str or None
        ``"cpu"``, ``"cuda"``, ``"mps"``, ``"/GPU:0"``, or ``None``
        for auto-detect via the active backend.

    Notes
    -----
    When neither PyTorch nor TensorFlow is available the imputer falls
    back to per-sample, per-channel linear interpolation along the
    frequency axis (:func:`_freq_axis_fill`) -- no spatial (cross-
    station) information, unlike the network path.

    Call :func:`~pycsamt.ai.processing.denoise.prepare_z_features` to
    convert a site collection to the
    ``(n_samples, n_components, n_freqs)`` array expected here, or use
    :meth:`apply` to go straight from a site collection to a filled
    one without handling the array form yourself.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.random.randn(200, 4, 32).astype("float32")
    >>> rng = np.random.default_rng(0)
    >>> gaps = rng.random(X.shape) < 0.05
    >>> X[gaps] = np.nan
    >>> imp = EMImputer()
    >>> imp.fit(X, epochs=5, verbose=False)  # doctest: +SKIP
    EMImputer(n_freqs=32, n_components=4, fitted)
    >>> X_filled = imp.transform(X)  # doctest: +SKIP

    Sites in, filled sites out -- the same convention as
    :meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply`:

    >>> complete_sites = imp.apply(sites)  # doctest: +SKIP
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
        mask_frac: float = 0.15,
        epochs: int = 80,
        batch_size: int = 64,
        lr: float = 1e-3,
        val_frac: float = 0.1,
        seed: int | None = None,
        verbose: bool = True,
    ) -> EMImputer:
        """
        Train the imputer on a (mostly) complete impedance array.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_components, n_freqs)
            Training data. ``NaN`` cells, if any, are excluded from
            the reconstruction loss -- only cells that are genuinely
            observed can be synthetically hidden and used as a
            self-supervised target.
        mask_frac : float, default 0.15
            Fraction of the *observed* cells to hide on every training
            step, supplying the reconstruction target (see the module
            docstring). Too low and the network sees little missing-
            data signal to learn from; too high and too little context
            remains for it to reconstruct anything.
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

        Raises
        ------
        ValueError
            If ``X`` is not 3-D, its frequency count does not match a
            previously-fitted ``n_freqs``, or it contains no finite
            values at all.
        """
        rng = np.random.default_rng(seed)
        X = np.asarray(X, dtype=np.float32)

        if X.ndim != 3:
            raise ValueError(
                "X must be 3-D (n_samples, n_components, n_freqs), "
                f"got shape {X.shape}"
            )

        if self.n_freqs is None:
            self.n_freqs = X.shape[2]
        elif X.shape[2] != self.n_freqs:
            raise ValueError(
                f"Expected n_freqs={self.n_freqs}, got {X.shape[2]}"
            )

        observed = np.isfinite(X)
        if not observed.any():
            raise ValueError(
                "No finite values found in X -- EMImputer has nothing "
                "to learn from."
            )

        masked_X = np.where(observed, X, np.nan)
        with np.errstate(invalid="ignore"):
            x_mean = np.nanmean(masked_X, axis=(0, 2), keepdims=True)
            x_std = np.nanstd(masked_X, axis=(0, 2), keepdims=True)
        # A channel with zero observed values anywhere in the batch has
        # no mean/std to compute; fall back to an identity normalisation
        # for it rather than propagating NaN into every weight update.
        x_mean = np.where(np.isfinite(x_mean), x_mean, 0.0).astype(
            np.float32
        )
        x_std = np.where(
            np.isfinite(x_std) & (x_std > 1e-8), x_std, 1.0
        ).astype(np.float32)
        self._x_mean, self._x_std = x_mean, x_std

        Xn = np.where(observed, (X - x_mean) / x_std, 0.0).astype(
            np.float32
        )

        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        val_idx, trn_idx = idx[:n_val], idx[n_val:]
        Xtr, Xva = Xn[trn_idx], Xn[val_idx]
        Otr, Ova = observed[trn_idx], observed[val_idx]

        try:
            self._backend_name = active_backend()
            if self._backend_name == "tensorflow":
                self._fit_tensorflow(
                    Xtr,
                    Xva,
                    Otr,
                    Ova,
                    mask_frac,
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
                    Otr,
                    Ova,
                    mask_frac,
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
                    "  EMImputer (frequency-interpolation fallback) "
                    "— no DL backend found."
                )

        self._is_fitted = True
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Fill missing (``NaN``) cells in ``X``.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_components, n_freqs)
            Data with ``NaN`` marking missing cells.

        Returns
        -------
        X_filled : ndarray, same shape as input
            Observed cells returned unchanged (byte-identical to the
            input); missing cells replaced by the network's
            reconstruction.
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")

        X = np.asarray(X, dtype=np.float32)
        missing = ~np.isfinite(X)

        if self._use_numpy:
            return _freq_axis_fill(X, self._x_mean)

        Xn = np.where(
            np.isfinite(X), (X - self._x_mean) / self._x_std, 0.0
        ).astype(np.float32)
        mask = np.isfinite(X).astype(np.float32)
        net_in = np.concatenate([Xn, mask], axis=1)  # (B, 2*n_comp, F)

        if self._backend_name == "tensorflow":
            net_in_tf = net_in.transpose(0, 2, 1)  # → (B, F, 2*n_comp)
            out_tf = self._network.predict(net_in_tf, verbose=0).astype(
                np.float32
            )
            out = out_tf.transpose(0, 2, 1)  # → (B, n_comp, F)
        else:
            import torch

            dev = next(self._network.parameters()).device
            self._network.eval()
            with torch.no_grad():
                t = torch.from_numpy(net_in).to(dev)
                out = self._network(t).cpu().numpy()

        pred = out * self._x_std + self._x_mean
        X_filled = X.copy()
        X_filled[missing] = pred[missing]
        return X_filled

    # ─── sites-in / sites-out ──────────────────────────────────────────────

    def apply(
        self,
        sites: Any,
        *,
        freq_ref: np.ndarray | None = None,
        log_amp: bool = True,
        inplace: bool = False,
        recursive: bool = True,
        on_dup: str = "replace",
        strict: bool = False,
        verbose: int = 0,
    ) -> Any:
        """
        Fill a site collection's missing impedance cells in place.

        Extracts features with
        :func:`~pycsamt.ai.processing.denoise.prepare_z_features`,
        imputes them with :meth:`transform`, then writes the
        reconstruction back onto each site's own frequency grid --
        but only at frequency rows that were genuinely missing on that
        site to begin with; every observed row is left byte-identical
        to the input. Mirrors
        :meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply`'s
        sites-in / sites-out convention.

        Parameters
        ----------
        sites : SiteCollection or compatible
        freq_ref : ndarray or None
            Common frequency grid (Hz) the network was trained on.
            Must match the grid used to prepare the training data --
            when ``None`` (the default), the first site's own grid is
            used, matching
            :func:`~pycsamt.ai.processing.denoise.prepare_z_features`'s
            own default.
        log_amp : bool, default ``True``
            Must match the ``log_amp`` used to build the training
            features.
        inplace : bool, default ``False``
            When ``False`` (the default), *sites* is left untouched
            and a filled copy is returned. When ``True``, *sites* is
            modified in place and returned.
        recursive, on_dup, strict, verbose
            Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

        Returns
        -------
        filled : Sites
            A site collection with the same stations. Stations with no
            missing cells at all are returned untouched.

        Warns
        -----
        With the default ``n_components=4``, only :math:`Z_{xy}` and
        :math:`Z_{yx}` are modelled -- a real missing measurement
        usually means the *whole* frequency row (all four components)
        is absent, and the diagonal :math:`Z_{xx}, Z_{yy}` at those
        rows is left exactly as it was (``NaN``, if it was ``NaN``),
        the same "not modelled -> untouched" contract
        :meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply` uses.
        Fit and apply with ``n_components=8`` to fill every component.
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before apply().")

        try:
            from pycsamt.emtools._core import (
                _apply_each,
                _get_z_block,
                _iter_items,
                ensure_sites,
            )
        except ImportError as exc:
            raise ImportError(
                "emtools is required for EMImputer.apply()"
            ) from exc

        S = ensure_sites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )

        if freq_ref is None:
            first_ed = next(_iter_items(S))
            _, _, freq_ref = _get_z_block(first_ed, with_errors=False)[:3]

        X = prepare_z_features(
            S,
            n_components=self.n_components,
            log_amp=log_amp,
            freq_ref=freq_ref,
        )
        X_filled = self.transform(X)

        site_idx = iter(range(len(X_filled)))

        def _one(Si):
            i = next(site_idx)
            ed = next(_iter_items(Si))
            Z, z, fr = _get_z_block(ed, with_errors=False)[:3]
            if z is None:
                return Si
            missing = np.isnan(z).any(axis=(1, 2))
            if not missing.any():
                return Si
            z_recon = _reconstruct_z_block(
                z, fr, freq_ref, X_filled[i], self.n_components, log_amp
            )
            z2 = z.copy()
            z2[missing] = z_recon[missing]
            Z.z = z2
            return Si

        return _apply_each(S, _one, inplace=inplace, verbose=verbose)

    # ─── internal training paths ──────────────────────────────────────────

    def _fit_torch(
        self,
        Xtr: np.ndarray,
        Xva: np.ndarray,
        Otr: np.ndarray,
        Ova: np.ndarray,
        mask_frac: float,
        epochs: int,
        batch_size: int,
        lr: float,
        rng: np.random.Generator,
        verbose: bool,
    ) -> None:
        import torch
        from torch.utils.data import DataLoader, TensorDataset

        dev = resolve_device(self.device)
        self._network = _build_masked_cae_torch(
            self.n_components, self.n_freqs, self.channels, self.dropout
        ).to(dev)

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=7, min_lr=1e-6
        )

        def _make_masked(batch_np: np.ndarray, observed_np: np.ndarray):
            syn = (rng.random(observed_np.shape) < mask_frac) & observed_np
            input_mask = observed_np & ~syn
            x_in = np.where(input_mask, batch_np, 0.0).astype(np.float32)
            mask_ch = input_mask.astype(np.float32)
            net_in = np.concatenate([x_in, mask_ch], axis=1)
            loss_mask = syn.astype(np.float32)
            return (
                torch.from_numpy(net_in).to(dev),
                torch.from_numpy(loss_mask).to(dev),
            )

        def _masked_mse(pred, target_t, loss_mask_t):
            num = ((pred - target_t) ** 2 * loss_mask_t).sum()
            den = loss_mask_t.sum().clamp_min(1.0)
            return num / den

        tr_ds = TensorDataset(
            torch.from_numpy(Xtr), torch.from_numpy(Otr.astype(np.float32))
        )
        # A fixed validation synthetic mask, drawn once, keeps the
        # reported val loss comparable across epochs (re-drawing it
        # every epoch would add its own noise to the early-stopping
        # signal below).
        Xva_in, Xva_loss_mask = _make_masked(Xva, Ova)
        Xva_t = torch.from_numpy(Xva).to(dev)

        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []

        for ep in range(1, epochs + 1):
            self._network.train()
            ep_loss, ep_n = 0.0, 0
            for xb, ob in DataLoader(
                tr_ds, batch_size=batch_size, shuffle=True
            ):
                xb_np = xb.numpy()
                ob_np = ob.numpy().astype(bool)
                net_in, loss_mask = _make_masked(xb_np, ob_np)
                target = xb.to(dev)
                pred = self._network(net_in)
                loss = _masked_mse(pred, target, loss_mask)
                opt.zero_grad()
                loss.backward()
                opt.step()
                ep_loss += loss.item() * len(xb)
                ep_n += len(xb)
            ep_loss /= max(ep_n, 1)

            self._network.eval()
            with torch.no_grad():
                v_pred = self._network(Xva_in)
                v_loss = _masked_mse(v_pred, Xva_t, Xva_loss_mask).item()

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
        Otr: np.ndarray,
        Ova: np.ndarray,
        mask_frac: float,
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
        Otr_tf = Otr.transpose(0, 2, 1)
        Ova_tf = Ova.transpose(0, 2, 1)

        def _make_masked(x_np: np.ndarray, o_np: np.ndarray):
            syn = (rng.random(o_np.shape) < mask_frac) & o_np
            input_mask = o_np & ~syn
            x_in = np.where(input_mask, x_np, 0.0).astype(np.float32)
            mask_ch = input_mask.astype(np.float32)
            net_in = np.concatenate([x_in, mask_ch], axis=-1)
            loss_mask = syn.astype(np.float32)
            return net_in, loss_mask

        def _masked_mse(pred, target, loss_mask):
            num = tf.reduce_sum(((pred - target) ** 2) * loss_mask)
            den = tf.maximum(tf.reduce_sum(loss_mask), 1.0)
            return num / den

        dev = resolve_device(self.device)
        train_losses: list[float] = []
        val_losses: list[float] = []

        with tf.device(dev):
            self._network = _build_masked_cae_tf(
                self.n_components, self.n_freqs, self.channels, self.dropout
            )
            opt = tf.keras.optimizers.Adam(learning_rate=lr)

            Xva_in, Xva_loss_mask = _make_masked(Xva_tf, Ova_tf)
            Xva_in_t = tf.constant(Xva_in)
            Xva_t = tf.constant(Xva_tf)
            Xva_loss_mask_t = tf.constant(Xva_loss_mask)

            n = len(Xtr_tf)
            best_val = np.inf
            best_weights: list[np.ndarray] | None = None
            patience, bad_epochs = 7, 0

            for ep in range(1, epochs + 1):
                perm = rng.permutation(n)
                ep_loss, ep_n = 0.0, 0
                for start in range(0, n, batch_size):
                    bi = perm[start : start + batch_size]
                    xb, ob = Xtr_tf[bi], Otr_tf[bi]
                    net_in, loss_mask = _make_masked(xb, ob)
                    net_in_t = tf.constant(net_in)
                    target_t = tf.constant(xb)
                    loss_mask_t = tf.constant(loss_mask)
                    with tf.GradientTape() as tape:
                        pred = self._network(net_in_t, training=True)
                        loss = _masked_mse(pred, target_t, loss_mask_t)
                    grads = tape.gradient(
                        loss, self._network.trainable_variables
                    )
                    opt.apply_gradients(
                        zip(grads, self._network.trainable_variables)
                    )
                    ep_loss += float(loss) * len(bi)
                    ep_n += len(bi)
                ep_loss /= max(ep_n, 1)

                v_pred = self._network(Xva_in_t, training=False)
                v_loss = float(
                    _masked_mse(v_pred, Xva_t, Xva_loss_mask_t)
                )

                train_losses.append(ep_loss)
                val_losses.append(v_loss)

                if v_loss < best_val - 1e-6:
                    best_val = v_loss
                    best_weights = [
                        w.numpy() for w in self._network.weights
                    ]
                    bad_epochs = 0
                else:
                    bad_epochs += 1
                    if bad_epochs >= patience:
                        opt.learning_rate.assign(
                            float(opt.learning_rate) * 0.5
                        )
                        bad_epochs = 0

                if verbose and (
                    ep % max(1, epochs // 10) == 0 or ep == 1
                ):
                    print(
                        f"  Epoch {ep:>4d}/{epochs}  "
                        f"train={ep_loss:.5f}  val={v_loss:.5f}"
                    )

            if best_weights is not None:
                for w, val in zip(self._network.weights, best_weights):
                    w.assign(val)

        self._history = {"train_loss": train_losses, "val_loss": val_losses}

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "n_freqs": self.n_freqs,
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
            self._network = _build_masked_cae_tf(
                self.n_components, self.n_freqs, self.channels, self.dropout
            )
        else:
            self._network = _build_masked_cae_torch(
                self.n_components, self.n_freqs, self.channels, self.dropout
            )
        set_weights(self._network, weights)
        self._is_fitted = True

    @property
    def history_(self) -> dict[str, list]:
        """
        Training history recorded by the last :meth:`fit` call.

        Returns
        -------
        history : dict
            ``{"train_loss": [...], "val_loss": [...]}``, one value per
            epoch, both computed only over synthetically-hidden cells
            (see the module docstring). Empty when the frequency-
            interpolation fallback was used (no epoch loop) or before
            :meth:`fit` has been called. Pass directly to
            :func:`~pycsamt.ai.processing.plot.plot_training_history`.
        """
        return dict(self._history)

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"EMImputer(n_freqs={self.n_freqs}, "
            f"n_components={self.n_components}, {status})"
        )

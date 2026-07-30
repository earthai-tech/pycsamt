# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
EMInverter2D — high-level U-Net–based 2-D MT inversion pipeline.

Wraps :class:`~pycsamt.ai.nets.unet.UNet2DNet` with a complete
data-loading → normalisation → training → prediction workflow.

Input / output convention
-------------------------
* **Input** ``X``:
  ``ndarray (n_profiles, n_components, n_freqs, n_stations)`` — each
  entry is one profile of MT apparent-resistivity / phase maps.

* **Target** ``y``:
  ``ndarray (n_profiles, n_depth, n_stations)`` — 2-D log₁₀(ρ) sections.

Both ``n_freqs`` and ``n_depth`` may differ (the U-Net uses bilinear
upsampling); ``n_stations`` must match between ``X`` and ``y``.

Synthetic training data
-----------------------
When PyTorch is available and a 2-D pre-built dataset is not on hand,
use :func:`~pycsamt.forward.batch.generate_dataset` to build
pseudo-2-D profiles by running the 1-D MT forward solver independently
at each virtual station, then stack the per-station feature vectors into
the 2-D panel format expected here.

Example
-------
>>> from pycsamt.ai.inversion import EMInverter2D
>>> inv = EMInverter2D(n_components=4, n_depth=40, n_stations=20, n_freqs=32)
>>> inv.fit(X_train, y_train, epochs=30)  # doctest: +SKIP
EMInverter2D(arch='unet', fitted)
>>> rho_pred = inv.predict(X_test)  # doctest: +SKIP
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
from .._base import BaseEMNet

__all__ = ["EMInverter2D"]


def _staged_torch_loss(
    pred, target, *, lambda_x=0.0, lambda_z=0.0, lambda_tv=0.0
):
    """Mean-squared data fit plus optional spatial regularization.

    ``pred``/``target`` are normalized ``(N, 1, H, W)`` torch tensors
    with ``H`` the depth-like axis and ``W`` the station/lateral
    axis. Mirrors :mod:`pycsamt.ai.losses.model`/``spatial``'s L2/L1
    math (``L_model``/``L_grad_x``/``L_grad_z``/``L_TV``) but stays
    differentiable by operating on torch tensors directly — those
    NumPy implementations cannot participate in backpropagation.
    With ``lambda_x = lambda_z = lambda_tv = 0`` this is exactly
    ``nn.MSELoss()(pred, target)``.
    """
    loss = ((pred - target) ** 2).mean()
    if lambda_x:
        grad_x = pred[..., 1:] - pred[..., :-1]
        loss = loss + lambda_x * (grad_x**2).mean()
    if lambda_z:
        grad_z = pred[..., 1:, :] - pred[..., :-1, :]
        loss = loss + lambda_z * (grad_z**2).mean()
    if lambda_tv:
        tv_x = (pred[..., 1:] - pred[..., :-1]).abs().mean()
        tv_z = (pred[..., 1:, :] - pred[..., :-1, :]).abs().mean()
        loss = loss + lambda_tv * (tv_x + tv_z)
    return loss


# ─────────────────────────────────────────────────────────────────────────────
# EMInverter2D
# ─────────────────────────────────────────────────────────────────────────────


class EMInverter2D(BaseEMNet):
    """
    U-Net–based 2-D MT inversion estimator.

    Parameters
    ----------
    n_components : int, default 4
        Number of input channels (EM data components per frequency and
        station).
    n_depth : int, default 40
        Number of depth cells in the target resistivity section.
    n_stations : int, default 20
        Number of stations along the profile.
    n_freqs : int, default 32
        Number of frequency channels in the input data.
    arch : str, default ``'unet'``
        Network architecture.  Only ``'unet'`` is supported in Phase 4.
    device : str or None
        Compute device.
    log_rho_out : bool, default True
        If ``True``, targets and predictions are in log₁₀(ρ) scale.
    **net_kwargs
        Additional keyword arguments forwarded to the architecture
        factory (e.g. ``channels``, ``dropout``).
    """

    # Base channel widths for each encoder stage (bridge appended adaptively)
    _BASE_CHANNELS: tuple[int, ...] = (32, 64, 128, 256, 512)

    def __init__(
        self,
        n_components: int = 4,
        n_depth: int = 40,
        n_stations: int = 20,
        n_freqs: int = 32,
        *,
        arch: str = "unet",
        unet_depth: int | None = None,
        channels: tuple[int, ...] | None = None,
        dropout: float = 0.2,
        device: str | None = None,
        log_rho_out: bool = True,
        **net_kwargs,
    ) -> None:
        super().__init__(
            arch=arch, n_layers=n_depth, solver="mt2d", device=device
        )
        self.n_components = int(n_components)
        self.n_depth = int(n_depth)
        self.n_stations = int(n_stations)
        self.n_freqs = int(n_freqs)
        self.log_rho_out = bool(log_rho_out)
        self.dropout = float(dropout)
        self._net_kwargs = net_kwargs

        # Compute adaptive channel spec if not supplied explicitly
        self._channels = self._resolve_channels(channels, unet_depth)

        self._x_mean: float | None = None
        self._x_std: float | None = None
        self._y_mean: float | None = None
        self._y_std: float | None = None
        self._backend_name: str | None = None

    def _resolve_channels(
        self,
        channels: tuple[int, ...] | None,
        unet_depth: int | None,
    ) -> tuple[int, ...]:
        """Return channel tuple safe for current (n_freqs, n_stations)."""
        import math

        if channels is not None:
            return tuple(channels)
        # Max safe pooling depth: each MaxPool2d(2) halves both H and W.
        # We need min_dim / 2^depth >= 1, so depth <= log2(min_dim).
        min_dim = min(self.n_freqs, self.n_stations)
        max_safe = int(math.floor(math.log2(max(min_dim, 1))))
        n_stages = max(1, min(max_safe, 4))
        if unet_depth is not None:
            n_stages = max(1, min(int(unet_depth), max_safe))
        # channels has n_stages encoder widths + 1 bridge width
        base = self._BASE_CHANNELS
        return base[:n_stages] + (base[min(n_stages, len(base) - 1)],)

    # ─── BaseEMNet interface ──────────────────────────────────────────────

    def _build_network(self) -> Any:
        from pycsamt.backends import get_backend_instance

        spec = {
            "arch": "unet2d",
            "n_in": self.n_components,
            "n_out": 1,
            "channels": self._channels,
            "dropout": self.dropout,
            **self._net_kwargs,
        }
        return get_backend_instance().build(spec)

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
        *,
        epochs: int = 100,
        batch_size: int = 16,
        lr: float = 1e-3,
        patience: int = 15,
        val_frac: float = 0.1,
        grad_clip: float | None = 1.0,
        seed: int | None = None,
        verbose: bool = True,
        lambda_x: float = 0.0,
        lambda_z: float = 0.0,
        lambda_tv: float = 0.0,
    ) -> EMInverter2D:
        """
        Train the 2-D inversion network.

        Parameters
        ----------
        X : ndarray (n_profiles, n_components, n_freqs, n_stations)
            Input data panels.
        y : ndarray (n_profiles, n_depth, n_stations)
            Target 2-D log₁₀(ρ) sections.
        epochs, batch_size, lr, patience, val_frac, grad_clip, seed, verbose
            Standard training hyper-parameters.
        lambda_x, lambda_z, lambda_tv : float, default 0.0
            Weights for lateral-gradient, depth-gradient, and total-
            variation spatial regularization terms added to the plain
            data-fit MSE, matching the staged objective in
            :mod:`pycsamt.ai.losses` (``L_grad_x``/``L_grad_z``/
            ``L_TV``) but computed directly on normalized torch
            tensors so the combined loss stays differentiable — the
            NumPy implementations in :mod:`pycsamt.ai.losses` cannot
            participate in backpropagation. All zero (the default)
            reproduces the exact plain-MSE behaviour of earlier
            versions. Only supported on the PyTorch backend: passing
            a nonzero value while TensorFlow is active raises
            :class:`NotImplementedError`.

        Returns
        -------
        self
        """
        self._backend_name = active_backend()
        if self._backend_name == "tensorflow" and (
            lambda_x or lambda_z or lambda_tv
        ):
            raise NotImplementedError(
                "lambda_x/lambda_z/lambda_tv staged regularization is "
                "only implemented for the PyTorch backend; the "
                "TensorFlow path still trains with plain MSE."
            )
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.float32)

        self._x_mean = float(np.nanmean(X))
        self._x_std = float(np.nanstd(X)) + 1e-8
        self._y_mean = float(np.nanmean(y))
        self._y_std = float(np.nanstd(y)) + 1e-8

        Xn = (X - self._x_mean) / self._x_std
        yn = (y - self._y_mean) / self._y_std

        rng = np.random.default_rng(seed)
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        vi, ti = idx[:n_val], idx[n_val:]

        self._network = self._build_network()

        if self._backend_name == "tensorflow":
            hist, best_val = self._fit_tensorflow(
                Xn[ti],
                yn[ti],
                Xn[vi],
                yn[vi],
                epochs=epochs,
                batch_size=batch_size,
                lr=lr,
                patience=patience,
                verbose=verbose,
            )
        else:
            hist, best_val = self._fit_torch(
                Xn[ti],
                yn[ti],
                Xn[vi],
                yn[vi],
                epochs=epochs,
                batch_size=batch_size,
                lr=lr,
                patience=patience,
                grad_clip=grad_clip,
                verbose=verbose,
                lambda_x=lambda_x,
                lambda_z=lambda_z,
                lambda_tv=lambda_tv,
            )

        self._history = hist
        self._meta["best_val_loss"] = float(best_val)
        self._meta["loss_weights"] = {
            "lambda_x": lambda_x,
            "lambda_z": lambda_z,
            "lambda_tv": lambda_tv,
        }
        self._is_fitted = True
        return self

    # ─── internal training paths ──────────────────────────────────────────

    @staticmethod
    def _resize(pred, target_hw):
        """Bilinear-resize *pred* to (H, W) = *target_hw* when shapes differ."""
        import torch.nn.functional as F

        if pred.shape[-2:] == target_hw:
            return pred
        return F.interpolate(
            pred, size=target_hw, mode="bilinear", align_corners=False
        )

    def _fit_torch(
        self,
        Xtr,
        ytr,
        Xva,
        yva,
        *,
        epochs,
        batch_size,
        lr,
        patience,
        grad_clip,
        verbose,
        lambda_x=0.0,
        lambda_z=0.0,
        lambda_tv=0.0,
    ):
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        # channels-first: (n, n_comp, n_freqs, n_sta) — already correct
        # target: add channel dim → (n, 1, n_depth, n_sta)
        ytr = ytr[:, np.newaxis, :, :]
        yva = yva[:, np.newaxis, :, :]
        target_hw = (ytr.shape[2], ytr.shape[3])  # (n_depth, n_stations)

        dev = resolve_device(self.device)
        self._network = self._network.to(dev)

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=max(5, patience // 3), min_lr=1e-6
        )

        def loss_fn(pred, target):
            return _staged_torch_loss(
                pred,
                target,
                lambda_x=lambda_x,
                lambda_z=lambda_z,
                lambda_tv=lambda_tv,
            )

        tr_ds = TensorDataset(torch.from_numpy(Xtr), torch.from_numpy(ytr))
        # Drop a trailing batch of size 1: BatchNorm layers cannot
        # compute per-channel statistics from a single sample.
        drop_last = len(tr_ds) > batch_size and len(tr_ds) % batch_size == 1
        Xva_t = torch.from_numpy(Xva).to(dev)
        yva_t = torch.from_numpy(yva).to(dev)

        best_val, best_state, no_improve = np.inf, None, 0
        train_losses, val_losses = [], []

        loader = DataLoader(
            tr_ds, batch_size=batch_size, shuffle=True, drop_last=drop_last
        )
        n_seen = len(tr_ds) - (1 if drop_last else 0)
        for ep in range(1, epochs + 1):
            self._network.train()
            ep_loss = 0.0
            for xb, yb in loader:
                xb, yb = xb.to(dev), yb.to(dev)
                pred = self._resize(self._network(xb), yb.shape[-2:])
                loss = loss_fn(pred, yb)
                opt.zero_grad()
                loss.backward()
                if grad_clip:
                    nn.utils.clip_grad_norm_(
                        self._network.parameters(), grad_clip
                    )
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= n_seen

            self._network.eval()
            with torch.no_grad():
                v_pred = self._resize(self._network(Xva_t), target_hw)
                v_loss = loss_fn(v_pred, yva_t).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val - 1e-6:
                best_val = v_loss
                best_state = copy.deepcopy(self._network.state_dict())
                no_improve = 0
            else:
                no_improve += 1

            if verbose and (ep % max(1, epochs // 10) == 0 or ep == 1):
                print(
                    f"  EMInverter2D  ep {ep:>4d}/{epochs}  "
                    f"train={ep_loss:.5f}  val={v_loss:.5f}"
                )

            if no_improve >= patience:
                if verbose:
                    print(f"  Early stop at epoch {ep} (patience={patience})")
                break

        if best_state is not None:
            self._network.load_state_dict(best_state)

        return {"train_loss": train_losses, "val_loss": val_losses}, best_val

    def _fit_tensorflow(
        self, Xtr, ytr, Xva, yva, *, epochs, batch_size, lr, patience, verbose
    ):
        import tensorflow as tf

        # TF UNet2D expects channels-last: (n, n_freqs, n_sta, n_comp)
        Xtr_tf = Xtr.transpose(0, 2, 3, 1)  # (n, n_freqs, n_sta, n_comp)
        Xva_tf = Xva.transpose(0, 2, 3, 1)
        # Resize targets to match UNet output (n_freqs, n_sta) via tf.image.resize
        target_size = [Xtr.shape[2], Xtr.shape[3]]  # [n_freqs, n_sta]
        ytr_up = tf.image.resize(
            ytr[:, :, :, np.newaxis], target_size
        ).numpy()  # (n, n_freqs, n_sta, 1)
        yva_up = tf.image.resize(yva[:, :, :, np.newaxis], target_size).numpy()

        dev = resolve_device(self.device)
        with tf.device(dev):
            self._network.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
                loss="mse",
            )
            hist = self._network.fit(
                Xtr_tf,
                ytr_up,
                validation_data=(Xva_tf, yva_up),
                epochs=epochs,
                batch_size=batch_size,
                callbacks=[
                    tf.keras.callbacks.EarlyStopping(
                        monitor="val_loss",
                        patience=patience,
                        restore_best_weights=True,
                        min_delta=1e-6,
                    ),
                    tf.keras.callbacks.ReduceLROnPlateau(
                        monitor="val_loss",
                        factor=0.5,
                        patience=max(5, patience // 3),
                        min_lr=1e-6,
                    ),
                ],
                verbose=1 if verbose else 0,
            )
        best_val = min(hist.history["val_loss"])
        return {
            "train_loss": hist.history["loss"],
            "val_loss": hist.history["val_loss"],
        }, best_val

    def predict(
        self,
        X: np.ndarray,
        *,
        as_log_rho: bool = True,
    ) -> np.ndarray:
        """
        Predict 2-D resistivity sections.

        Parameters
        ----------
        X : ndarray (n_profiles, n_components, n_freqs, n_stations)
        as_log_rho : bool
            If ``True`` (default) output is log₁₀(ρ); otherwise
            linear ρ.

        Returns
        -------
        rho_2d : ndarray (n_profiles, n_depth, n_stations)
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before predict().")

        X = np.asarray(X, dtype=np.float32)
        Xn = (X - self._x_mean) / self._x_std

        target_hw = (self.n_depth, self.n_stations)

        if self._backend_name == "tensorflow":
            import tensorflow as tf

            # channels-last: (n, n_freqs, n_sta, n_comp) → (n, n_depth, n_sta, 1)
            Xn_tf = Xn.transpose(0, 2, 3, 1)
            out_tf = self._network.predict(
                Xn_tf, verbose=0
            )  # (n, ?, n_sta, 1)
            y_norm = out_tf.squeeze(-1)  # (n, ?, n_sta)
            # Resize depth axis if UNet output height != n_depth
            if y_norm.shape[1] != self.n_depth:
                y_t = (
                    tf.image.resize(
                        y_norm[:, :, :, np.newaxis],
                        [self.n_depth, self.n_stations],
                    )
                    .numpy()
                    .squeeze(-1)
                )
                y_norm = y_t
        else:
            import torch

            dev = next(self._network.parameters()).device
            self._network.eval()
            batch_size = 16
            outs = []
            for i in range(0, len(Xn), batch_size):
                xb = torch.from_numpy(Xn[i : i + batch_size]).to(dev)
                with torch.no_grad():
                    raw = self._network(xb)  # (b, 1, n_freqs, n_sta)
                    pred = self._resize(raw, target_hw).squeeze(
                        1
                    )  # (b, n_depth, n_sta)
                outs.append(pred.cpu().numpy())
            y_norm = np.concatenate(outs, axis=0)

        y_log = y_norm * self._y_std + self._y_mean  # (n, n_depth, n_sta)

        if as_log_rho:
            return y_log
        return 10.0**y_log

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        p = {
            "n_components": self.n_components,
            "n_depth": self.n_depth,
            "n_stations": self.n_stations,
            "n_freqs": self.n_freqs,
            "arch": self.arch,
            "channels": list(self._channels),
            "dropout": self.dropout,
            "device": self.device,
            "log_rho_out": self.log_rho_out,
        }
        p.update(self._net_kwargs)
        return p

    def _get_weights(self) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {}
        if self._network is not None:
            out.update(get_weights(self._network))
        for attr in ("_x_mean", "_x_std", "_y_mean", "_y_std"):
            val = getattr(self, attr, None)
            if val is not None:
                out[attr] = np.array([val])
        if self._backend_name:
            out["_backend"] = np.array(self._backend_name)
        out["_channels"] = np.array(list(self._channels))
        return out

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        backend_name = str(weights.pop("_backend", np.array("torch")))
        self._backend_name = backend_name

        from pycsamt.backends import set_backend

        set_backend(backend_name)

        if "_channels" in weights:
            self._channels = tuple(
                int(c) for c in weights.pop("_channels").tolist()
            )

        for attr in ("_x_mean", "_x_std", "_y_mean", "_y_std"):
            if attr in weights:
                setattr(self, attr, float(weights.pop(attr)[0]))
        self._network = self._build_network()
        set_weights(self._network, weights)
        self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"EMInverter2D(arch={self.arch!r}, {status})"

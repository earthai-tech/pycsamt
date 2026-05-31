# -*- coding: utf-8 -*-
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
>>> inv.fit(X_train, y_train, epochs=30)              # doctest: +SKIP
EMInverter2D(arch='unet', fitted)
>>> rho_pred = inv.predict(X_test)                    # doctest: +SKIP
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .._base import BaseEMNet, EMCheckpoint

__all__ = ["EMInverter2D"]

# ─────────────────────────────────────────────────────────────────────────────
# Architecture registry (populated lazily)
# ─────────────────────────────────────────────────────────────────────────────

_ARCH_REGISTRY_2D: Dict[str, Any] = {}


def _lazy_register_2d() -> None:
    if _ARCH_REGISTRY_2D:
        return
    from ..nets.unet import UNet2DNet

    def _unet(n_comp, n_out, **kw):
        return UNet2DNet(n_in=n_comp, n_out=n_out, **kw).build()

    _ARCH_REGISTRY_2D["unet"] = _unet


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

    def __init__(
        self,
        n_components: int = 4,
        n_depth: int = 40,
        n_stations: int = 20,
        n_freqs: int = 32,
        *,
        arch: str = "unet",
        device: Optional[str] = None,
        log_rho_out: bool = True,
        **net_kwargs,
    ) -> None:
        super().__init__(arch=arch, n_layers=n_depth, solver="mt2d", device=device)
        self.n_components = int(n_components)
        self.n_depth = int(n_depth)
        self.n_stations = int(n_stations)
        self.n_freqs = int(n_freqs)
        self.log_rho_out = bool(log_rho_out)
        self._net_kwargs = net_kwargs

        self._x_mean: Optional[float] = None
        self._x_std: Optional[float] = None
        self._y_mean: Optional[float] = None
        self._y_std: Optional[float] = None

    # ─── BaseEMNet interface ──────────────────────────────────────────────

    def _build_network(self) -> Any:
        _lazy_register_2d()
        return _ARCH_REGISTRY_2D[self.arch](
            self.n_components, 1, **self._net_kwargs
        )

    def fit(
        self,
        X: np.ndarray,
        y: Optional[np.ndarray] = None,
        *,
        epochs: int = 100,
        batch_size: int = 16,
        lr: float = 1e-3,
        patience: int = 15,
        val_frac: float = 0.1,
        grad_clip: Optional[float] = 1.0,
        seed: Optional[int] = None,
        verbose: bool = True,
    ) -> "EMInverter2D":
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

        Returns
        -------
        self
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import TensorDataset, DataLoader
        except ImportError as exc:
            raise ImportError("PyTorch is required for EMInverter2D.fit") from exc

        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.float32)

        # Normalise inputs/targets globally (simpler than per-channel for 2D)
        self._x_mean = float(np.nanmean(X))
        self._x_std = float(np.nanstd(X)) + 1e-8
        self._y_mean = float(np.nanmean(y))
        self._y_std = float(np.nanstd(y)) + 1e-8

        Xn = (X - self._x_mean) / self._x_std
        yn = (y - self._y_mean) / self._y_std

        # Add channel dim to y: (n, n_depth, n_sta) → (n, 1, n_depth, n_sta)
        yn = yn[:, np.newaxis, :, :]

        rng = np.random.default_rng(seed)
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        vi, ti = idx[:n_val], idx[n_val:]

        dev = self._resolve_device()
        self._network = self._build_network().to(dev)

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=max(5, patience // 3), min_lr=1e-6
        )
        mse = nn.MSELoss()

        Xtr = torch.from_numpy(Xn[ti])
        ytr = torch.from_numpy(yn[ti])
        tr_ds = TensorDataset(Xtr, ytr)

        Xva = torch.from_numpy(Xn[vi]).to(dev)
        yva = torch.from_numpy(yn[vi]).to(dev)

        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []
        no_improve = 0

        for ep in range(1, epochs + 1):
            self._network.train()
            loader = DataLoader(tr_ds, batch_size=batch_size, shuffle=True)
            ep_loss = 0.0
            for xb, yb in loader:
                xb, yb = xb.to(dev), yb.to(dev)
                pred = self._network(xb)
                loss = mse(pred, yb)
                opt.zero_grad()
                loss.backward()
                if grad_clip:
                    torch.nn.utils.clip_grad_norm_(
                        self._network.parameters(), grad_clip
                    )
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= len(ti)

            self._network.eval()
            with torch.no_grad():
                v_loss = mse(self._network(Xva), yva).item()

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

        self._history = {"train_loss": train_losses, "val_loss": val_losses}
        self._meta["best_val_loss"] = float(best_val)
        self._is_fitted = True
        return self

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
        try:
            import torch
        except ImportError as exc:
            raise ImportError("PyTorch is required for EMInverter2D.predict") from exc

        X = np.asarray(X, dtype=np.float32)
        Xn = (X - self._x_mean) / self._x_std

        dev = next(self._network.parameters()).device
        self._network.eval()

        batch_size = 16
        outs = []
        for i in range(0, len(Xn), batch_size):
            xb = torch.from_numpy(Xn[i:i + batch_size]).to(dev)
            with torch.no_grad():
                pred = self._network(xb).squeeze(1).cpu().numpy()
            outs.append(pred)

        y_norm = np.concatenate(outs, axis=0)
        y_log = y_norm * self._y_std + self._y_mean  # (n, n_depth, n_sta)

        if as_log_rho:
            return y_log
        return 10.0 ** y_log

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> Dict[str, Any]:
        p = {
            "n_components": self.n_components,
            "n_depth": self.n_depth,
            "n_stations": self.n_stations,
            "n_freqs": self.n_freqs,
            "arch": self.arch,
            "device": self.device,
            "log_rho_out": self.log_rho_out,
        }
        p.update(self._net_kwargs)
        return p

    def _get_weights(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        if self._network is not None:
            for k, v in self._network.state_dict().items():
                out[k] = v.cpu().numpy()
        for attr in ("_x_mean", "_x_std", "_y_mean", "_y_std"):
            val = getattr(self, attr, None)
            if val is not None:
                out[attr] = np.array([val])
        return out

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        import torch
        for attr in ("_x_mean", "_x_std", "_y_mean", "_y_std"):
            if attr in weights:
                setattr(self, attr, float(weights.pop(attr)[0]))
        self._network = self._build_network()
        state = {k: torch.from_numpy(v) for k, v in weights.items()}
        self._network.load_state_dict(state)
        self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        return f"EMInverter2D(arch={self.arch!r}, {status})"

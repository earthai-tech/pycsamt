# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PyTorch backend implementation.

All torch imports are deferred to method bodies — ``import pycsamt``
never fails without PyTorch installed.
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np

from ._base import BackendSpec, NeuralBackend
from ._detect import probe_backend

__all__ = ["TorchBackend"]


class TorchBackend(NeuralBackend):
    """PyTorch backend for pycsamt AI/ML modules."""

    name = "torch"

    @classmethod
    def is_available(cls) -> bool:
        return probe_backend("torch")

    # ─── device ───────────────────────────────────────────────────────────

    def resolve_device(self, device: str | None = None) -> str:
        if device is not None:
            return device
        try:
            import torch

            if torch.cuda.is_available():
                return "cuda"
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
        except ImportError:
            pass
        return "cpu"

    # ─── model building ───────────────────────────────────────────────────

    def build(self, spec: BackendSpec) -> Any:
        self.require()
        arch = spec["arch"].lower()
        builders = {
            "cnn1d": self._build_cnn1d,
            "cnn": self._build_cnn1d,
            "resnet1d": self._build_resnet1d,
            "resnet": self._build_resnet1d,
            "fcn1d": self._build_fcn1d,
            "fcn": self._build_fcn1d,
            "unet2d": self._build_unet2d,
            "unet": self._build_unet2d,
            "drcnn": self._build_drcnn,
        }
        if arch not in builders:
            raise ValueError(
                f"Unknown arch {arch!r} for torch backend.  "
                f"Available: {list(builders)}"
            )
        return builders[arch](spec)

    def _build_cnn1d(self, spec):
        from ..ai.nets.cnn1d import CNN1DNet

        return CNN1DNet(
            spec["n_features"],
            spec["n_out"],
            channels=spec.get("channels", (32, 64, 128)),
            kernel_size=spec.get("kernel_size", 5),
            dropout=spec.get("dropout", 0.3),
        ).build()

    def _build_resnet1d(self, spec):
        from ..ai.nets.resnet import ResNet1DNet

        return ResNet1DNet(
            spec["n_features"],
            spec["n_out"],
            channels=spec.get("channels", (64, 128, 256)),
            dropout=spec.get("dropout", 0.3),
            n_blocks=spec.get("n_blocks", 2),
        ).build()

    def _build_fcn1d(self, spec):
        from ..ai.nets.fcn import FCN1DNet

        return FCN1DNet(
            spec["n_features"],
            spec["n_out"],
            channels=spec.get("channels", (32, 64, 128, 64)),
            dropout=spec.get("dropout", 0.2),
        ).build()

    def _build_unet2d(self, spec):
        from ..ai.nets.unet import UNet2DNet

        return UNet2DNet(
            n_in=spec["n_in"],
            n_out=spec.get("n_out", 1),
            channels=spec.get("channels", (32, 64, 128, 256, 512)),
            dropout=spec.get("dropout", 0.2),
        ).build()

    def _build_drcnn(self, spec):
        from ..ai.nets.drcnn import DRCNNNet

        return DRCNNNet(
            n_features_list=spec["n_features_list"],
            n_out=spec["n_out"],
            growth_rate=spec.get("growth_rate", 32),
            n_layers=spec.get("n_layers", 6),
            hidden_dim=spec.get("hidden_dim", 256),
            dropout=spec.get("dropout", 0.2),
        ).build()

    # ─── training ─────────────────────────────────────────────────────────

    def train(
        self,
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        *,
        epochs=100,
        batch_size=256,
        lr=1e-3,
        patience=20,
        grad_clip=None,
        device=None,
        verbose=True,
        loss="mse",
    ):
        self.require()
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        dev = self.resolve_device(device)
        model = model.to(dev)

        loss_fn = self._get_loss_fn(loss)
        opt = torch.optim.Adam(model.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=max(5, patience // 3), min_lr=1e-6
        )

        Xtr = torch.from_numpy(X_train.astype(np.float32))
        ytr = torch.from_numpy(y_train.astype(np.float32))
        Xva = torch.from_numpy(X_val.astype(np.float32)).to(dev)
        yva = torch.from_numpy(y_val.astype(np.float32)).to(dev)

        loader = DataLoader(
            TensorDataset(Xtr, ytr),
            batch_size=batch_size,
            shuffle=True,
            pin_memory=(dev != "cpu"),
        )

        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []
        no_improve = 0

        for ep in range(1, epochs + 1):
            model.train()
            ep_loss = 0.0
            for xb, yb in loader:
                xb, yb = xb.to(dev), yb.to(dev)
                pred = model(xb)
                l = loss_fn(pred, yb)
                opt.zero_grad()
                l.backward()
                if grad_clip:
                    nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
                opt.step()
                ep_loss += l.item() * len(xb)
            ep_loss /= len(X_train)

            model.eval()
            with torch.no_grad():
                v_loss = loss_fn(model(Xva), yva).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val - 1e-6:
                best_val = v_loss
                best_state = copy.deepcopy(model.state_dict())
                no_improve = 0
            else:
                no_improve += 1

            if verbose and (ep % max(1, epochs // 10) == 0 or ep == 1):
                print(
                    f"  [torch]  ep {ep:>4d}/{epochs}  "
                    f"train={ep_loss:.5f}  val={v_loss:.5f}"
                )
            if no_improve >= patience:
                if verbose:
                    print(f"  Early stop at epoch {ep}")
                break

        if best_state is not None:
            model.load_state_dict(best_state)

        return {"train_loss": train_losses, "val_loss": val_losses}

    # ─── inference ────────────────────────────────────────────────────────

    def predict(self, model, X, *, device=None, batch_size=256):
        self.require()
        import torch

        dev = next(model.parameters()).device
        model.eval()
        X = np.asarray(X, dtype=np.float32)
        outs = []
        for i in range(0, len(X), batch_size):
            xb = torch.from_numpy(X[i : i + batch_size]).to(dev)
            with torch.no_grad():
                outs.append(model(xb).cpu().numpy())
        return np.concatenate(outs, axis=0)

    # ─── serialisation ────────────────────────────────────────────────────

    def get_weights(self, model) -> dict[str, np.ndarray]:
        return {k: v.cpu().numpy() for k, v in model.state_dict().items()}

    def set_weights(self, model, weights) -> None:
        import torch

        state = {k: torch.from_numpy(v) for k, v in weights.items()}
        model.load_state_dict(state)

    # ─── internal ─────────────────────────────────────────────────────────

    @staticmethod
    def _get_loss_fn(loss: str):
        import torch.nn as nn

        if loss == "mse":
            return nn.MSELoss()
        if loss == "masked_mse":
            from ..ai.training.metrics import masked_mse_loss

            return masked_mse_loss
        if loss == "cross_entropy":
            return nn.CrossEntropyLoss()
        raise ValueError(f"Unknown loss {loss!r}")

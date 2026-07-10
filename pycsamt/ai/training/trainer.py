# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Training loop for EM neural networks.

:class:`EMTrainer` wraps a PyTorch ``nn.Module`` and handles:

* Mini-batch training with ``DataLoader``
* Masked MSE loss (ignores NaN padding for variable-depth datasets)
* ``ReduceLROnPlateau`` learning-rate scheduling
* Early stopping with configurable patience
* Per-epoch progress display (``tqdm`` if available)
* History dict for :func:`~pycsamt.ai.plot.convergence.plot_convergence`
"""
from __future__ import annotations

import time

import numpy as np

__all__ = ["EMTrainer"]

_TQDM_AVAILABLE = False
try:
    from tqdm.auto import tqdm as _tqdm
    _TQDM_AVAILABLE = True
except ImportError:
    pass


class EMTrainer:
    """
    Training loop manager for EM 1-D inversion networks.

    Parameters
    ----------
    model : nn.Module
        The PyTorch network to train.
    lr : float
        Initial learning rate.  Default 1e-3.
    weight_decay : float
        L2 regularisation.  Default 1e-5.
    patience : int
        Early-stopping patience (epochs without val-loss improvement).
    min_delta : float
        Minimum val-loss improvement to reset patience counter.
    batch_size : int
        Mini-batch size.
    device : str
        Compute device (``'cpu'``, ``'cuda'``, ``'mps'``).
    grad_clip : float or None
        If set, clip gradient norms to this value.
    verbose : bool
        Print per-epoch summary.

    Attributes
    ----------
    history : dict
        Keys ``'train_loss'`` and ``'val_loss'`` — lists of per-epoch
        averages.  Also ``'lr'`` and ``'epoch_time'``.
    best_epoch : int
        Epoch index at which the best validation loss was achieved.
    best_val_loss : float
    """

    def __init__(
        self,
        model,
        *,
        lr: float = 1e-3,
        weight_decay: float = 1e-5,
        patience: int = 20,
        min_delta: float = 1e-5,
        batch_size: int = 256,
        device: str = "cpu",
        grad_clip: float | None = None,
        verbose: bool = True,
    ):
        self.model = model
        self.lr = lr
        self.weight_decay = weight_decay
        self.patience = patience
        self.min_delta = min_delta
        self.batch_size = batch_size
        self.device = device
        self.grad_clip = grad_clip
        self.verbose = verbose

        self.history: dict[str, list] = {
            "train_loss": [], "val_loss": [], "lr": [], "epoch_time": []
        }
        self.best_epoch: int = 0
        self.best_val_loss: float = float("inf")
        self._best_state: dict | None = None

    # ─── main fit ──────────────────────────────────────────────────────────

    def fit(
        self,
        train_ds,
        val_ds,
        epochs: int = 100,
    ) -> EMTrainer:
        """
        Train the network.

        Parameters
        ----------
        train_ds : EMDataset
            Training split.
        val_ds : EMDataset
            Validation split.
        epochs : int
            Maximum number of epochs.

        Returns
        -------
        self
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import DataLoader
        except ImportError:
            raise ImportError("PyTorch required for EMTrainer.fit().")

        from .metrics import masked_mse_loss

        device = torch.device(self.device)
        self.model = self.model.to(device)

        optimiser = torch.optim.Adam(
            self.model.parameters(),
            lr=self.lr,
            weight_decay=self.weight_decay,
        )
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimiser, mode="min", factor=0.5, patience=10
        )

        train_loader = DataLoader(
            train_ds, batch_size=self.batch_size, shuffle=True,
            drop_last=False, pin_memory=(self.device != "cpu"),
        )
        val_loader = DataLoader(
            val_ds, batch_size=self.batch_size * 2, shuffle=False,
        )

        patience_ctr = 0

        for epoch in range(1, epochs + 1):
            t0 = time.perf_counter()
            train_loss = self._train_epoch(
                train_loader, optimiser, device, masked_mse_loss
            )
            val_loss = self._eval_epoch(val_loader, device, masked_mse_loss)
            elapsed = time.perf_counter() - t0

            scheduler.step(val_loss)
            current_lr = optimiser.param_groups[0]["lr"]

            self.history["train_loss"].append(train_loss)
            self.history["val_loss"].append(val_loss)
            self.history["lr"].append(current_lr)
            self.history["epoch_time"].append(elapsed)

            # Early stopping
            if val_loss < self.best_val_loss - self.min_delta:
                self.best_val_loss = val_loss
                self.best_epoch = epoch
                self._best_state = {
                    k: v.cpu().clone()
                    for k, v in self.model.state_dict().items()
                }
                patience_ctr = 0
            else:
                patience_ctr += 1

            if self.verbose and (epoch % max(1, epochs // 20) == 0 or epoch == 1):
                print(
                    f"  Epoch {epoch:4d}/{epochs} | "
                    f"train={train_loss:.5f}  val={val_loss:.5f}  "
                    f"lr={current_lr:.2e}  [{elapsed:.1f}s]"
                )

            if patience_ctr >= self.patience:
                if self.verbose:
                    print(
                        f"  Early stop at epoch {epoch} "
                        f"(best val={self.best_val_loss:.5f} @ epoch {self.best_epoch})"
                    )
                break

        # Restore best weights
        if self._best_state is not None:
            self.model.load_state_dict(self._best_state)
            self.model = self.model.to(device)

        return self

    # ─── internal epoch helpers ────────────────────────────────────────────

    def _train_epoch(self, loader, optimiser, device, loss_fn) -> float:
        import torch
        self.model.train()
        total = 0.0
        n_batch = 0
        for x_batch, y_batch in loader:
            x_batch = x_batch.to(device, non_blocking=True)
            y_batch = y_batch.to(device, non_blocking=True)
            optimiser.zero_grad(set_to_none=True)
            pred = self.model(x_batch)
            loss = loss_fn(pred, y_batch)
            loss.backward()
            if self.grad_clip is not None:
                torch.nn.utils.clip_grad_norm_(
                    self.model.parameters(), self.grad_clip
                )
            optimiser.step()
            total += loss.item()
            n_batch += 1
        return total / max(1, n_batch)

    def _eval_epoch(self, loader, device, loss_fn) -> float:
        import torch
        self.model.eval()
        total = 0.0
        n_batch = 0
        with torch.no_grad():
            for x_batch, y_batch in loader:
                x_batch = x_batch.to(device, non_blocking=True)
                y_batch = y_batch.to(device, non_blocking=True)
                pred = self.model(x_batch)
                total += loss_fn(pred, y_batch).item()
                n_batch += 1
        return total / max(1, n_batch)

    # ─── weight helpers ────────────────────────────────────────────────────

    def get_weights(self) -> dict[str, np.ndarray]:
        """Return model weights as numpy dict."""
        return {
            k: v.cpu().numpy()
            for k, v in self.model.state_dict().items()
        }

    def load_weights(self, weights: dict[str, np.ndarray]) -> None:
        """Restore weights from numpy dict."""
        import torch
        state = {k: torch.from_numpy(v) for k, v in weights.items()}
        self.model.load_state_dict(state)

    def __repr__(self) -> str:
        return (
            f"EMTrainer(lr={self.lr}, batch_size={self.batch_size}, "
            f"patience={self.patience}, device={self.device!r})"
        )

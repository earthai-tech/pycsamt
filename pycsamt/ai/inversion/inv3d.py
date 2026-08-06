# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
GCNInverter3D — graph-convolutional 3-D MT spatial inversion.

Exploits inter-station spatial relationships by representing the survey
network as a graph whose adjacency encodes station proximity.  Each
node receives features from its neighbours via Kipf & Welling (2017)
spectral graph convolutions before predicting a per-station 1-D
subsurface model, yielding a spatially consistent 3-D resistivity image.

Input / output convention
--------------------------
* **Input** ``X``: ``ndarray (n_samples, n_stations, n_features)`` or
  ``(n_stations, n_features)`` for a single survey.
  Typical features: ``[log10(rho_a_f1), phi_f1, …, log10(rho_a_fK), phi_fK]``
  giving ``n_features = 2 × K``.

* **Target** ``y``: ``ndarray (n_samples, n_stations, 2*n_layers-1)`` —
  per-station model parameters (log10(ρ) for each layer concatenated with
  log10(thickness) for each interface).

* **Adjacency** ``A``: ``ndarray (n_stations, n_stations)`` — symmetric
  normalised adjacency :math:`\\tilde{D}^{-1/2}\\tilde{A}\\tilde{D}^{-1/2}`.
  Build from station coordinates with
  :func:`~pycsamt.ai.nets.gcn.build_adjacency`.

Message-passing rule (Kipf & Welling 2017)
------------------------------------------
.. math::

    H^{(l+1)} = \\sigma\\!\\bigl(
        \\tilde{D}^{-1/2}\\tilde{A}\\tilde{D}^{-1/2}\\,
        H^{(l)} W^{(l)}\\bigr)

where :math:`\\tilde{A} = A + I` (self-loops) and
:math:`\\tilde{D}_{ii} = \\sum_j \\tilde{A}_{ij}`.

References
----------
.. [1] Kipf, T. N. & Welling, M. (2017). Semi-supervised classification
       with graph convolutional networks. *ICLR 2017*.

Example
-------
>>> import numpy as np
>>> from pycsamt.ai.nets.gcn import build_adjacency
>>> from pycsamt.ai.inversion.inv3d import GCNInverter3D
>>> rng = np.random.default_rng(0)
>>> coords = rng.uniform(0, 10_000, (30, 2))
>>> A = build_adjacency(coords, radius=3_000)
>>> X = rng.standard_normal((200, 30, 40)).astype("f4")
>>> y = rng.standard_normal((200, 30, 9)).astype("f4")
>>> inv = GCNInverter3D(n_features=40, n_layers=5)
>>> inv.fit(X, y, adjacency=A, epochs=5, verbose=False)  # doctest: +SKIP
GCNInverter3D(n_features=40, n_layers=5, ..., n_stations=30, fitted)
"""

from __future__ import annotations

import copy
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from .._backend_utils import (
    active_backend,
    get_weights,
    resolve_device,
    set_weights,
)
from .._base import BaseEMNet

__all__ = ["GCNInverter3D"]


class GCNInverter3D(BaseEMNet):
    """
    Graph-convolutional 3-D MT inversion estimator.

    Models the spatial context of a multi-station survey as a graph and
    applies spectral GCN message-passing before regressing per-station
    1-D subsurface models.  The result is a spatially coherent 3-D
    resistivity volume without any external graph library dependency.

    Parameters
    ----------
    n_features : int, default 40
        Per-station input feature dimension.  A typical choice is
        ``2 × n_freqs`` (log10(ρ_a) + phase at each frequency).
    n_layers : int, default 5
        Number of depth layers per station.  Output size per station is
        ``2 * n_layers - 1`` (n_layers log10(ρ) + (n_layers-1) log10(h)).
    hidden : sequence of int, default (256, 128, 64)
        Width of each GCN message-passing layer.
    dropout : float, default 0.1
        Dropout probability applied between GCN layers during training.
    device : str or None
        Compute device (``'cpu'``, ``'cuda'``, ``'gpu:0'``, …).
        ``None`` auto-detects.
    **net_kwargs
        Extra keyword arguments forwarded to
        :class:`~pycsamt.ai.nets.gcn.GCNNet`.
    """

    def __init__(
        self,
        n_features: int = 40,
        n_layers: int = 5,
        hidden: Sequence[int] = (256, 128, 64),
        dropout: float = 0.1,
        device: str | None = None,
        **net_kwargs,
    ) -> None:
        super().__init__(
            arch="gcn", n_layers=n_layers, solver="mt3d", device=device
        )
        self.n_features = int(n_features)
        self.n_out = 2 * n_layers - 1
        self.hidden = tuple(int(h) for h in hidden)
        self.dropout = float(dropout)
        self._net_kwargs = net_kwargs

        self._x_mean: np.ndarray | None = None
        self._x_std: np.ndarray | None = None
        self._y_mean: np.ndarray | None = None
        self._y_std: np.ndarray | None = None
        self._backend_name: str | None = None
        self._A_stored: np.ndarray | None = None  # adjacency from training

    # ── internal helpers ──────────────────────────────────────────────────────

    def _build_network(self) -> Any:
        from ..nets.gcn import GCNNet

        factory = GCNNet(
            n_features=self.n_features,
            n_out=self.n_out,
            hidden=self.hidden,
            dropout=self.dropout,
            **self._net_kwargs,
        )
        if self._backend_name == "tensorflow":
            return factory.build_tf()
        return factory.build()

    @staticmethod
    def _prepare_adjacency(
        adjacency: np.ndarray | None,
        coords: np.ndarray | None,
        radius: float,
        n_stations: int,
    ) -> np.ndarray:
        """Return a float32 normalised adjacency ``(n_stations, n_stations)``."""
        from ..nets.gcn import build_adjacency

        if adjacency is not None:
            A = np.asarray(adjacency, dtype=np.float32)
            if A.ndim != 2 or A.shape[0] != A.shape[1]:
                raise ValueError(
                    f"adjacency must be square (n_stations, n_stations); "
                    f"got shape {A.shape}."
                )
            return A
        if coords is not None:
            c = np.asarray(coords, dtype=np.float64)
            if c.ndim != 2 or c.shape[1] != 2:
                raise ValueError(
                    f"coords must be shape (n_stations, 2); got {c.shape}."
                )
            return build_adjacency(c, radius=radius)
        # Identity fallback — no inter-station coupling, with warning
        warnings.warn(
            "Neither adjacency nor coords supplied; using identity adjacency "
            "(no inter-station coupling).  Pass coords= or adjacency= for "
            "full GCN spatial modelling benefits.",
            UserWarning,
            stacklevel=3,
        )
        return np.eye(n_stations, dtype=np.float32)

    # ── fit ───────────────────────────────────────────────────────────────────

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray | None = None,
        adjacency: np.ndarray | None = None,
        *,
        coords: np.ndarray | None = None,
        radius: float = 5_000.0,
        epochs: int = 100,
        batch_size: int = 16,
        lr: float = 1e-3,
        patience: int = 15,
        val_frac: float = 0.15,
        grad_clip: float | None = 1.0,
        seed: int | None = None,
        verbose: bool = True,
    ) -> GCNInverter3D:
        """
        Train the 3-D GCN inversion network.

        Parameters
        ----------
        X : ndarray (n_samples, n_stations, n_features) or
            (n_stations, n_features)
            Per-station MT feature matrices.
        y : ndarray (n_samples, n_stations, 2*n_layers-1) or
            (n_stations, 2*n_layers-1)
            Target per-station model parameters (log10 scale recommended).
        adjacency : ndarray (n_stations, n_stations), optional
            Pre-computed normalised adjacency matrix.  If ``None``,
            *coords* and *radius* must be provided.
        coords : ndarray (n_stations, 2), optional
            Station (x, y) positions used to build the adjacency when
            *adjacency* is not given.
        radius : float
            Maximum inter-station edge distance in the same units as
            *coords*; ignored when *adjacency* is supplied.
        epochs, batch_size, lr, patience, val_frac, grad_clip, seed, verbose
            Standard training hyper-parameters.

        Returns
        -------
        self
        """
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.float32)

        if X.ndim == 2:  # single survey → add sample dim
            X, y = X[np.newaxis], y[np.newaxis]
        if X.ndim != 3:
            raise ValueError(f"X must be 2-D or 3-D; got {X.ndim}-D.")

        _, n_stations, _ = X.shape

        A = self._prepare_adjacency(adjacency, coords, radius, n_stations)
        self._A_stored = A.copy()

        # Per-feature/per-output z-score normalisation (over the sample and
        # station/triangle axes, keeping the last axis separate). A single
        # global scalar here would let whichever feature happens to have
        # the largest raw magnitude (e.g. phase in degrees, O(10-100))
        # dominate the shared std and flatten every other feature -- e.g.
        # log10(rho_a) at O(1) or a normalized [0, 1] position feature --
        # to a nearly constant value after normalisation, effectively
        # hiding it from the network regardless of its real information
        # content. See the mt2d_tri position-feature addition in
        # pycsamt.agents.inv2d_agent for the case this was found from.
        # np.nanmean/nanstd upcast float32 input to float64 for a vector
        # result (unlike a Python-float scalar, which numpy's weak-scalar
        # promotion rules leave float32 arithmetic alone) -- cast back
        # explicitly so Xn/yn below, and every later predict() call, stay
        # float32 and match the float32 network weights.
        self._x_mean = np.nanmean(X, axis=(0, 1)).astype(np.float32)
        self._x_std = (np.nanstd(X, axis=(0, 1)) + 1e-8).astype(np.float32)
        self._y_mean = np.nanmean(y, axis=(0, 1)).astype(np.float32)
        self._y_std = (np.nanstd(y, axis=(0, 1)) + 1e-8).astype(np.float32)

        Xn = (X - self._x_mean) / self._x_std
        yn = (y - self._y_mean) / self._y_std

        # Train / validation split
        rng = np.random.default_rng(seed)
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        vi, ti = idx[:n_val], idx[n_val:]

        self._backend_name = active_backend()
        if self._backend_name == "none":
            raise ImportError(
                "No DL backend found.  Install PyTorch or TensorFlow:\n"
                "    pip install torch\n"
                "    pip install tensorflow"
            )

        self._network = self._build_network()

        if self._backend_name == "tensorflow":
            hist, best_val = self._fit_tensorflow(
                Xn[ti],
                yn[ti],
                Xn[vi],
                yn[vi],
                A,
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
                A,
                epochs=epochs,
                batch_size=batch_size,
                lr=lr,
                patience=patience,
                grad_clip=grad_clip,
                verbose=verbose,
            )

        self._history = hist
        self._meta["best_val_loss"] = float(best_val)
        self._is_fitted = True
        return self

    # ── internal training loops ───────────────────────────────────────────────

    def _fit_torch(
        self,
        Xtr,
        ytr,
        Xva,
        yva,
        A,
        *,
        epochs,
        batch_size,
        lr,
        patience,
        grad_clip,
        verbose,
    ):
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        from pycsamt.api.view.progress import get_progress_bar

        dev = resolve_device(self.device)
        self._network = self._network.to(dev)
        A_t = torch.from_numpy(A).to(dev)  # (n_sta, n_sta) — fixed

        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt,
            factor=0.5,
            patience=max(5, patience // 3),
            min_lr=1e-6,
        )
        mse = nn.MSELoss()

        tr_ds = TensorDataset(torch.from_numpy(Xtr), torch.from_numpy(ytr))
        Xva_t = torch.from_numpy(Xva).to(dev)
        yva_t = torch.from_numpy(yva).to(dev)

        best_val, best_state, no_improve = np.inf, None, 0
        train_losses, val_losses = [], []

        with get_progress_bar(
            total=epochs, desc="GCNInverter3D", unit="epoch", verbose=verbose
        ) as bar:
            for ep in range(1, epochs + 1):
                self._network.train()
                ep_loss = 0.0

                for xb, yb in DataLoader(
                    tr_ds, batch_size=batch_size, shuffle=True
                ):
                    xb, yb = xb.to(dev), yb.to(dev)  # (b, n_sta, n_feat/n_out)
                    b = xb.shape[0]
                    # Process each survey in the mini-batch; n_sta is
                    # typically small (< 300) so per-survey iteration has
                    # negligible overhead
                    preds = torch.stack(
                        [self._network(xb[i], A_t) for i in range(b)], dim=0
                    )  # (b, n_sta, n_out)
                    loss = mse(preds, yb)
                    opt.zero_grad()
                    loss.backward()
                    if grad_clip:
                        nn.utils.clip_grad_norm_(
                            self._network.parameters(), grad_clip
                        )
                    opt.step()
                    ep_loss += loss.item() * b

                ep_loss /= len(Xtr)

                self._network.eval()
                with torch.no_grad():
                    val_preds = torch.stack(
                        [
                            self._network(Xva_t[i], A_t)
                            for i in range(len(Xva_t))
                        ],
                        dim=0,
                    )
                    v_loss = mse(val_preds, yva_t).item()

                sched.step(v_loss)
                train_losses.append(ep_loss)
                val_losses.append(v_loss)

                if v_loss < best_val - 1e-6:
                    best_val = v_loss
                    best_state = copy.deepcopy(self._network.state_dict())
                    no_improve = 0
                else:
                    no_improve += 1

                bar.update(1, metrics={"train": ep_loss, "val": v_loss})

                if no_improve >= patience:
                    if bar.verbose:
                        print(
                            f"  Early stop at epoch {ep} (patience={patience})"
                        )
                    break

        if best_state is not None:
            self._network.load_state_dict(best_state)

        return {"train_loss": train_losses, "val_loss": val_losses}, best_val

    def _fit_tensorflow(
        self,
        Xtr,
        ytr,
        Xva,
        yva,
        A,
        *,
        epochs,
        batch_size,
        lr,
        patience,
        verbose,
    ):
        import tensorflow as tf

        from pycsamt.api.view.progress import get_progress_bar

        A_t = tf.constant(A)  # (n_sta, n_sta)

        def _forward_batch(X_batch, A_c, training):
            b = X_batch.shape[0]
            return tf.stack(
                [
                    self._network([X_batch[i], A_c], training=training)
                    for i in range(b)
                ],
                axis=0,
            )

        def make_dataset(X, y, shuffle=False):
            ds = tf.data.Dataset.from_tensor_slices(
                (X.astype(np.float32), y.astype(np.float32))
            )
            if shuffle:
                ds = ds.shuffle(len(X), reshuffle_each_iteration=True)
            return ds.batch(batch_size)

        tr_ds = make_dataset(Xtr, ytr, shuffle=True)
        va_ds = make_dataset(Xva, yva)

        opt = tf.keras.optimizers.Adam(learning_rate=lr)
        mse = tf.keras.losses.MeanSquaredError()
        best_val, best_weights, no_improve = np.inf, None, 0
        train_losses, val_losses = [], []

        with get_progress_bar(
            total=epochs,
            desc="GCNInverter3D (TF)",
            unit="epoch",
            verbose=verbose,
        ) as bar:
            for ep in range(1, epochs + 1):
                ep_loss, n_batches = 0.0, 0

                for xb, yb in tr_ds:
                    with tf.GradientTape() as tape:
                        preds = _forward_batch(xb, A_t, training=True)
                        loss = mse(yb, preds)
                    grads = tape.gradient(
                        loss, self._network.trainable_variables
                    )
                    opt.apply_gradients(
                        zip(grads, self._network.trainable_variables)
                    )
                    ep_loss += float(loss)
                    n_batches += 1
                ep_loss /= max(n_batches, 1)

                v_loss, n_vb = 0.0, 0
                for xb, yb in va_ds:
                    preds = _forward_batch(xb, A_t, training=False)
                    v_loss += float(mse(yb, preds))
                    n_vb += 1
                v_loss /= max(n_vb, 1)

                train_losses.append(ep_loss)
                val_losses.append(v_loss)

                if v_loss < best_val - 1e-6:
                    best_val = v_loss
                    best_weights = self._network.get_weights()
                    no_improve = 0
                else:
                    no_improve += 1

                bar.update(1, metrics={"train": ep_loss, "val": v_loss})

                if no_improve >= patience:
                    if bar.verbose:
                        print(
                            f"  Early stop at epoch {ep} (patience={patience})"
                        )
                    break

        if best_weights is not None:
            self._network.set_weights(best_weights)

        return {"train_loss": train_losses, "val_loss": val_losses}, best_val

    # ── predict ───────────────────────────────────────────────────────────────

    def predict(
        self,
        X: np.ndarray,
        adjacency: np.ndarray | None = None,
        *,
        coords: np.ndarray | None = None,
        radius: float = 5_000.0,
        as_log_rho: bool = True,
    ) -> np.ndarray:
        """
        Predict per-station subsurface models.

        Parameters
        ----------
        X : ndarray (n_stations, n_features) or
            (n_samples, n_stations, n_features)
        adjacency : ndarray (n_stations, n_stations), optional
            If ``None``, uses the adjacency stored from :meth:`fit`.
        coords : ndarray (n_stations, 2), optional
            Build adjacency from coordinates when *adjacency* is absent.
        radius : float
            Edge radius (used only when *coords* is supplied).
        as_log_rho : bool
            Return log10(ρ) when ``True``; linear-scale ρ otherwise.

        Returns
        -------
        y_pred : ndarray (n_stations, n_out) or (n_samples, n_stations, n_out)
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before predict().")

        X = np.asarray(X, dtype=np.float32)
        squeeze = X.ndim == 2
        if squeeze:
            X = X[np.newaxis]

        n_samples, n_stations, _ = X.shape
        Xn = (X - self._x_mean) / self._x_std

        if adjacency is not None or coords is not None:
            A = self._prepare_adjacency(adjacency, coords, radius, n_stations)
        elif self._A_stored is not None:
            A = self._A_stored
        else:
            A = np.eye(n_stations, dtype=np.float32)

        if self._backend_name == "tensorflow":
            import tensorflow as tf

            A_t = tf.constant(A)
            out = np.stack(
                [
                    self._network([Xn[i], A_t], training=False).numpy()
                    for i in range(n_samples)
                ],
                axis=0,
            )
        else:
            import torch

            dev = next(self._network.parameters()).device
            A_t = torch.from_numpy(A).to(dev)
            self._network.eval()
            with torch.no_grad():
                out = np.stack(
                    [
                        self._network(torch.from_numpy(Xn[i]).to(dev), A_t)
                        .cpu()
                        .numpy()
                        for i in range(n_samples)
                    ],
                    axis=0,
                )

        y_log = out * self._y_std + self._y_mean  # de-normalise

        if squeeze:
            y_log = y_log[0]

        return y_log if as_log_rho else 10.0**y_log

    def predict_with_uncertainty(
        self,
        X: np.ndarray,
        adjacency: np.ndarray | None = None,
        *,
        coords: np.ndarray | None = None,
        radius: float = 5_000.0,
        n_mc: int = 30,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        MC-dropout uncertainty estimate for 3-D predictions.

        Runs *n_mc* stochastic forward passes with dropout active and
        returns the mean and pointwise standard deviation.

        Parameters
        ----------
        X : ndarray (..., n_stations, n_features)
        adjacency : ndarray, optional
        coords : ndarray, optional
        radius : float
        n_mc : int
            Number of Monte-Carlo dropout samples.

        Returns
        -------
        mean : ndarray — same shape as :meth:`predict` output
        std  : ndarray — same shape
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before predict_with_uncertainty().")

        X = np.asarray(X, dtype=np.float32)
        squeeze = X.ndim == 2
        if squeeze:
            X = X[np.newaxis]
        n_samples, n_stations, _ = X.shape
        Xn = (X - self._x_mean) / self._x_std

        if adjacency is not None or coords is not None:
            A = self._prepare_adjacency(adjacency, coords, radius, n_stations)
        elif self._A_stored is not None:
            A = self._A_stored
        else:
            A = np.eye(n_stations, dtype=np.float32)

        mc_samples = []

        if self._backend_name == "tensorflow":
            import tensorflow as tf

            A_t = tf.constant(A)
            for _ in range(n_mc):
                out = np.stack(
                    [
                        self._network([Xn[i], A_t], training=True).numpy()
                        for i in range(n_samples)
                    ],
                    axis=0,
                )
                mc_samples.append(out * self._y_std + self._y_mean)
        else:
            import torch

            dev = next(self._network.parameters()).device
            A_t = torch.from_numpy(A).to(dev)
            self._network.train()  # activate dropout for MC passes
            with torch.no_grad():
                for _ in range(n_mc):
                    out = np.stack(
                        [
                            self._network(torch.from_numpy(Xn[i]).to(dev), A_t)
                            .cpu()
                            .numpy()
                            for i in range(n_samples)
                        ],
                        axis=0,
                    )
                    mc_samples.append(out * self._y_std + self._y_mean)
            self._network.eval()

        stacked = np.stack(
            mc_samples, axis=0
        )  # (n_mc, n_samples, n_sta, n_out)
        mean = stacked.mean(axis=0)
        std = stacked.std(axis=0)

        if squeeze:
            mean, std = mean[0], std[0]

        return mean, std

    # ── serialisation ─────────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "n_features": self.n_features,
            "n_layers": self.n_layers,
            "hidden": list(self.hidden),
            "dropout": self.dropout,
            "device": self.device,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {}
        if self._network is not None:
            out.update(get_weights(self._network))
        for attr in ("_x_mean", "_x_std", "_y_mean", "_y_std"):
            val = getattr(self, attr, None)
            if val is not None:
                out[attr] = np.atleast_1d(np.asarray(val, dtype=np.float32))
        if self._backend_name:
            out["_backend"] = np.array(self._backend_name)
        if self._A_stored is not None:
            out["_A_stored"] = self._A_stored
        return out

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        backend_name = str(weights.pop("_backend", np.array("torch")))
        self._backend_name = backend_name

        from pycsamt.backends import set_backend

        set_backend(backend_name)

        for attr in ("_x_mean", "_x_std", "_y_mean", "_y_std"):
            if attr in weights:
                setattr(self, attr, np.asarray(weights.pop(attr), dtype=np.float32))
        if "_A_stored" in weights:
            self._A_stored = weights.pop("_A_stored")
        self._network = self._build_network()
        set_weights(self._network, weights)
        self._is_fitted = True

    def __repr__(self) -> str:
        status = "fitted" if self._is_fitted else "unfitted"
        n_sta = self._A_stored.shape[0] if self._A_stored is not None else "?"
        return (
            f"GCNInverter3D(n_features={self.n_features}, "
            f"n_layers={self.n_layers}, hidden={self.hidden}, "
            f"n_stations={n_sta}, {status})"
        )

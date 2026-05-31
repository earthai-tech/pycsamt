# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
JointInverter — multi-modal / multi-physics joint inversion.

Uses the :class:`~pycsamt.ai.nets.drcnn.DRCNNNet` architecture to fuse
two or more EM datasets (or EM + non-EM methods) into a single
subsurface model prediction.

Supported combinations
----------------------
* MT + TEM — complementary depth sensitivity
* MT + CSAMT — frequency overlap to constrain near-surface
* MT + gravity or seismic — cross-property constraints

Joint training
--------------
All modalities must be observed at the same sites.  Use
:func:`~pycsamt.forward.batch.generate_dataset` with different
``solver`` values and the same :class:`~pycsamt.forward.synthetic.LayeredModel`
to generate correlated synthetic datasets for supervised training.

Example
-------
>>> from pycsamt.ai.inversion import JointInverter
>>> inv = JointInverter(n_features_list=(120, 48), n_layers=5)
>>> inv.fit([X_mt, X_tem], y, epochs=50)              # doctest: +SKIP
JointInverter(modalities=2, fitted)
>>> y_pred = inv.predict([X_mt_test, X_tem_test])     # doctest: +SKIP
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .._base import BaseEMNet, EMCheckpoint
from ..training.dataset import Normalizer

__all__ = ["JointInverter"]


class JointInverter(BaseEMNet):
    """
    Multi-modal joint inversion estimator based on DRCNN.

    Parameters
    ----------
    n_features_list : sequence of int
        Feature vector lengths for each modality.
        E.g. ``(120, 48)`` for two modalities.
    n_layers : int, default 5
        Number of earth layers to invert for.
    growth_rate : int, default 32
        Dense block growth rate.
    n_dense_layers : int, default 6
        Sub-layers per dense block.
    hidden_dim : int, default 256
        Encoded-feature dimension for each modality and fusion stage.
    dropout : float, default 0.2
    device : str or None
    log_thickness : bool, default True
        Apply log₁₀ transform to thickness targets before normalisation.
    **net_kwargs
        Forwarded to :class:`~pycsamt.ai.nets.drcnn.DRCNNNet`.
    """

    def __init__(
        self,
        n_features_list: Sequence[int],
        n_layers: int = 5,
        *,
        growth_rate: int = 32,
        n_dense_layers: int = 6,
        hidden_dim: int = 256,
        dropout: float = 0.2,
        device: Optional[str] = None,
        log_thickness: bool = True,
        **net_kwargs,
    ) -> None:
        super().__init__(arch="drcnn", n_layers=n_layers, solver="joint",
                         device=device)
        self.n_features_list = tuple(int(n) for n in n_features_list)
        self.growth_rate = int(growth_rate)
        self.n_dense_layers = int(n_dense_layers)
        self.hidden_dim = int(hidden_dim)
        self.dropout = float(dropout)
        self.log_thickness = bool(log_thickness)
        self._net_kwargs = net_kwargs

        self._n_out = 2 * n_layers - 1  # n_layers rho + (n_layers-1) thick
        self._x_norms: List[Optional[Normalizer]] = [None] * len(n_features_list)
        self._y_norm: Optional[Normalizer] = None

    # ─── BaseEMNet interface ──────────────────────────────────────────────

    def _build_network(self) -> Any:
        from ..nets.drcnn import DRCNNNet
        return DRCNNNet(
            n_features_list=self.n_features_list,
            n_out=self._n_out,
            growth_rate=self.growth_rate,
            n_layers=self.n_dense_layers,
            hidden_dim=self.hidden_dim,
            dropout=self.dropout,
            **self._net_kwargs,
        ).build()

    def fit(
        self,
        X_list: Sequence[np.ndarray],
        y: np.ndarray,
        *,
        epochs: int = 100,
        batch_size: int = 256,
        lr: float = 1e-3,
        patience: int = 20,
        val_frac: float = 0.1,
        grad_clip: Optional[float] = 1.0,
        seed: Optional[int] = None,
        verbose: bool = True,
    ) -> "JointInverter":
        """
        Train the joint inverter.

        Parameters
        ----------
        X_list : list of ndarray
            Feature matrices for each modality, each
            ``(n_samples, n_features_i)``.  Lengths must match.
        y : ndarray (n_samples, 2*n_layers-1)
            Target model parameters: first ``n_layers`` columns are
            log₁₀(ρ) or ρ; last ``n_layers-1`` columns are thicknesses.
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
            raise ImportError("PyTorch is required for JointInverter.fit") from exc

        X_arrs = [np.asarray(Xm, dtype=np.float32) for Xm in X_list]
        y_arr = np.asarray(y, dtype=np.float32)

        # Normalise each modality independently
        self._x_norms = []
        X_normed = []
        for Xm in X_arrs:
            norm = Normalizer().fit(Xm)
            self._x_norms.append(norm)
            X_normed.append(norm.transform(Xm).astype(np.float32))

        # Normalise targets (apply log10 to thickness columns)
        y_proc = y_arr.copy()
        if self.log_thickness and y_proc.shape[1] > self.n_layers:
            y_proc[:, self.n_layers:] = np.log10(
                np.maximum(y_proc[:, self.n_layers:], 1e-3)
            )
        self._y_norm = Normalizer().fit(y_proc)
        yn = self._y_norm.transform(y_proc).astype(np.float32)

        rng = np.random.default_rng(seed)
        n = len(yn)
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

        tr_tensors = [torch.from_numpy(Xm[ti]) for Xm in X_normed]
        tr_tensors.append(torch.from_numpy(yn[ti]))
        tr_ds = TensorDataset(*tr_tensors)

        va_inputs = [torch.from_numpy(Xm[vi]).to(dev) for Xm in X_normed]
        yva = torch.from_numpy(yn[vi]).to(dev)

        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []
        no_improve = 0

        for ep in range(1, epochs + 1):
            self._network.train()
            loader = DataLoader(tr_ds, batch_size=batch_size, shuffle=True)
            ep_loss = 0.0
            for *xbs, yb in loader:
                xbs = [t.to(dev) for t in xbs]
                yb = yb.to(dev)
                pred = self._network(*xbs)
                loss = mse(pred, yb)
                opt.zero_grad()
                loss.backward()
                if grad_clip:
                    torch.nn.utils.clip_grad_norm_(
                        self._network.parameters(), grad_clip
                    )
                opt.step()
                ep_loss += loss.item() * len(yb)
            ep_loss /= len(ti)

            self._network.eval()
            with torch.no_grad():
                v_pred = self._network(*va_inputs)
                v_loss = mse(v_pred, yva).item()

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
                    f"  JointInverter  ep {ep:>4d}/{epochs}  "
                    f"train={ep_loss:.5f}  val={v_loss:.5f}"
                )
            if no_improve >= patience:
                if verbose:
                    print(f"  Early stop at epoch {ep}")
                break

        if best_state is not None:
            self._network.load_state_dict(best_state)

        self._history = {"train_loss": train_losses, "val_loss": val_losses}
        self._meta["best_val_loss"] = float(best_val)
        self._is_fitted = True
        return self

    def predict(
        self,
        X_list: Union[Sequence[np.ndarray], np.ndarray],
        *,
        as_log_rho: bool = True,
    ) -> np.ndarray:
        """
        Predict subsurface model parameters.

        Parameters
        ----------
        X_list : list of ndarray, each (n_samples, n_features_i)
        as_log_rho : bool
            If ``True`` (default), resistivity columns are returned in
            log₁₀(ρ); otherwise linear.

        Returns
        -------
        y_pred : ndarray (n_samples, 2*n_layers-1)
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before predict().")
        try:
            import torch
        except ImportError as exc:
            raise ImportError("PyTorch is required for JointInverter.predict") from exc

        if isinstance(X_list, np.ndarray):
            X_list = [X_list]

        X_normed = [
            norm.transform(np.asarray(Xm, dtype=np.float32))
            for norm, Xm in zip(self._x_norms, X_list)
        ]

        dev = next(self._network.parameters()).device
        self._network.eval()

        batch_size = 256
        n = len(X_normed[0])
        outs = []
        for i in range(0, n, batch_size):
            inputs = [
                torch.from_numpy(Xm[i:i + batch_size]).to(dev)
                for Xm in X_normed
            ]
            with torch.no_grad():
                pred = self._network(*inputs).cpu().numpy()
            outs.append(pred)

        y_norm = np.concatenate(outs, axis=0)
        y_pred = self._y_norm.inverse_transform(y_norm)

        # Undo log10 on thickness columns
        if self.log_thickness and y_pred.shape[1] > self.n_layers:
            y_pred[:, self.n_layers:] = 10.0 ** y_pred[:, self.n_layers:]

        if not as_log_rho:
            y_pred[:, : self.n_layers] = 10.0 ** y_pred[:, : self.n_layers]

        return y_pred

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> Dict[str, Any]:
        return {
            "n_features_list": list(self.n_features_list),
            "n_layers": self.n_layers,
            "growth_rate": self.growth_rate,
            "n_dense_layers": self.n_dense_layers,
            "hidden_dim": self.hidden_dim,
            "dropout": self.dropout,
            "device": self.device,
            "log_thickness": self.log_thickness,
        }

    def _get_weights(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        if self._network is not None:
            for k, v in self._network.state_dict().items():
                out[k] = v.cpu().numpy()
        for i, norm in enumerate(self._x_norms):
            if norm is not None:
                d = norm.to_dict()
                out[f"_xnorm_{i}_mean"] = np.array(d["mean"])
                out[f"_xnorm_{i}_std"] = np.array(d["std"])
        if self._y_norm is not None:
            d = self._y_norm.to_dict()
            out["_ynorm_mean"] = np.array(d["mean"])
            out["_ynorm_std"] = np.array(d["std"])
        return out

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        import torch
        self._x_norms = []
        for i in range(len(self.n_features_list)):
            km, ks = f"_xnorm_{i}_mean", f"_xnorm_{i}_std"
            if km in weights:
                norm = Normalizer()
                norm.mean = weights.pop(km)
                norm.std = weights.pop(ks)
                self._x_norms.append(norm)
            else:
                self._x_norms.append(None)
        if "_ynorm_mean" in weights:
            self._y_norm = Normalizer()
            self._y_norm.mean = weights.pop("_ynorm_mean")
            self._y_norm.std = weights.pop("_ynorm_std")
        self._network = self._build_network()
        state = {k: torch.from_numpy(v) for k, v in weights.items()}
        self._network.load_state_dict(state)
        self._is_fitted = True

    def __repr__(self) -> str:
        n = len(self.n_features_list)
        status = "fitted" if self._is_fitted else "unfitted"
        return f"JointInverter(modalities={n}, {status})"

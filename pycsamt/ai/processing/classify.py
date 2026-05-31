# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
DimensionalityClassifier — MLP-based MT data dimensionality classifier.

Replaces the heuristic skew/ellipticity thresholds in
:func:`~pycsamt.emtools.dimensionality.classify_dimensionality` with a
trained multi-layer perceptron that operates on per-(site, frequency)
phase-tensor-derived features.

Classes
-------
* 0 — 1-D  (low skew, low ellipticity)
* 1 — 2-D  (low skew, high ellipticity)
* 2 — 3-D  (high skew / arbitrary ellipticity)

The network is also equipped with a regression head that predicts the
geoelectric strike direction :math:`\\alpha` (in degrees,
:math:`-90^\\circ` to :math:`90^\\circ`) for 2-D observations.

Feature vector (per site × frequency)
--------------------------------------
``[β_abs, ellipt_abs, Φmin, Φmax, log_{10}ρ_{det}, tip_amp]``

where:

* :math:`|\\beta|` — phase tensor skew (Caldwell 2004)
* ``ellipt_abs`` — phase tensor ellipticity :math:`(\\Phi_{max}-\\Phi_{min})/(\\Phi_{max}+\\Phi_{min})`
* :math:`\\Phi_{min}`, :math:`\\Phi_{max}` — eigenvalues of the phase tensor
* :math:`\\log_{10}\\rho_{det}` — determinant apparent resistivity
* ``tip_amp`` — tipper amplitude :math:`\\sqrt{|T_x|^2 + |T_y|^2}`

Integration with emtools
------------------------
Use :meth:`from_features_table` to construct a ready-to-classify
instance from the DataFrame returned by
:func:`~pycsamt.emtools.dimensionality.phase_features_table`.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .._base import BaseEMProcessor

__all__ = ["DimensionalityClassifier"]

_FEATURE_COLS = ["beta_abs", "ellipt_abs", "phi_min", "phi_max",
                 "logrho_det", "tip_amp"]
_N_FEATURES = len(_FEATURE_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# Network builders
# ─────────────────────────────────────────────────────────────────────────────

def _build_dim_mlp(
    n_features: int,
    n_classes: int,
    hidden: Tuple[int, ...],
    dropout: float,
) -> Any:
    """MLP classification + strike regression head (shared backbone)."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError("PyTorch is required for DimensionalityClassifier") from exc

    dims = [n_features] + list(hidden)

    class _DimMLP(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            backbone_layers: list = []
            for i in range(len(dims) - 1):
                backbone_layers += [
                    nn.Linear(dims[i], dims[i + 1]),
                    nn.BatchNorm1d(dims[i + 1]),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                ]
            self.backbone = nn.Sequential(*backbone_layers)
            self.cls_head = nn.Linear(dims[-1], n_classes)
            # Strike: single scalar in (-90, 90) degrees
            self.strike_head = nn.Sequential(
                nn.Linear(dims[-1], 32),
                nn.ReLU(),
                nn.Linear(32, 1),
            )

        def forward(self, x):
            feat = self.backbone(x)
            return self.cls_head(feat), self.strike_head(feat).squeeze(-1)

    return _DimMLP()


# ─────────────────────────────────────────────────────────────────────────────
# Label generation from rule-based classifier (for self-training)
# ─────────────────────────────────────────────────────────────────────────────

def _rule_labels(
    beta_abs: np.ndarray,
    ellipt_abs: np.ndarray,
    skew_th: float = 3.0,
    ellipt_th: float = 0.2,
) -> np.ndarray:
    """Assign dimensionality labels using the rule-based approach."""
    labels = np.full(len(beta_abs), 2, dtype=int)
    ok2 = beta_abs <= skew_th
    labels[ok2 & (ellipt_abs <= ellipt_th)] = 0
    labels[ok2 & (ellipt_abs > ellipt_th)] = 1
    return labels


# ─────────────────────────────────────────────────────────────────────────────
# DimensionalityClassifier
# ─────────────────────────────────────────────────────────────────────────────

class DimensionalityClassifier(BaseEMProcessor):
    """
    MLP classifier for MT data dimensionality (1-D / 2-D / 3-D).

    Parameters
    ----------
    hidden : tuple of int, default (128, 64)
        Hidden layer widths of the shared MLP backbone.
    dropout : float, default 0.2
        Dropout probability in each hidden layer.
    n_classes : int, default 3
        Number of dimensionality classes (0=1D, 1=2D, 2=3D).
    lr : float, default 1e-3
        Default learning rate; can be overridden in :meth:`fit`.
    device : str or None

    Notes
    -----
    When PyTorch is not available the classifier uses a random-forest
    fallback (scikit-learn).  Both backends expose the same interface.

    Examples
    --------
    >>> from pycsamt.ai.processing import DimensionalityClassifier
    >>> clf = DimensionalityClassifier()
    >>> clf.fit(X_train, y_train, epochs=30)         # doctest: +SKIP
    DimensionalityClassifier(n_classes=3, torch)
    >>> clf.predict(X_test)                           # doctest: +SKIP
    array([0, 1, 2, ...])
    >>> clf.predict_strike(X_2d)                      # doctest: +SKIP
    array([ 35., -12., ...])
    """

    def __init__(
        self,
        hidden: Tuple[int, ...] = (128, 64),
        dropout: float = 0.2,
        n_classes: int = 3,
        lr: float = 1e-3,
        device: Optional[str] = None,
    ) -> None:
        self.hidden = tuple(hidden)
        self.dropout = float(dropout)
        self.n_classes = int(n_classes)
        self.lr = float(lr)
        self.device = device

        self._network: Any = None
        self._rf: Any = None  # sklearn fallback
        self._use_rf: bool = False
        self._x_mean: Optional[np.ndarray] = None
        self._x_std: Optional[np.ndarray] = None
        self._is_fitted: bool = False
        self._history: Dict[str, list] = {}

    # ─── factory from emtools ─────────────────────────────────────────────

    @classmethod
    def from_features_table(
        cls,
        df: pd.DataFrame,
        *,
        label_col: Optional[str] = "dim",
        strike_col: Optional[str] = None,
        skew_th: float = 3.0,
        ellipt_th: float = 0.2,
        **fit_kwargs,
    ) -> "DimensionalityClassifier":
        """
        Construct and train a classifier from a
        :func:`~pycsamt.emtools.dimensionality.phase_features_table`
        DataFrame.

        If the DataFrame has a ``dim`` column the labels are taken from
        there; otherwise they are generated by the rule-based approach.

        Parameters
        ----------
        df : DataFrame
            Output of ``phase_features_table`` or ``classify_dimensionality``.
        label_col : str or None
            Column name for pre-computed labels.  ``None`` → use rules.
        strike_col : str or None
            Column name for strike direction in degrees.
        skew_th : float
            Rule-based skew threshold (used when ``label_col`` is absent).
        ellipt_th : float
            Rule-based ellipticity threshold.
        **fit_kwargs
            Passed to :meth:`fit` (e.g. ``epochs=50``).

        Returns
        -------
        DimensionalityClassifier
        """
        X, y, strike = _df_to_Xy(df, label_col, strike_col, skew_th, ellipt_th)
        obj = cls()
        obj.fit(X, y, strike=strike, **fit_kwargs)
        return obj

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(
        self,
        X: Union[np.ndarray, pd.DataFrame],
        y: Optional[np.ndarray] = None,
        *,
        strike: Optional[np.ndarray] = None,
        epochs: int = 80,
        batch_size: int = 256,
        lr: Optional[float] = None,
        val_frac: float = 0.15,
        seed: Optional[int] = None,
        verbose: bool = True,
    ) -> "DimensionalityClassifier":
        """
        Train the dimensionality classifier.

        Parameters
        ----------
        X : ndarray (n_samples, 6) or DataFrame
            Feature matrix.  If a DataFrame, columns matching
            ``_FEATURE_COLS`` are used automatically.
        y : int ndarray (n_samples,) or None
            Dimensionality labels (0/1/2).  Required when ``X`` is an
            array; inferred from rules when ``X`` is a site-feature
            DataFrame.
        strike : float ndarray (n_samples,) or None
            Ground-truth strike angles in degrees.  ``NaN`` for non-2D
            observations.
        epochs : int
        batch_size : int
        lr : float or None
            Overrides the constructor ``lr`` when given.
        val_frac : float
        seed : int or None
        verbose : bool

        Returns
        -------
        self
        """
        X_arr, y_arr, strike_arr = self._coerce_Xy(X, y, strike)

        self._x_mean = X_arr.mean(axis=0, keepdims=True)
        self._x_std = X_arr.std(axis=0, keepdims=True) + 1e-8
        Xn = (X_arr - self._x_mean) / self._x_std

        _lr = float(lr) if lr is not None else self.lr

        try:
            self._fit_torch(Xn, y_arr, strike_arr,
                            epochs=epochs, batch_size=batch_size,
                            lr=_lr, val_frac=val_frac, seed=seed,
                            verbose=verbose)
            self._use_rf = False
        except ImportError:
            self._fit_rf(Xn, y_arr, verbose=verbose)
            self._use_rf = True

        self._is_fitted = True
        return self

    def transform(self, X: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
        """
        Compute class probabilities.

        Parameters
        ----------
        X : ndarray (n_samples, 6) or DataFrame

        Returns
        -------
        proba : ndarray, shape (n_samples, n_classes)
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")
        X_arr = self._coerce_X(X)
        Xn = (X_arr - self._x_mean) / self._x_std
        return self._predict_proba(Xn)

    def predict(self, X: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
        """
        Predict dimensionality class (0=1D, 1=2D, 2=3D).

        Parameters
        ----------
        X : ndarray (n_samples, 6) or DataFrame

        Returns
        -------
        labels : int ndarray, shape (n_samples,)
        """
        proba = self.transform(X)
        return proba.argmax(axis=1)

    def predict_strike(
        self, X: Union[np.ndarray, pd.DataFrame]
    ) -> np.ndarray:
        """
        Predict geoelectric strike direction.

        Only meaningful for 2-D sites (label == 1).  For 1-D / 3-D
        observations the returned value is ``NaN``.

        Parameters
        ----------
        X : ndarray (n_samples, 6) or DataFrame

        Returns
        -------
        strike : float ndarray, shape (n_samples,)
            Estimated strike angles in degrees.  ``NaN`` for non-2-D sites.
        """
        if self._use_rf or self._network is None:
            # RF backend has no strike head — return NaN
            n = len(self._coerce_X(X))
            return np.full(n, np.nan)

        try:
            import torch
        except ImportError:
            n = len(self._coerce_X(X))
            return np.full(n, np.nan)

        X_arr = self._coerce_X(X)
        Xn = (X_arr - self._x_mean) / self._x_std
        dev = next(self._network.parameters()).device
        self._network.eval()
        with torch.no_grad():
            t = torch.from_numpy(Xn.astype(np.float32)).to(dev)
            _, strike_raw = self._network(t)
            strike = strike_raw.cpu().numpy()

        # Mask non-2D predictions with NaN
        labels = self._predict_proba(Xn).argmax(axis=1)
        strike[labels != 1] = np.nan
        return strike

    def predict_table(self, sites: Any) -> pd.DataFrame:
        """
        Classify an entire site collection and return a result DataFrame.

        Parameters
        ----------
        sites : SiteCollection or compatible

        Returns
        -------
        df : DataFrame
            Columns: station, freq, period, dim (0/1/2),
            dim_label (str), strike (°), confidence.
        """
        try:
            from pycsamt.emtools.dimensionality import phase_features_table
        except ImportError as exc:
            raise ImportError("emtools is required for predict_table") from exc

        df = phase_features_table(sites)
        if df.empty:
            return df

        X_arr = _df_to_feature_matrix(df)
        labels = self.predict(X_arr)
        proba = self.transform(X_arr)
        strike = self.predict_strike(X_arr)

        _label_map = {0: "1D", 1: "2D", 2: "3D"}
        out = df[["station", "freq", "period"]].copy()
        out["dim"] = labels
        out["dim_label"] = [_label_map.get(d, "?") for d in labels]
        out["strike"] = strike
        out["confidence"] = proba.max(axis=1)
        return out

    # ─── internal ─────────────────────────────────────────────────────────

    def _fit_torch(
        self, Xn, y, strike, *, epochs, batch_size, lr, val_frac, seed, verbose
    ) -> None:
        import torch
        import torch.nn as nn
        from torch.utils.data import TensorDataset, DataLoader

        rng = np.random.default_rng(seed)
        dev = self._resolve_device()
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        vi, ti = idx[:n_val], idx[n_val:]

        has_strike = strike is not None and np.any(np.isfinite(strike))

        def _mk_tensors(idx_):
            t_x = torch.from_numpy(Xn[idx_].astype(np.float32))
            t_y = torch.from_numpy(y[idx_].astype(np.int64))
            if has_strike:
                t_s = torch.from_numpy(strike[idx_].astype(np.float32))
            else:
                t_s = torch.zeros(len(idx_))
            return TensorDataset(t_x, t_y, t_s)

        tr_ds = _mk_tensors(ti)
        ce = nn.CrossEntropyLoss()
        mse = nn.MSELoss()

        self._network = _build_dim_mlp(
            _N_FEATURES, self.n_classes, self.hidden, self.dropout
        ).to(dev)
        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=8, min_lr=1e-6
        )

        Xva = torch.from_numpy(Xn[vi].astype(np.float32)).to(dev)
        yva = torch.from_numpy(y[vi].astype(np.int64)).to(dev)

        best_val, best_state = np.inf, None
        train_losses, val_losses = [], []

        for ep in range(1, epochs + 1):
            self._network.train()
            loader = DataLoader(tr_ds, batch_size=batch_size, shuffle=True)
            ep_loss = 0.0
            for xb, yb, sb in loader:
                xb, yb, sb = xb.to(dev), yb.to(dev), sb.to(dev)
                cls_out, str_out = self._network(xb)
                loss = ce(cls_out, yb)
                if has_strike:
                    is2d = (yb == 1)
                    if is2d.any():
                        s_mask = sb[is2d]
                        valid_s = torch.isfinite(s_mask)
                        if valid_s.any():
                            loss = loss + 0.1 * mse(
                                str_out[is2d][valid_s], s_mask[valid_s]
                            )
                opt.zero_grad()
                loss.backward()
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= len(ti)

            self._network.eval()
            with torch.no_grad():
                cls_v, _ = self._network(Xva)
                v_loss = ce(cls_v, yva).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val:
                best_val = v_loss
                import copy
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 5) == 0 or ep == 1):
                acc = (cls_v.argmax(dim=1) == yva).float().mean().item()
                print(f"  DimClassifier  ep {ep:>4d}/{epochs}  "
                      f"loss={ep_loss:.4f}  val_loss={v_loss:.4f}  "
                      f"val_acc={acc:.3f}")

        if best_state is not None:
            self._network.load_state_dict(best_state)
        self._history = {"train_loss": train_losses, "val_loss": val_losses}

    def _fit_rf(self, Xn: np.ndarray, y: np.ndarray, *, verbose: bool) -> None:
        try:
            from sklearn.ensemble import RandomForestClassifier
        except ImportError as exc:
            raise ImportError(
                "Either PyTorch or scikit-learn is required for "
                "DimensionalityClassifier"
            ) from exc

        self._rf = RandomForestClassifier(n_estimators=200, random_state=0)
        valid = np.all(np.isfinite(Xn), axis=1)
        self._rf.fit(Xn[valid], y[valid])
        if verbose:
            print("  DimClassifier (RandomForest fallback) fitted.")

    def _predict_proba(self, Xn: np.ndarray) -> np.ndarray:
        if self._use_rf:
            valid = np.all(np.isfinite(Xn), axis=1)
            out = np.full((len(Xn), self.n_classes), 1.0 / self.n_classes)
            if valid.any():
                out[valid] = self._rf.predict_proba(Xn[valid])
            return out
        import torch
        dev = next(self._network.parameters()).device
        self._network.eval()
        with torch.no_grad():
            t = torch.from_numpy(Xn.astype(np.float32)).to(dev)
            cls_out, _ = self._network(t)
            proba = torch.softmax(cls_out, dim=1).cpu().numpy()
        return proba

    def _coerce_X(self, X) -> np.ndarray:
        if isinstance(X, pd.DataFrame):
            return _df_to_feature_matrix(X)
        return np.asarray(X, dtype=np.float32)

    def _coerce_Xy(self, X, y, strike):
        if isinstance(X, pd.DataFrame):
            X_arr, y_arr, strike_arr = _df_to_Xy(X, "dim", None, 3.0, 0.2)
            if y is not None:
                y_arr = np.asarray(y, dtype=int)
            if strike is not None:
                strike_arr = np.asarray(strike, dtype=float)
        else:
            X_arr = np.asarray(X, dtype=np.float32)
            X_arr = np.where(np.isfinite(X_arr), X_arr, 0.0)
            if y is None:
                raise ValueError("y (labels) must be provided when X is an array.")
            y_arr = np.asarray(y, dtype=int)
            strike_arr = (
                np.asarray(strike, dtype=float)
                if strike is not None else None
            )
        return X_arr, y_arr, strike_arr

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

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> Dict[str, Any]:
        return {
            "hidden": list(self.hidden),
            "dropout": self.dropout,
            "n_classes": self.n_classes,
            "lr": self.lr,
            "device": self.device,
        }

    def _get_weights(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        if self._x_mean is not None:
            out["_x_mean"] = self._x_mean
            out["_x_std"] = self._x_std
        if self._network is not None:
            for k, v in self._network.state_dict().items():
                out[k] = v.cpu().numpy()
        elif self._rf is not None:
            try:
                import pickle, io
                buf = io.BytesIO()
                pickle.dump(self._rf, buf)
                out["_rf_pickle"] = np.frombuffer(buf.getvalue(), dtype=np.uint8)
            except Exception:
                pass
        return out

    def _load_weights(self, weights: Dict[str, np.ndarray]) -> None:
        self._x_mean = weights.pop("_x_mean", None)
        self._x_std = weights.pop("_x_std", None)
        rf_blob = weights.pop("_rf_pickle", None)
        if rf_blob is not None:
            try:
                import pickle, io
                self._rf = pickle.load(io.BytesIO(bytes(rf_blob)))
                self._use_rf = True
                self._is_fitted = True
                return
            except Exception:
                pass
        if weights:
            import torch
            self._network = _build_dim_mlp(
                _N_FEATURES, self.n_classes, self.hidden, self.dropout
            )
            state = {k: torch.from_numpy(v) for k, v in weights.items()}
            self._network.load_state_dict(state)
        self._is_fitted = True

    def __repr__(self) -> str:
        backend = "rf" if self._use_rf else "torch"
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"DimensionalityClassifier(n_classes={self.n_classes}, "
            f"{backend}, {status})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# DataFrame helpers
# ─────────────────────────────────────────────────────────────────────────────

def _df_to_feature_matrix(df: pd.DataFrame) -> np.ndarray:
    """Extract feature matrix from a phase_features_table DataFrame."""
    present = [c for c in _FEATURE_COLS if c in df.columns]
    missing = [c for c in _FEATURE_COLS if c not in df.columns]

    mat = np.full((len(df), _N_FEATURES), 0.0, dtype=np.float32)
    for ci, col in enumerate(_FEATURE_COLS):
        if col in df.columns:
            vals = df[col].to_numpy(dtype=float)
            mat[:, ci] = np.where(np.isfinite(vals), vals, 0.0)
    return mat


def _df_to_Xy(
    df: pd.DataFrame,
    label_col: Optional[str],
    strike_col: Optional[str],
    skew_th: float,
    ellipt_th: float,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    X = _df_to_feature_matrix(df)

    if label_col is not None and label_col in df.columns:
        y = df[label_col].to_numpy(dtype=int)
    else:
        beta = df["beta_abs"].to_numpy(dtype=float) if "beta_abs" in df.columns else np.zeros(len(df))
        ellipt = df["ellipt_abs"].to_numpy(dtype=float) if "ellipt_abs" in df.columns else np.zeros(len(df))
        y = _rule_labels(beta, ellipt, skew_th, ellipt_th)

    strike: Optional[np.ndarray] = None
    if strike_col is not None and strike_col in df.columns:
        strike = df[strike_col].to_numpy(dtype=float)

    return X, y, strike

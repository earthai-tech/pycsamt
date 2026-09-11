# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
DistortionTypeClassifier — learned galvanic-distortion triage.

Why this does not re-solve the physics
------------------------------------------
:mod:`pycsamt.emtools.ss` already estimates static shift with four
different spatial-statistics methods (AMA, LOESS, bilateral filter,
reference-median), and :mod:`pycsamt.emtools.gb` already fits a full
Groom-Bailey galvanic-distortion decomposition
(:func:`~pycsamt.emtools.gb.groom_bailey_table`,
:func:`~pycsamt.emtools.gb.apply_groom_bailey`) by nonlinear least
squares. ``DistortionTypeClassifier`` is not meant to replace either
-- it is meant to *triage*: given a station's phase-tensor invariants
and its already-fitted Groom-Bailey parameters, label it as clean,
static-shift-only, or in need of the full decomposition, the same way
:class:`~pycsamt.ai.processing.classify.DimensionalityClassifier`
refines :func:`~pycsamt.emtools.dimensionality.classify_dimensionality`'s
threshold rule into a smoother multi-feature classifier rather than
inventing a new definition of dimensionality. The value is routing a
survey's many stations toward the right *existing* tool quickly, not
a new distortion model.

Classes
-------
0. **clean** -- resistivity close to the along-line spatial trend,
   twist and shear close to 0; no correction needed.
1. **static-shift-only** -- resistivity departs from the spatial trend
   but twist and shear stay small: a pure multiplicative offset, the
   case :mod:`pycsamt.emtools.ss` targets directly.
2. **distorted** -- twist and/or shear are non-negligible: rotation-
   and anisotropy-like distortion that a static-shift correction alone
   cannot fix, the case :func:`~pycsamt.emtools.gb.apply_groom_bailey`
   targets.

Feature vector (per station)
--------------------------------
``[beta_abs, ellipt_abs, delta_log10_rho, twist_deg, shear,
anisotropy]`` -- ``beta_abs``/``ellipt_abs`` from
:func:`~pycsamt.emtools.dimensionality.phase_features_table` (median
over frequency, already reused by ``DimensionalityClassifier``);
``twist_deg``/``shear``/``anisotropy`` from
:func:`~pycsamt.emtools.gb.groom_bailey_table`, frequency-independent
by construction (Groom-Bailey assumes one real distortion matrix per
station); ``delta_log10_rho``, the station's log10-resistivity
deviation from the spatial trend along the line, from
:func:`~pycsamt.emtools.ss.estimate_ss_ama`.
:func:`build_distortion_features_table` builds this table directly
from a site collection.

.. note::

   Groom-Bailey's own fitted distortion matrix is normalised to unit
   determinant at every iteration
   (:func:`~pycsamt.emtools.gb._normalise_distortion`) -- the textbook
   convention, since the absolute gain is degenerate with the unknown
   regional resistivity and genuinely unrecoverable from one station's
   data alone. Its ``gain`` column is therefore always exactly ``1.0``
   and carries no information about static shift; ``delta_log10_rho``,
   a real cross-station spatial comparison, is used here instead.

Self-training labels
------------------------
Mirroring ``DimensionalityClassifier``'s ``_rule_labels`` /
``from_features_table`` pattern: default training labels come from
simple thresholds on ``|twist_deg|``, ``|shear|``, and
``|delta_log10_rho|`` when no ``label_col`` is supplied, so the
network starts out approximating a transparent rule and can be refined
with real labels later exactly the way ``DimensionalityClassifier`` is.
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np
import pandas as pd

from .._backend_utils import (
    active_backend,
    get_weights,
    resolve_device,
    set_weights,
)
from .._base import BaseEMProcessor

__all__ = ["DistortionTypeClassifier", "build_distortion_features_table"]

_FEATURE_COLS = [
    "beta_abs",
    "ellipt_abs",
    "delta_log10_rho",
    "twist_deg",
    "shear",
    "anisotropy",
]
_N_FEATURES = len(_FEATURE_COLS)
_DISTORTION_LABELS = ("clean", "static_shift_only", "distorted")


# ─────────────────────────────────────────────────────────────────────────────
# Feature-table builder (no torch / TF dependency)
# ─────────────────────────────────────────────────────────────────────────────


def build_distortion_features_table(
    sites: Any,
    *,
    band: tuple[float, float] | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Build the per-station feature table :class:`DistortionTypeClassifier`
    expects.

    Aggregates
    :func:`~pycsamt.emtools.dimensionality.phase_features_table`'s
    per-(station, frequency) ``beta_abs`` / ``ellipt_abs`` to a
    per-station median, then merges with
    :func:`~pycsamt.emtools.gb.groom_bailey_table`'s per-station
    ``twist_deg`` / ``shear`` / ``anisotropy`` and
    :func:`~pycsamt.emtools.ss.estimate_ss_ama`'s per-station
    ``delta_log10_rho`` -- one row per station present in all three.
    Groom-Bailey's own ``gain`` column is *not* used; see the module
    docstring for why.

    Parameters
    ----------
    sites : SiteCollection or compatible
    band : (f_lo, f_hi) or None
        Frequency band passed to
        :func:`~pycsamt.emtools.gb.groom_bailey_table`.
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    df : DataFrame
        Columns: ``station``, ``beta_abs``, ``ellipt_abs``,
        ``delta_log10_rho``, ``twist_deg``, ``shear``, ``anisotropy``.
        Empty (but correctly columned) when no station has estimates
        from all three sources.
    """
    try:
        from pycsamt.emtools.dimensionality import phase_features_table
        from pycsamt.emtools.gb import groom_bailey_table
        from pycsamt.emtools.ss import estimate_ss_ama
    except ImportError as exc:
        raise ImportError(
            "emtools is required for build_distortion_features_table"
        ) from exc

    cols = ["station", *_FEATURE_COLS]

    phase_df = phase_features_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if phase_df.empty:
        return pd.DataFrame(columns=cols)

    phase_agg = (
        phase_df.groupby("station")[["beta_abs", "ellipt_abs"]]
        .median()
        .reset_index()
    )

    gb_df = groom_bailey_table(
        sites,
        band=band,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if gb_df.empty or "status" not in gb_df.columns:
        return pd.DataFrame(columns=cols)
    gb_ok = gb_df[gb_df["status"] == "ok"]
    if gb_ok.empty:
        return pd.DataFrame(columns=cols)

    ss_df = estimate_ss_ama(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ss_df.empty:
        return pd.DataFrame(columns=cols)

    merged = phase_agg.merge(
        gb_ok[["station", "twist_deg", "shear", "anisotropy"]],
        on="station",
        how="inner",
    ).merge(
        ss_df[["station", "delta_log10_rho"]],
        on="station",
        how="inner",
    )
    return merged[cols]


# ─────────────────────────────────────────────────────────────────────────────
# Rule-based label generation (for self-training)
# ─────────────────────────────────────────────────────────────────────────────


def _rule_labels(
    shift: np.ndarray,
    twist_deg: np.ndarray,
    shear: np.ndarray,
    shift_th: float = 0.1,
    twist_th: float = 10.0,
    shear_th: float = 0.1,
) -> np.ndarray:
    """
    Threshold-based distortion-regime label, matching the module
    docstring's class definitions.

    ``shift_th`` is in log10 units on ``delta_log10_rho`` (0.1 ->
    roughly a 25% resistivity deviation from the spatial trend);
    ``twist_th`` in degrees; ``shear_th`` on the dimensionless shear
    ratio.
    """
    low_rotation = (np.abs(twist_deg) <= twist_th) & (
        np.abs(shear) <= shear_th
    )
    labels = np.full(len(shift), 2, dtype=int)  # distorted (default)
    labels[low_rotation & (np.abs(shift) <= shift_th)] = 0  # clean
    labels[low_rotation & (np.abs(shift) > shift_th)] = 1  # static shift
    return labels


# ─────────────────────────────────────────────────────────────────────────────
# Network builders
# ─────────────────────────────────────────────────────────────────────────────


def _build_distortion_mlp_torch(
    n_features: int,
    n_classes: int,
    hidden: tuple[int, ...],
    dropout: float,
) -> Any:
    """MLP classifier -- PyTorch."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError(
            "PyTorch is required for DistortionTypeClassifier"
        ) from exc

    dims = [n_features, *hidden]

    class _DistortionMLP(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            layers: list = []
            for i in range(len(dims) - 1):
                layers += [
                    nn.Linear(dims[i], dims[i + 1]),
                    nn.BatchNorm1d(dims[i + 1]),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                ]
            self.backbone = nn.Sequential(*layers)
            self.cls_head = nn.Linear(dims[-1], n_classes)

        def forward(self, x):
            return self.cls_head(self.backbone(x))

    return _DistortionMLP()


def _build_distortion_mlp_tf(
    n_features: int,
    n_classes: int,
    hidden: tuple[int, ...],
    dropout: float,
) -> Any:
    """MLP classifier -- TensorFlow/Keras."""
    try:
        import tensorflow as tf
        from tensorflow.keras import Model, layers
    except ImportError as exc:
        raise ImportError(
            "TensorFlow is required for DistortionTypeClassifier "
            "(TF backend)"
        ) from exc

    inp = tf.keras.Input(shape=(n_features,), name="input")
    x = inp
    for h in hidden:
        x = layers.Dense(h)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
        x = layers.Dropout(dropout)(x)
    out = layers.Dense(n_classes, name="cls")(x)

    return Model(inp, out, name="distortion_mlp")


# ─────────────────────────────────────────────────────────────────────────────
# DistortionTypeClassifier
# ─────────────────────────────────────────────────────────────────────────────


class DistortionTypeClassifier(BaseEMProcessor):
    """
    Learned triage classifier for galvanic-distortion regime.

    Parameters
    ----------
    hidden : tuple of int, default (32, 16)
        Hidden-layer widths of the shared MLP backbone -- deliberately
        smaller than
        :class:`~pycsamt.ai.processing.classify.DimensionalityClassifier`'s
        default, since this problem is station-level (tens to low
        hundreds of samples per survey) rather than
        station-frequency-level.
    dropout : float, default 0.2
        Dropout probability in each hidden layer.
    n_classes : int, default 3
        Number of distortion-regime classes; see the module docstring.
    lr : float, default 1e-3
        Default learning rate; can be overridden in :meth:`fit`.
    device : str or None

    Notes
    -----
    Falls back to a random-forest classifier (scikit-learn) when
    neither PyTorch nor TensorFlow is available; :meth:`transform`
    then returns the forest's own class probabilities.

    Examples
    --------
    >>> from pycsamt.ai.processing.distortion import (
    ...     DistortionTypeClassifier, build_distortion_features_table,
    ... )
    >>> feats = build_distortion_features_table(sites)  # doctest: +SKIP
    >>> clf = DistortionTypeClassifier.from_features_table(
    ...     feats, epochs=80,
    ... )  # doctest: +SKIP
    >>> table = clf.predict_table(sites)  # doctest: +SKIP
    """

    def __init__(
        self,
        hidden: tuple[int, ...] = (32, 16),
        dropout: float = 0.2,
        n_classes: int = 3,
        lr: float = 1e-3,
        device: str | None = None,
    ) -> None:
        self.hidden = tuple(hidden)
        self.dropout = float(dropout)
        self.n_classes = int(n_classes)
        self.lr = float(lr)
        self.device = device

        self._network: Any = None
        self._rf: Any = None
        self._use_rf: bool = False
        self._backend_name: str | None = None
        self._x_mean: np.ndarray | None = None
        self._x_std: np.ndarray | None = None
        self._is_fitted: bool = False
        self._history: dict[str, list] = {}

    # ─── factory from emtools ─────────────────────────────────────────────

    @classmethod
    def from_features_table(
        cls,
        df: pd.DataFrame,
        *,
        label_col: str | None = None,
        shift_th: float = 0.1,
        twist_th: float = 10.0,
        shear_th: float = 0.1,
        **fit_kwargs,
    ) -> DistortionTypeClassifier:
        """
        Construct and train a classifier from a
        :func:`build_distortion_features_table` DataFrame.

        Parameters
        ----------
        df : DataFrame
        label_col : str or None
            Column for pre-computed labels; ``None`` -- the default --
            falls back to rule-based self-training labels (see the
            module docstring).
        shift_th, twist_th, shear_th : float
            Rule-based thresholds (used when ``label_col`` is absent).
        **fit_kwargs
            Passed to :meth:`fit` (e.g. ``epochs=80``).

        Returns
        -------
        DistortionTypeClassifier
        """
        X, y = _df_to_Xy(df, label_col, shift_th, twist_th, shear_th)
        obj = cls()
        obj.fit(X, y, **fit_kwargs)
        return obj

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(
        self,
        X: np.ndarray | pd.DataFrame,
        y: np.ndarray | None = None,
        *,
        epochs: int = 80,
        batch_size: int = 64,
        lr: float | None = None,
        val_frac: float = 0.15,
        seed: int | None = None,
        verbose: bool = True,
    ) -> DistortionTypeClassifier:
        """
        Train the distortion-triage classifier.

        Parameters
        ----------
        X : ndarray (n_samples, 6) or DataFrame
        y : int ndarray (n_samples,) or None
            Class labels (0=clean, 1=static-shift-only, 2=distorted).
            ``None`` falls back to the rule-based self-training labels
            described in the module docstring.
        epochs, batch_size, lr, val_frac, seed, verbose
            Training hyper-parameters. ``epochs`` defaults low (``80``)
            for a quick look; on a real, small survey (tens of
            stations, as most are) this is often too few for the
            self-trained network to converge reliably -- observed on
            28-station data, repeated fits at the default can disagree
            on which regime is the *majority* class, not just on a
            minority label appearing or not. Prefer more epochs (e.g.
            ``150``-``200``) and treat a single run's
            :meth:`predict_table` output as one plausible smoothing of
            the rule boundary rather than a converged answer -- see
            the user guide for a concrete before/after comparison.

        Returns
        -------
        self
        """
        X_arr, y_arr = self._coerce_Xy(X, y)

        self._x_mean = X_arr.mean(axis=0, keepdims=True)
        self._x_std = X_arr.std(axis=0, keepdims=True) + 1e-8
        Xn = (X_arr - self._x_mean) / self._x_std
        _lr = float(lr) if lr is not None else self.lr

        try:
            self._backend_name = active_backend()
            if self._backend_name == "tensorflow":
                self._fit_tensorflow(
                    Xn,
                    y_arr,
                    epochs=epochs,
                    batch_size=batch_size,
                    lr=_lr,
                    val_frac=val_frac,
                    seed=seed,
                    verbose=verbose,
                )
            else:
                self._fit_torch(
                    Xn,
                    y_arr,
                    epochs=epochs,
                    batch_size=batch_size,
                    lr=_lr,
                    val_frac=val_frac,
                    seed=seed,
                    verbose=verbose,
                )
            self._use_rf = False
        except (RuntimeError, ImportError):
            self._backend_name = "rf"
            self._fit_rf(Xn, y_arr, verbose=verbose)
            self._use_rf = True

        self._is_fitted = True
        return self

    def transform(self, X: np.ndarray | pd.DataFrame) -> np.ndarray:
        """
        Compute class probabilities.

        Returns
        -------
        proba : ndarray, shape (n_samples, n_classes)
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")
        Xn = (self._coerce_X(X) - self._x_mean) / self._x_std
        return self._predict_proba(Xn)

    def predict(self, X: np.ndarray | pd.DataFrame) -> np.ndarray:
        """Predict the distortion regime (0=clean, 1=static-shift-only,
        2=distorted)."""
        return self.transform(X).argmax(axis=1)

    def predict_table(
        self,
        sites: Any,
        *,
        band: tuple[float, float] | None = None,
        recursive: bool = True,
        on_dup: str = "replace",
        strict: bool = False,
        verbose: int = 0,
    ) -> pd.DataFrame:
        """
        Classify an entire site collection and return a result
        DataFrame.

        Parameters
        ----------
        sites : SiteCollection or compatible
        band : (f_lo, f_hi) or None
            Passed to :func:`build_distortion_features_table`.
        recursive, on_dup, strict, verbose
            Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

        Returns
        -------
        df : DataFrame
            Columns: ``station``, the six input features, ``regime``
            (0/1/2), ``regime_label`` (str), ``confidence``.
        """
        df = build_distortion_features_table(
            sites,
            band=band,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        if df.empty:
            return df

        X_arr = _df_to_feature_matrix(df)
        labels = self.predict(X_arr)
        proba = self.transform(X_arr)

        out = df.copy()
        out["regime"] = labels
        out["regime_label"] = [_DISTORTION_LABELS[label] for label in labels]
        out["confidence"] = proba.max(axis=1)
        return out

    # ─── internal training paths ──────────────────────────────────────────

    def _fit_torch(
        self,
        Xn: np.ndarray,
        y: np.ndarray,
        *,
        epochs: int,
        batch_size: int,
        lr: float,
        val_frac: float,
        seed: int | None,
        verbose: bool,
    ) -> None:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset

        rng = np.random.default_rng(seed)
        dev = resolve_device(self.device)
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        vi, ti = idx[:n_val], idx[n_val:]

        tr_ds = TensorDataset(
            torch.from_numpy(Xn[ti].astype(np.float32)),
            torch.from_numpy(y[ti].astype(np.int64)),
        )
        ce = nn.CrossEntropyLoss()

        self._network = _build_distortion_mlp_torch(
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
            ep_loss = 0.0
            for xb, yb in DataLoader(
                tr_ds, batch_size=batch_size, shuffle=True
            ):
                xb, yb = xb.to(dev), yb.to(dev)
                out = self._network(xb)
                loss = ce(out, yb)
                opt.zero_grad()
                loss.backward()
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= len(ti)

            self._network.eval()
            with torch.no_grad():
                v_out = self._network(Xva)
                v_loss = ce(v_out, yva).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val:
                best_val = v_loss
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 5) == 0 or ep == 1):
                acc = (v_out.argmax(dim=1) == yva).float().mean().item()
                print(
                    f"  DistortionClassifier  ep {ep:>4d}/{epochs}  "
                    f"loss={ep_loss:.4f}  val_loss={v_loss:.4f}  "
                    f"val_acc={acc:.3f}"
                )

        if best_state is not None:
            self._network.load_state_dict(best_state)
        self._history = {"train_loss": train_losses, "val_loss": val_losses}

    def _fit_tensorflow(
        self,
        Xn: np.ndarray,
        y: np.ndarray,
        *,
        epochs: int,
        batch_size: int,
        lr: float,
        val_frac: float,
        seed: int | None,
        verbose: bool,
    ) -> None:
        import tensorflow as tf

        rng = np.random.default_rng(seed)
        n = len(Xn)
        idx = rng.permutation(n)
        n_val = max(1, int(n * val_frac))
        vi, ti = idx[:n_val], idx[n_val:]

        Xtr = Xn[ti].astype(np.float32)
        ytr = y[ti].astype(np.int64)
        Xva = Xn[vi].astype(np.float32)
        yva = y[vi].astype(np.int64)

        dev = resolve_device(self.device)
        with tf.device(dev):
            self._network = _build_distortion_mlp_tf(
                _N_FEATURES, self.n_classes, self.hidden, self.dropout
            )
            self._network.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
                loss=tf.keras.losses.SparseCategoricalCrossentropy(
                    from_logits=True
                ),
                metrics=["accuracy"],
            )
            hist = self._network.fit(
                Xtr,
                ytr,
                validation_data=(Xva, yva),
                epochs=epochs,
                batch_size=batch_size,
                callbacks=[
                    tf.keras.callbacks.EarlyStopping(
                        monitor="val_loss",
                        patience=12,
                        restore_best_weights=True,
                        min_delta=1e-6,
                    ),
                    tf.keras.callbacks.ReduceLROnPlateau(
                        monitor="val_loss",
                        factor=0.5,
                        patience=8,
                        min_lr=1e-6,
                    ),
                ],
                verbose=1 if verbose else 0,
            )
        self._history = {
            "train_loss": hist.history["loss"],
            "val_loss": hist.history.get("val_loss", []),
        }

    def _fit_rf(self, Xn: np.ndarray, y: np.ndarray, *, verbose: bool) -> None:
        try:
            from sklearn.ensemble import RandomForestClassifier
        except ImportError as exc:
            raise ImportError(
                "PyTorch, TensorFlow, or scikit-learn is required for "
                "DistortionTypeClassifier"
            ) from exc

        self._rf = RandomForestClassifier(n_estimators=200, random_state=0)
        valid = np.all(np.isfinite(Xn), axis=1)
        self._rf.fit(Xn[valid], y[valid])
        if verbose:
            print("  DistortionTypeClassifier (RandomForest fallback) fitted.")

    def _predict_proba(self, Xn: np.ndarray) -> np.ndarray:
        if self._use_rf:
            valid = np.all(np.isfinite(Xn), axis=1)
            out = np.full((len(Xn), self.n_classes), 1.0 / self.n_classes)
            if valid.any():
                # RandomForestClassifier.predict_proba only returns a
                # column per class actually present in the training
                # labels -- with this few stations and three classes,
                # one regime (usually "clean") is often entirely
                # absent, so map columns back by self._rf.classes_
                # rather than assuming they span range(n_classes).
                proba = self._rf.predict_proba(Xn[valid])
                filled = np.zeros((proba.shape[0], self.n_classes))
                for j, c in enumerate(self._rf.classes_):
                    filled[:, int(c)] = proba[:, j]
                out[valid] = filled
            return out

        if self._backend_name == "tensorflow":
            logits = self._network.predict(Xn.astype(np.float32), verbose=0)
            e = np.exp(logits - logits.max(axis=1, keepdims=True))
            return e / e.sum(axis=1, keepdims=True)

        import torch

        dev = next(self._network.parameters()).device
        self._network.eval()
        with torch.no_grad():
            t = torch.from_numpy(Xn.astype(np.float32)).to(dev)
            out = self._network(t)
            proba = torch.softmax(out, dim=1).cpu().numpy()
        return proba

    def _coerce_X(self, X) -> np.ndarray:
        if isinstance(X, pd.DataFrame):
            return _df_to_feature_matrix(X)
        return np.asarray(X, dtype=np.float32)

    def _coerce_Xy(self, X, y):
        if isinstance(X, pd.DataFrame):
            X_arr, y_arr = _df_to_Xy(X, None)
            if y is not None:
                y_arr = np.asarray(y, dtype=int)
        else:
            X_arr = np.asarray(X, dtype=np.float32)
            X_arr = np.where(np.isfinite(X_arr), X_arr, 0.0)
            if y is None:
                if X_arr.shape[1] >= 5:
                    # columns: beta_abs, ellipt_abs, delta_log10_rho,
                    # twist_deg, shear
                    y_arr = _rule_labels(
                        X_arr[:, 2], X_arr[:, 3], X_arr[:, 4]
                    )
                else:
                    y_arr = np.zeros(len(X_arr), dtype=int)
            else:
                y_arr = np.asarray(y, dtype=int)
        return X_arr, y_arr

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "hidden": list(self.hidden),
            "dropout": self.dropout,
            "n_classes": self.n_classes,
            "lr": self.lr,
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
        elif self._rf is not None:
            try:
                import io
                import pickle

                buf = io.BytesIO()
                pickle.dump(self._rf, buf)
                out["_rf_pickle"] = np.frombuffer(
                    buf.getvalue(), dtype=np.uint8
                )
            except Exception:
                pass
        return out

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        self._x_mean = weights.pop("_x_mean", None)
        self._x_std = weights.pop("_x_std", None)
        backend_blob = weights.pop("_backend", None)
        self._backend_name = (
            str(backend_blob) if backend_blob is not None else "torch"
        )

        rf_blob = weights.pop("_rf_pickle", None)
        if rf_blob is not None:
            try:
                import io
                import pickle

                self._rf = pickle.load(io.BytesIO(bytes(rf_blob)))
                self._use_rf = True
                self._is_fitted = True
                return
            except Exception:
                pass

        if weights:
            if self._backend_name == "tensorflow":
                self._network = _build_distortion_mlp_tf(
                    _N_FEATURES, self.n_classes, self.hidden, self.dropout
                )
            else:
                self._network = _build_distortion_mlp_torch(
                    _N_FEATURES, self.n_classes, self.hidden, self.dropout
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
            epoch. Empty when the random-forest fallback was used (no
            epoch loop) or before :meth:`fit` has been called. Pass
            directly to
            :func:`~pycsamt.ai.processing.plot.plot_training_history`.
        """
        return dict(self._history)

    def __repr__(self) -> str:
        backend = self._backend_name or ("rf" if self._use_rf else "torch")
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"DistortionTypeClassifier(n_classes={self.n_classes}, "
            f"{backend}, {status})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# DataFrame helpers
# ─────────────────────────────────────────────────────────────────────────────


def _df_to_feature_matrix(df: pd.DataFrame) -> np.ndarray:
    mat = np.full((len(df), _N_FEATURES), 0.0, dtype=np.float32)
    for ci, col in enumerate(_FEATURE_COLS):
        if col in df.columns:
            vals = df[col].to_numpy(dtype=float)
            mat[:, ci] = np.where(np.isfinite(vals), vals, 0.0)
    return mat


def _df_to_Xy(
    df: pd.DataFrame,
    label_col: str | None,
    shift_th: float = 0.1,
    twist_th: float = 10.0,
    shear_th: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    X = _df_to_feature_matrix(df)

    if label_col is not None and label_col in df.columns:
        y = df[label_col].to_numpy(dtype=int)
    else:
        shift = (
            df["delta_log10_rho"].to_numpy(dtype=float)
            if "delta_log10_rho" in df.columns
            else np.zeros(len(df))
        )
        twist = (
            df["twist_deg"].to_numpy(dtype=float)
            if "twist_deg" in df.columns
            else np.zeros(len(df))
        )
        shear = (
            df["shear"].to_numpy(dtype=float)
            if "shear" in df.columns
            else np.zeros(len(df))
        )
        y = _rule_labels(shift, twist, shear, shift_th, twist_th, shear_th)

    return X, y

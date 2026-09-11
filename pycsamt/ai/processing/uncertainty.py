# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
UncertaintyCalibrator — learned per-cell error-floor re-estimation.

Why this is not EMQCScorer
----------------------------
:class:`~pycsamt.ai.processing.qc.EMQCScorer` answers "is this
(station, frequency) cell trustworthy?" with a score in ``[0, 1]``. It
does not answer "how large should the error bar on this cell be?" --
yet that second question is exactly what an inversion error floor
needs, and today pyCSAMT (like most MT software) uses whatever
``z_err`` field the processing software wrote, or a flat
percentage-of-|Z| floor. Both are frequently miscalibrated: field
error propagation is often too optimistic in quiet bands and too
pessimistic near noise spikes. ``UncertaintyCalibrator`` re-estimates
that error directly, using the same per-cell diagnostic features
:class:`~pycsamt.ai.processing.qc.EMQCScorer` already extracts (SNR,
Swift skew, off-diagonal amplitude asymmetry, both off-diagonal
phases) as regression inputs instead of classification inputs.

Training target: recalibration mode
--------------------------------------
There is no pre-existing "correct error bar" to regress against, the
way ``EMQCScorer``'s rule-based labels or ``DimensionalityClassifier``'s
skew/ellipticity thresholds exist. Two options were considered:

1. **Recalibration mode** (implemented here) -- train against the
   field-processing ``z_err`` itself, expressed as a *fractional*
   error :math:`|z_\\text{err}| / |Z|` rather than an absolute value,
   so the target is comparable across stations and frequencies with
   very different :math:`|Z|` (the raw ``z_err`` column on L18PLT
   spans three decades, 1.3 to 1625, purely from :math:`|Z|` scale;
   the fractional error is a tight, physically meaningful 1.1% to
   54%, median 4.9%, genuinely varying across 1428 distinct values --
   confirmed non-degenerate before committing to it, the same check
   that caught Groom-Bailey's always-1.0 ``gain`` column being unsafe
   as a feature in
   :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`).
   The network is not memorising ``z_err`` verbatim -- it only sees
   the five diagnostic features, so its output is a smoothed,
   internally consistent re-estimate driven by signal quality, not a
   pass-through copy.
2. **Empirical mode** -- train against a genuine bootstrap or
   repeat-measurement spread computed from raw cross-power spectra
   (:mod:`pycsamt.seg.spectra`) when available, rather than a
   recalibrated copy of the field estimate. Not implemented here: none
   of pyCSAMT's bundled *field-survey* EDI data carries an ``>SPECTRA``
   section (only four small vendor demo files under
   ``data/MT/SPECTRA/`` do, with placeholder ``COUNTRY="EXAMPLE"``
   metadata -- not a real profile, so not a fair verification target
   for this package's real-data standard). :meth:`fit` still accepts
   any ``y`` explicitly, so an empirical target computed independently
   -- from :mod:`pycsamt.seg.spectra`, or anywhere else -- plugs into
   the exact same regressor without an API change; only the
   *automatic* ``y=None`` default is recalibration-mode.

Feature vector (per station-frequency cell)
------------------------------------------------
``[swift_skew, asym, phase_xy, phase_yx]`` -- four of
:class:`~pycsamt.ai.processing.qc.EMQCScorer`'s five features, built
by the same :func:`~pycsamt.ai.processing.qc._extract_qc_features`.

.. note::

   ``EMQCScorer``'s fifth feature, ``snr``, is deliberately **not**
   used as a model input here, even though
   :func:`build_uncertainty_features_table` still reports it as a
   diagnostic column. ``snr = amp / (z_err + eps)`` and the
   recalibration target ``z_err_frac = z_err / (amp + eps)`` are built
   from the exact same pooled ``amp``/``z_err`` terms -- verified on
   real L18PLT data, ``corr(1/snr, z_err_frac) = 1 - 2e-16``, i.e.
   identical to floating-point precision. Feeding ``snr`` in would let
   the network "recalibrate" by trivially inverting one input column
   rather than genuinely re-estimating error from signal *shape* -- a
   *target-leakage* pitfall (a feature that deterministically encodes
   the answer), distinct from the *degenerate*-feature pitfall caught
   in
   :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`
   (its original, always-1.0 ``gain`` feature carried no information
   at all), but both were caught the same way -- by checking a
   candidate feature's actual behaviour against real data before
   trusting it.

:func:`build_uncertainty_features_table` builds this table (plus the
``z_err_frac`` recalibration target) directly from a site collection.

Working directly with site collections
---------------------------------------
:meth:`UncertaintyCalibrator.apply` accepts and returns a site
collection, sites-in / sites-out like
:meth:`~pycsamt.ai.processing.qc.EMQCScorer.apply` -- but it *rescales*
``Z.z_err`` rather than replacing it outright, preserving whatever
per-component structure the field error already has (e.g. a larger
diagonal-term error) while shifting the overall level to match the
calibrated estimate.
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

__all__ = ["UncertaintyCalibrator", "build_uncertainty_features_table"]

_TABLE_COLS = ["snr", "swift_skew", "asym", "phase_xy", "phase_yx"]
# Model inputs deliberately exclude "snr": snr = amp / (err + eps) and
# z_err_frac = err / (amp + eps) share the exact same amp/err terms
# from _pooled_fractional_error, so snr is (to floating-point
# precision) the reciprocal of the training target, not an
# independent diagnostic -- verified on real L18PLT data:
# corr(1/snr, z_err_frac) = 1 - 2e-16. Feeding it in would let the
# network "recalibrate" by trivially inverting one input column
# rather than genuinely re-estimating error from signal shape; see
# the module docstring.
_FEATURE_COLS = ["swift_skew", "asym", "phase_xy", "phase_yx"]
_N_FEATURES = len(_FEATURE_COLS)


# ─────────────────────────────────────────────────────────────────────────────
# Feature-table builder (no torch / TF dependency)
# ─────────────────────────────────────────────────────────────────────────────


def build_uncertainty_features_table(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    """
    Build the per-(station, frequency) feature table
    :class:`UncertaintyCalibrator` expects.

    Parameters
    ----------
    sites : SiteCollection or compatible
    recursive, on_dup, strict, verbose
        Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

    Returns
    -------
    df : DataFrame
        Columns: ``station``, ``freq``, ``snr``, ``swift_skew``,
        ``asym``, ``phase_xy``, ``phase_yx`` (identical to
        :class:`~pycsamt.ai.processing.qc.EMQCScorer`'s feature table
        -- ``snr`` is reported for inspection but is *not* one of
        :class:`UncertaintyCalibrator`'s own model inputs, see the
        module docstring) plus ``z_err_frac`` -- the field-processing
        fractional error :math:`|z_\\text{err}|/|Z|`, pooled the same
        way as the ``snr`` column's own numerator/denominator. Rows
        whose site has no ``z_err`` at all are dropped (``z_err_frac``
        cannot be computed); empty (but correctly columned) when no
        station has error estimates.
    """
    try:
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )
    except ImportError as exc:
        raise ImportError(
            "emtools is required for build_uncertainty_features_table"
        ) from exc
    from .qc import _extract_qc_features

    cols = ["station", "freq", *_TABLE_COLS, "z_err_frac"]

    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    rows: list[dict[str, Any]] = []

    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        result = _get_z_block(ed, with_errors=True)
        if len(result) == 4:
            _, z, fr, ze = result
        else:
            _, z, fr = result[:3]
            ze = None
        if z is None or ze is None:
            continue

        F = _extract_qc_features(z, ze)
        frac = _pooled_fractional_error(z, ze)

        for fi, freq in enumerate(fr):
            row: dict[str, Any] = dict(station=st, freq=float(freq))
            row["snr"] = F[fi, 0]
            row["swift_skew"] = F[fi, 1]
            row["asym"] = F[fi, 2]
            row["phase_xy"] = F[fi, 3]
            row["phase_yx"] = F[fi, 4]
            row["z_err_frac"] = frac[fi]
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame.from_records(rows)[cols]


def _pooled_fractional_error(z: np.ndarray, ze: np.ndarray) -> np.ndarray:
    """
    Off-diagonal-pooled fractional error :math:`|z_\\text{err}|/|Z|`,
    one value per frequency -- the same Zxy/Zyx pooling
    :func:`~pycsamt.ai.processing.qc._extract_qc_features` uses for
    its ``snr`` column, so the two stay directly comparable
    (``snr`` :math:`\\approx` 1 / ``z_err_frac``).
    """
    zxy, zyx = z[:, 0, 1], z[:, 1, 0]
    zex, zey = ze[:, 0, 1], ze[:, 1, 0]
    amp = np.sqrt(0.5 * (np.abs(zxy) ** 2 + np.abs(zyx) ** 2))
    err = np.sqrt(0.5 * (np.abs(zex) ** 2 + np.abs(zey) ** 2))
    return err / (amp + 1e-24)


# ─────────────────────────────────────────────────────────────────────────────
# Network builders
# ─────────────────────────────────────────────────────────────────────────────


def _build_uncertainty_mlp_torch(
    n_features: int,
    hidden: tuple[int, ...],
    dropout: float,
) -> Any:
    """Regression MLP (single continuous output) -- PyTorch."""
    try:
        import torch.nn as nn
    except ImportError as exc:
        raise ImportError(
            "PyTorch is required for UncertaintyCalibrator"
        ) from exc

    dims = [n_features, *hidden]

    class _UncertaintyMLP(nn.Module):
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
            self.reg_head = nn.Linear(dims[-1], 1)

        def forward(self, x):
            return self.reg_head(self.backbone(x)).squeeze(-1)

    return _UncertaintyMLP()


def _build_uncertainty_mlp_tf(
    n_features: int,
    hidden: tuple[int, ...],
    dropout: float,
) -> Any:
    """Regression MLP (single continuous output) -- TensorFlow/Keras."""
    try:
        import tensorflow as tf
        from tensorflow.keras import Model, layers
    except ImportError as exc:
        raise ImportError(
            "TensorFlow is required for UncertaintyCalibrator "
            "(TF backend)"
        ) from exc

    inp = tf.keras.Input(shape=(n_features,), name="input")
    x = inp
    for h in hidden:
        x = layers.Dense(h)(x)
        x = layers.BatchNormalization()(x)
        x = layers.ReLU()(x)
        x = layers.Dropout(dropout)(x)
    out = layers.Dense(1, name="reg")(x)

    return Model(inp, out, name="uncertainty_mlp")


# ─────────────────────────────────────────────────────────────────────────────
# UncertaintyCalibrator
# ─────────────────────────────────────────────────────────────────────────────


class UncertaintyCalibrator(BaseEMProcessor):
    """
    Learned per-(station, frequency) uncertainty / error-floor estimator.

    Parameters
    ----------
    n_features : int, default 4
        Length of the per-cell feature vector. :meth:`fit`/
        :meth:`transform` on a
        :func:`build_uncertainty_features_table` DataFrame always use
        the 4-column ``[swift_skew, asym, phase_xy, phase_yx]``
        layout -- deliberately excluding ``snr`` (see the module
        docstring); raise this only when fitting on a custom feature
        matrix with more columns of your own.
    hidden : tuple of int, default (64, 32)
        Hidden-layer widths of the regressor backbone.
    dropout : float, default 0.1
        Dropout probability in each hidden layer.
    lr : float, default 1e-3
        Default learning rate; can be overridden in :meth:`fit`.
    device : str or None
        ``"cpu"``, ``"cuda"``, ``"mps"``, ``"/GPU:0"``, or ``None``
        for auto-detect via the active backend.

    Notes
    -----
    Falls back to a random-forest regressor (scikit-learn) when
    neither PyTorch nor TensorFlow is available.

    Examples
    --------
    >>> from pycsamt.ai.processing.uncertainty import (
    ...     UncertaintyCalibrator, build_uncertainty_features_table,
    ... )
    >>> feats = build_uncertainty_features_table(sites)  # doctest: +SKIP
    >>> cal = UncertaintyCalibrator.from_features_table(
    ...     feats, epochs=150,
    ... )  # doctest: +SKIP
    >>> table = cal.predict_table(sites)  # doctest: +SKIP
    """

    def __init__(
        self,
        n_features: int = 4,
        *,
        hidden: tuple[int, ...] = (64, 32),
        dropout: float = 0.1,
        lr: float = 1e-3,
        device: str | None = None,
    ) -> None:
        self.n_features = int(n_features)
        self.hidden = tuple(hidden)
        self.dropout = float(dropout)
        self.lr = float(lr)
        self.device = device

        self._network: Any = None
        self._rf: Any = None
        self._use_rf: bool = False
        self._backend_name: str | None = None
        self._x_mean: np.ndarray | None = None
        self._x_std: np.ndarray | None = None
        self._y_mean: float = 0.0
        self._y_std: float = 1.0
        self._log_target: bool = True
        self._is_fitted: bool = False
        self._history: dict[str, list] = {}

    # ─── factory from emtools ─────────────────────────────────────────────

    @classmethod
    def from_features_table(
        cls,
        df: pd.DataFrame,
        *,
        log_target: bool = True,
        **fit_kwargs,
    ) -> UncertaintyCalibrator:
        """
        Construct and train a calibrator from a
        :func:`build_uncertainty_features_table` DataFrame.

        Parameters
        ----------
        df : DataFrame
            Must contain a ``z_err_frac`` column (self-training
            recalibration target); pass ``y=`` directly to
            :meth:`fit` instead for an empirical/bootstrap target.
        log_target : bool, default True
            Train in log10 space (recommended -- fractional error is
            right-skewed); :meth:`transform` always returns the
            linear-space value regardless.
        **fit_kwargs
            Passed to :meth:`fit` (e.g. ``epochs=150``).

        Returns
        -------
        UncertaintyCalibrator
        """
        obj = cls()
        obj.fit(df, log_target=log_target, **fit_kwargs)
        return obj

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(
        self,
        X: np.ndarray | pd.DataFrame,
        y: np.ndarray | None = None,
        *,
        log_target: bool = True,
        epochs: int = 150,
        batch_size: int = 64,
        lr: float | None = None,
        val_frac: float = 0.15,
        seed: int | None = None,
        verbose: bool = True,
    ) -> UncertaintyCalibrator:
        """
        Train the uncertainty calibrator.

        Parameters
        ----------
        X : ndarray (n_samples, n_features) or DataFrame
            A DataFrame from :func:`build_uncertainty_features_table`
            self-trains against its own ``z_err_frac`` column when
            ``y`` is not given.
        y : ndarray (n_samples,) or None
            Target error/uncertainty values, any positive-valued
            scale (recalibration fractional error, an empirical
            bootstrap spread, ...). Required when ``X`` is a plain
            ndarray -- there is no rule to synthesise a numeric
            target from features alone the way the classifiers in
            this package synthesise class labels.
        log_target : bool, default True
            Train in log10 space. :meth:`transform` output is always
            de-logged back to the original target's linear scale.
        epochs, batch_size, lr, val_frac, seed, verbose
            Training hyper-parameters.

        Returns
        -------
        self

        Raises
        ------
        ValueError
            ``y`` is ``None`` and no target can be recovered from
            ``X`` (see above), or no finite training rows remain.
        """
        X_arr, y_arr = self._coerce_Xy(X, y)
        if X_arr.shape[1] != self.n_features:
            raise ValueError(
                f"Feature width mismatch: got {X_arr.shape[1]} "
                f"columns, expected n_features={self.n_features}."
            )

        finite_y = np.isfinite(y_arr) & (y_arr > 0)
        valid = np.all(np.isfinite(X_arr), axis=1) & finite_y
        X_arr, y_arr = X_arr[valid], y_arr[valid]
        if len(X_arr) == 0:
            raise ValueError(
                "No valid (finite, positive) training rows remain "
                "after filtering."
            )

        self._log_target = bool(log_target)
        y_fit = np.log10(y_arr) if self._log_target else y_arr

        self._x_mean = X_arr.mean(axis=0, keepdims=True)
        self._x_std = X_arr.std(axis=0, keepdims=True) + 1e-8
        Xn = (X_arr - self._x_mean) / self._x_std

        self._y_mean = float(y_fit.mean())
        self._y_std = float(y_fit.std() + 1e-8)
        yn = (y_fit - self._y_mean) / self._y_std
        _lr = float(lr) if lr is not None else self.lr

        try:
            self._backend_name = active_backend()
            if self._backend_name == "tensorflow":
                self._fit_tensorflow(
                    Xn,
                    yn,
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
                    yn,
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
            self._fit_rf(Xn, yn, verbose=verbose)
            self._use_rf = True

        self._is_fitted = True
        return self

    def transform(self, X: np.ndarray | pd.DataFrame) -> np.ndarray:
        """
        Predict calibrated per-cell uncertainty.

        Parameters
        ----------
        X : ndarray (n_samples, n_features) or DataFrame

        Returns
        -------
        y_pred : ndarray, shape (n_samples,)
            Calibrated error/uncertainty estimate, de-logged back to
            the same linear scale as the :meth:`fit` target (e.g.
            fractional error, in ``[0, 1]``-ish units, for the
            default recalibration target).
        """
        if not self._is_fitted:
            raise RuntimeError("Call fit() before transform().")
        X_arr = self._coerce_X(X)
        Xn = (X_arr - self._x_mean) / self._x_std
        yn = self._predict(Xn)
        y_fit = yn * self._y_std + self._y_mean
        return 10.0**y_fit if self._log_target else y_fit

    def predict_table(
        self,
        sites: Any,
        *,
        recursive: bool = True,
        on_dup: str = "replace",
        strict: bool = False,
        verbose: int = 0,
    ) -> pd.DataFrame:
        """
        Calibrate uncertainty for an entire site collection.

        Parameters
        ----------
        sites : SiteCollection or compatible
        recursive, on_dup, strict, verbose
            Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

        Returns
        -------
        df : DataFrame
            :func:`build_uncertainty_features_table`'s columns plus
            ``z_err_frac_calibrated``.
        """
        df = build_uncertainty_features_table(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
        if df.empty:
            return df

        out = df.copy()
        out["z_err_frac_calibrated"] = self.transform(df)
        return out

    # ─── sites-in / sites-out ──────────────────────────────────────────────

    def apply(
        self,
        sites: Any,
        *,
        inplace: bool = False,
        floor_frac: float | None = None,
        recursive: bool = True,
        on_dup: str = "replace",
        strict: bool = False,
        verbose: int = 0,
    ) -> Any:
        """
        Rescale ``Z.z_err`` to the calibrated level, sites-in /
        sites-out.

        Unlike :meth:`~pycsamt.ai.processing.qc.EMQCScorer.apply`,
        this does not replace ``z_err`` outright: it computes the
        ratio between the calibrated fractional error and the
        station's own original pooled fractional error at each
        frequency, then scales *all four* tensor components' existing
        ``z_err`` by that one ratio -- preserving whatever
        per-component structure the field error already carries
        (e.g. a larger diagonal-term error) while shifting the
        overall level to the calibrated estimate.

        Parameters
        ----------
        sites : SiteCollection or compatible
        inplace : bool, default ``False``
            When ``False`` (the default), *sites* is left untouched
            and a corrected copy is returned. When ``True``, *sites*
            is modified in place and returned.
        floor_frac : float or None, default None
            Minimum fractional error to enforce after calibration
            (e.g. ``0.02`` for a 2% floor) -- applied with
            ``np.maximum`` after the network prediction, before
            rescaling. ``None`` applies no floor.
        recursive, on_dup, strict, verbose
            Passed to :func:`~pycsamt.emtools._core.ensure_sites`.

        Returns
        -------
        corrected : Sites
            A site collection with the same stations; ``Z.z_err`` is
            rescaled station-by-station, frequency-by-frequency.
            Stations with no usable ``z_err`` to begin with are left
            untouched (there is nothing to rescale from).

        Examples
        --------
        >>> cal = UncertaintyCalibrator().fit(feats)  # doctest: +SKIP
        >>> calibrated = cal.apply(sites)  # doctest: +SKIP
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
                "emtools is required for UncertaintyCalibrator.apply()"
            ) from exc
        from .qc import _extract_qc_features

        S = ensure_sites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )

        def _one(Si):
            ed = next(_iter_items(Si))
            result = _get_z_block(ed, with_errors=True)
            if len(result) != 4:
                return Si
            Z, z, _fr, ze = result
            if z is None or ze is None:
                return Si

            # _extract_qc_features columns: snr, swift_skew, asym,
            # phase_xy, phase_yx -- drop snr (column 0), see
            # _FEATURE_COLS.
            F = _extract_qc_features(z, ze)[:, 1:]
            Xn = (F - self._x_mean) / self._x_std
            yn = self._predict(Xn)
            y_fit = yn * self._y_std + self._y_mean
            frac_pred = 10.0**y_fit if self._log_target else y_fit
            if floor_frac is not None:
                frac_pred = np.maximum(frac_pred, floor_frac)

            frac_orig = _pooled_fractional_error(z, ze)
            scale = frac_pred / np.maximum(frac_orig, 1e-12)
            Z.z_err = ze * scale[:, None, None]
            return Si

        return _apply_each(S, _one, inplace=inplace, verbose=verbose)

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
            torch.from_numpy(y[ti].astype(np.float32)),
        )
        mse = nn.MSELoss()

        self._network = _build_uncertainty_mlp_torch(
            self.n_features, self.hidden, self.dropout
        ).to(dev)
        opt = torch.optim.Adam(self._network.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ReduceLROnPlateau(
            opt, factor=0.5, patience=8, min_lr=1e-6
        )

        Xva = torch.from_numpy(Xn[vi].astype(np.float32)).to(dev)
        yva = torch.from_numpy(y[vi].astype(np.float32)).to(dev)

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
                loss = mse(out, yb)
                opt.zero_grad()
                loss.backward()
                opt.step()
                ep_loss += loss.item() * len(xb)
            ep_loss /= len(ti)

            self._network.eval()
            with torch.no_grad():
                v_out = self._network(Xva)
                v_loss = mse(v_out, yva).item()

            sched.step(v_loss)
            train_losses.append(ep_loss)
            val_losses.append(v_loss)

            if v_loss < best_val:
                best_val = v_loss
                best_state = copy.deepcopy(self._network.state_dict())

            if verbose and (ep % max(1, epochs // 5) == 0 or ep == 1):
                print(
                    f"  UncertaintyCalibrator  ep {ep:>4d}/{epochs}  "
                    f"loss={ep_loss:.4f}  val_loss={v_loss:.4f}"
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
        ytr = y[ti].astype(np.float32)
        Xva = Xn[vi].astype(np.float32)
        yva = y[vi].astype(np.float32)

        dev = resolve_device(self.device)
        with tf.device(dev):
            self._network = _build_uncertainty_mlp_tf(
                self.n_features, self.hidden, self.dropout
            )
            self._network.compile(
                optimizer=tf.keras.optimizers.Adam(learning_rate=lr),
                loss=tf.keras.losses.MeanSquaredError(),
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
            from sklearn.ensemble import RandomForestRegressor
        except ImportError as exc:
            raise ImportError(
                "PyTorch, TensorFlow, or scikit-learn is required for "
                "UncertaintyCalibrator"
            ) from exc

        self._rf = RandomForestRegressor(n_estimators=200, random_state=0)
        valid = np.all(np.isfinite(Xn), axis=1) & np.isfinite(y)
        self._rf.fit(Xn[valid], y[valid])
        if verbose:
            print("  UncertaintyCalibrator (RandomForest fallback) fitted.")

    def _predict(self, Xn: np.ndarray) -> np.ndarray:
        if self._use_rf:
            valid = np.all(np.isfinite(Xn), axis=1)
            out = np.zeros(len(Xn))
            if valid.any():
                out[valid] = self._rf.predict(Xn[valid])
            return out

        if self._backend_name == "tensorflow":
            pred = self._network.predict(Xn.astype(np.float32), verbose=0)
            return np.asarray(pred).reshape(-1)

        import torch

        dev = next(self._network.parameters()).device
        self._network.eval()
        with torch.no_grad():
            t = torch.from_numpy(Xn.astype(np.float32)).to(dev)
            out = self._network(t).cpu().numpy()
        return out

    def _coerce_X(self, X) -> np.ndarray:
        if isinstance(X, pd.DataFrame):
            missing = [c for c in _FEATURE_COLS if c not in X.columns]
            if missing:
                raise ValueError(
                    f"X is missing required feature columns: {missing}."
                )
            return X[_FEATURE_COLS].to_numpy(dtype=float)
        return np.asarray(X, dtype=np.float32)

    def _coerce_Xy(self, X, y):
        if isinstance(X, pd.DataFrame):
            X_arr = self._coerce_X(X)
            if y is not None:
                y_arr = np.asarray(y, dtype=float)
            elif "z_err_frac" in X.columns:
                y_arr = X["z_err_frac"].to_numpy(dtype=float)
            else:
                raise ValueError(
                    "fit(X, y=None) needs a 'z_err_frac' column "
                    "(from build_uncertainty_features_table) to "
                    "self-train a recalibration target; pass y= "
                    "explicitly for an empirical/bootstrap target "
                    "instead."
                )
        else:
            X_arr = np.asarray(X, dtype=np.float32)
            if y is None:
                raise ValueError(
                    "fit(X, y=None) requires either a features "
                    "DataFrame with a 'z_err_frac' column (see "
                    "build_uncertainty_features_table) or an "
                    "explicit y= target array."
                )
            y_arr = np.asarray(y, dtype=float)
        return X_arr, y_arr

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "n_features": self.n_features,
            "hidden": list(self.hidden),
            "dropout": self.dropout,
            "lr": self.lr,
            "device": self.device,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {}
        if self._x_mean is not None:
            out["_x_mean"] = self._x_mean
            out["_x_std"] = self._x_std
        out["_y_mean"] = np.array([self._y_mean], dtype=np.float64)
        out["_y_std"] = np.array([self._y_std], dtype=np.float64)
        out["_log_target"] = np.array([int(self._log_target)])
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
        y_mean_blob = weights.pop("_y_mean", None)
        self._y_mean = (
            float(y_mean_blob[0]) if y_mean_blob is not None else 0.0
        )
        y_std_blob = weights.pop("_y_std", None)
        self._y_std = (
            float(y_std_blob[0]) if y_std_blob is not None else 1.0
        )
        log_target_blob = weights.pop("_log_target", None)
        self._log_target = (
            bool(int(log_target_blob[0]))
            if log_target_blob is not None
            else True
        )
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
                self._network = _build_uncertainty_mlp_tf(
                    self.n_features, self.hidden, self.dropout
                )
            else:
                self._network = _build_uncertainty_mlp_torch(
                    self.n_features, self.hidden, self.dropout
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
            ``{"train_loss": [...], "val_loss": [...]}``, one value
            per epoch, in normalised log-target loss units. Empty
            when the random-forest fallback was used (no epoch loop)
            or before :meth:`fit` has been called. Pass directly to
            :func:`~pycsamt.ai.processing.plot.plot_training_history`.
        """
        return dict(self._history)

    def __repr__(self) -> str:
        backend = self._backend_name or ("rf" if self._use_rf else "torch")
        status = "fitted" if self._is_fitted else "unfitted"
        return (
            f"UncertaintyCalibrator(n_features={self.n_features}, "
            f"{backend}, {status})"
        )

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
EMQCScorer — ML-based per-frequency quality control scorer.

Extends the rule-based QC tools in :mod:`pycsamt.emtools.qc` with an
Isolation Forest anomaly model trained on signal-quality features.

Features extracted per (site, frequency) observation
-----------------------------------------------------
* SNR: :math:`|\\bar{Z}| / \\sigma(Z)` (signal-to-noise ratio)
* Phase stability: coefficient of variation of the off-diagonal phases
* Swift skew: :math:`|\\beta_\\text{Swift}| = |(Z_{xx}+Z_{yy})/(Z_{xy}-Z_{yx})|`
* Phase tensor skew :math:`|\\beta|` (requires full tensor)
* Off-diagonal amplitude asymmetry:
  :math:`\\log_{10}(|Z_{xy}|/|Z_{yx}|)`

The combined score lies in ``[0, 1]`` — 1 means good quality, 0 means
flagged bad.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .._base import BaseEMProcessor

__all__ = ["EMQCScorer"]


# ─────────────────────────────────────────────────────────────────────────────
# Internal feature extraction
# ─────────────────────────────────────────────────────────────────────────────


def _extract_qc_features(z: np.ndarray, ze: np.ndarray | None) -> np.ndarray:
    """
    Extract a (n_freqs, 5) feature matrix from a single site's Z data.

    Parameters
    ----------
    z : ndarray, shape (n_freqs, 2, 2), complex
    ze : ndarray or None, shape (n_freqs, 2, 2), real errors

    Returns
    -------
    F : ndarray, shape (n_freqs, 5)
        Columns: SNR, |β_swift|, |Zxy/Zyx|, phase_xy, phase_yx
    """
    n = z.shape[0]
    F = np.full((n, 5), np.nan, dtype=float)

    zxy = z[:, 0, 1]
    zyx = z[:, 1, 0]
    zxx = z[:, 0, 0]
    zyy = z[:, 1, 1]

    # SNR from error array
    if ze is not None and ze.shape == z.shape:
        amp = np.sqrt(0.5 * (np.abs(zxy) ** 2 + np.abs(zyx) ** 2))
        err = np.sqrt(0.5 * (np.abs(ze[:, 0, 1]) ** 2 + np.abs(ze[:, 1, 0]) ** 2))
        F[:, 0] = amp / (err + 1e-24)
    else:
        # Estimate SNR from local spectral smoothness
        amp_xy = np.abs(zxy)
        if len(amp_xy) > 3:
            from scipy.signal import medfilt

            smooth = medfilt(amp_xy, kernel_size=3)
            residual = np.abs(amp_xy - smooth)
            F[:, 0] = amp_xy / (residual + 1e-24)
        else:
            F[:, 0] = np.nan

    # Swift skew  |β| = |(Zxx+Zyy) / (Zxy-Zyx)|
    denom = np.abs(zxy - zyx)
    swift = np.abs(zxx + zyy) / (denom + 1e-24)
    F[:, 1] = swift

    # Off-diagonal amplitude asymmetry
    F[:, 2] = np.log10(np.maximum(np.abs(zxy), 1e-24) / np.maximum(np.abs(zyx), 1e-24))

    # Phase (degrees)
    F[:, 3] = np.degrees(np.angle(zxy))
    F[:, 4] = np.degrees(np.angle(zyx))

    return F


def _sites_to_feature_df(sites: Any) -> pd.DataFrame:
    """
    Convert a site collection to a DataFrame of per-(site, freq) features.
    """
    try:
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
            _name,
            ensure_sites,
        )
    except ImportError as exc:
        raise ImportError("emtools is required for site-based QC scoring") from exc

    S = ensure_sites(sites, recursive=True, on_dup="replace")
    rows: list[dict[str, Any]] = []

    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        result = _get_z_block(ed, with_errors=True)
        if len(result) == 4:
            _, z, fr, ze = result
        else:
            _, z, fr = result[:3]
            ze = None

        if z is None:
            continue

        F = _extract_qc_features(z, ze)
        for fi, freq in enumerate(fr):
            row: dict[str, Any] = dict(station=st, freq=float(freq))
            row["snr"] = F[fi, 0]
            row["swift_skew"] = F[fi, 1]
            row["asym"] = F[fi, 2]
            row["phase_xy"] = F[fi, 3]
            row["phase_yx"] = F[fi, 4]
            rows.append(row)

    return pd.DataFrame.from_records(rows)


# ─────────────────────────────────────────────────────────────────────────────
# EMQCScorer
# ─────────────────────────────────────────────────────────────────────────────


class EMQCScorer(BaseEMProcessor):
    """
    ML-based per-frequency QC scorer for MT impedance data.

    Combines hard thresholds on SNR and Swift skew with an
    ``IsolationForest`` anomaly model fitted on extracted signal-quality
    features.  The final quality score for each (site, frequency) cell
    is :math:`s \\in [0, 1]`; observations below ``score_threshold``
    should be rejected before inversion.

    Parameters
    ----------
    contamination : float, default 0.05
        Expected fraction of contaminated samples, passed directly to
        :class:`sklearn.ensemble.IsolationForest`.
    snr_threshold : float, default 3.0
        Observations with SNR below this value are hard-flagged as bad
        regardless of the ML score.
    skew_threshold : float, default 0.3
        Observations with Swift skew above this value are hard-flagged.
    score_threshold : float, default 0.5
        Quality scores below this value are considered bad.
    use_ml : bool, default True
        When ``False``, only the rule-based thresholds are applied and
        the IsolationForest is not fitted.
    n_estimators : int, default 100
        Number of trees in the IsolationForest.
    random_state : int or None

    Examples
    --------
    >>> from pycsamt.ai.processing import EMQCScorer
    >>> scorer = EMQCScorer(snr_threshold=5.0, use_ml=False)
    >>> scorer.fit(sites)  # doctest: +SKIP
    EMQCScorer(rule_only)
    >>> tbl = scorer.score_table(sites)  # doctest: +SKIP
    """

    def __init__(
        self,
        contamination: float = 0.05,
        snr_threshold: float = 3.0,
        skew_threshold: float = 0.3,
        score_threshold: float = 0.5,
        use_ml: bool = True,
        n_estimators: int = 100,
        random_state: int | None = None,
    ) -> None:
        self.contamination = float(contamination)
        self.snr_threshold = float(snr_threshold)
        self.skew_threshold = float(skew_threshold)
        self.score_threshold = float(score_threshold)
        self.use_ml = bool(use_ml)
        self.n_estimators = int(n_estimators)
        self.random_state = random_state

        self._model: Any = None  # IsolationForest
        self._feat_cols: list[str] = [
            "snr",
            "swift_skew",
            "asym",
            "phase_xy",
            "phase_yx",
        ]
        self._is_fitted: bool = False

    # ─── BaseEMProcessor interface ────────────────────────────────────────

    def fit(self, X: Any, **kwargs) -> EMQCScorer:
        """
        Fit the QC model on a training set.

        Parameters
        ----------
        X : SiteCollection, dict, list of Z, or ndarray (n_samples, 5)
            Training data.  Site collections are converted automatically.
            Pass a precomputed feature matrix to skip extraction.

        Returns
        -------
        self
        """
        if not self.use_ml:
            self._is_fitted = True
            return self

        try:
            from sklearn.ensemble import IsolationForest
        except ImportError as exc:
            raise ImportError(
                "scikit-learn is required for EMQCScorer with use_ml=True"
            ) from exc

        feat = self._to_feature_matrix(X)
        valid = np.all(np.isfinite(feat), axis=1)
        feat_clean = feat[valid]
        if len(feat_clean) == 0:
            raise ValueError("No valid (finite) feature rows found in training data.")

        self._model = IsolationForest(
            n_estimators=self.n_estimators,
            contamination=self.contamination,
            random_state=self.random_state,
        )
        self._model.fit(feat_clean)
        self._is_fitted = True
        return self

    def transform(self, X: Any) -> np.ndarray:
        """
        Compute quality scores.

        Parameters
        ----------
        X : SiteCollection or ndarray (n_samples, 5)

        Returns
        -------
        scores : ndarray, shape (n_samples,)
            Quality scores in ``[0, 1]``; 1 = good, 0 = bad.
        """
        feat = self._to_feature_matrix(X)
        return self._score_features(feat)

    def score_table(self, sites: Any) -> pd.DataFrame:
        """
        Return per-(site, frequency) QC table.

        Parameters
        ----------
        sites : SiteCollection or compatible

        Returns
        -------
        df : DataFrame
            Columns: station, freq, snr, swift_skew, asym,
            phase_xy, phase_yx, score, flag (0=bad, 1=good).
        """
        df = _sites_to_feature_df(sites)
        if df.empty:
            return df

        feat = df[self._feat_cols].to_numpy(dtype=float)
        scores = self._score_features(feat)
        df = df.copy()
        df["score"] = scores
        df["flag"] = (scores >= self.score_threshold).astype(int)
        return df

    # ─── internal ─────────────────────────────────────────────────────────

    def _to_feature_matrix(self, X: Any) -> np.ndarray:
        if isinstance(X, np.ndarray):
            return X.astype(float)
        if isinstance(X, pd.DataFrame):
            return X[self._feat_cols].to_numpy(dtype=float)
        # assume site collection
        df = _sites_to_feature_df(X)
        return df[self._feat_cols].to_numpy(dtype=float)

    def _score_features(self, feat: np.ndarray) -> np.ndarray:
        n = len(feat)
        scores = np.ones(n, dtype=float)

        # Hard rules
        snr_col = feat[:, 0]
        skew_col = feat[:, 1]
        bad_snr = np.isfinite(snr_col) & (snr_col < self.snr_threshold)
        bad_skew = np.isfinite(skew_col) & (skew_col > self.skew_threshold)
        scores[bad_snr | bad_skew] = 0.0

        # ML scores (IsolationForest returns -1=anomaly, +1=normal)
        if self.use_ml and self._model is not None:
            valid = np.all(np.isfinite(feat), axis=1)
            if valid.any():
                raw = self._model.decision_function(feat[valid])
                # Normalise to [0, 1]: higher decision_function = more normal
                rmin, rmax = raw.min(), raw.max()
                if rmax > rmin:
                    ml_score = (raw - rmin) / (rmax - rmin)
                else:
                    ml_score = np.ones_like(raw)
                # Combine: geometric mean of rule-based (binary) and ML score
                full_ml = np.full(n, 1.0)
                full_ml[valid] = ml_score
                scores = np.sqrt(scores * full_ml)

        return scores

    # ─── serialisation ────────────────────────────────────────────────────

    def _get_params(self) -> dict[str, Any]:
        return {
            "contamination": self.contamination,
            "snr_threshold": self.snr_threshold,
            "skew_threshold": self.skew_threshold,
            "score_threshold": self.score_threshold,
            "use_ml": self.use_ml,
            "n_estimators": self.n_estimators,
            "random_state": self.random_state,
        }

    def _get_weights(self) -> dict[str, np.ndarray]:
        if self._model is None:
            return {}
        try:
            import io
            import pickle

            buf = io.BytesIO()
            pickle.dump(self._model, buf)
            return {"_iso_model": np.frombuffer(buf.getvalue(), dtype=np.uint8)}
        except Exception:
            return {}

    def _load_weights(self, weights: dict[str, np.ndarray]) -> None:
        if "_iso_model" in weights:
            try:
                import io
                import pickle

                buf = io.BytesIO(bytes(weights["_iso_model"]))
                self._model = pickle.load(buf)
                self._is_fitted = True
            except Exception:
                pass

    def __repr__(self) -> str:
        mode = "ml+rules" if self.use_ml else "rule_only"
        return f"EMQCScorer({mode})"

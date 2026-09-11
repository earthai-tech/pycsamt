# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Focused coverage-gap tests for pycsamt.ai.processing.

Complements test_processing.py and test_ai_processing_api_contracts.py by
exercising branches those files don't reach: the ML/torch training paths
(this environment has PyTorch but not TensorFlow installed), sites-based
feature extraction (via the bundled ``data/3edis`` EDI files), verbose
logging branches, and save/load round trips.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_EDI_DIR = _PROJECT_ROOT / "data" / "3edis"
_HAS_EDIS = _EDI_DIR.exists() and any(_EDI_DIR.glob("*.edi"))


@pytest.fixture(scope="module")
def sites():
    if not _HAS_EDIS:
        pytest.skip(f"3edis dataset not found: {_EDI_DIR}")
    from pycsamt.agents import MTLoaderAgent

    r = MTLoaderAgent().execute({"path": str(_EDI_DIR)})
    if r.status != "success":
        pytest.skip("Could not load 3edis sites.")
    return r["sites"]


# ─────────────────────────────────────────────────────────────────────────────
# qc.py
# ─────────────────────────────────────────────────────────────────────────────


def test_extract_qc_features_without_errors_uses_smoothness_snr():
    from pycsamt.ai.processing.qc import _extract_qc_features

    rng = np.random.default_rng(0)
    n = 12
    zxy = (1.0 + 0.1 * rng.standard_normal(n)) + 1j * (
        1.0 + 0.1 * rng.standard_normal(n)
    )
    zyx = -zxy + 0.05 * rng.standard_normal(n)
    z = np.zeros((n, 2, 2), dtype=complex)
    z[:, 0, 1] = zxy
    z[:, 1, 0] = zyx
    z[:, 0, 0] = 0.02
    z[:, 1, 1] = 0.02

    F = _extract_qc_features(z, None)

    assert F.shape == (n, 5)
    assert np.all(np.isfinite(F[:, 0]))
    assert np.all(F[:, 0] > 0)


def test_sites_to_feature_df_and_prepare_z_features(sites):
    from pycsamt.ai.processing.denoise import prepare_z_features
    from pycsamt.ai.processing.qc import _sites_to_feature_df

    df = _sites_to_feature_df(sites)
    assert set(
        [
            "station",
            "freq",
            "snr",
            "swift_skew",
            "asym",
            "phase_xy",
            "phase_yx",
        ]
    ).issubset(df.columns)
    assert len(df) > 0

    X4 = prepare_z_features(sites, n_components=4)
    assert X4.ndim == 3 and X4.shape[1] == 4
    assert X4.dtype == np.float32

    X8 = prepare_z_features(sites, n_components=8)
    assert X8.shape[1] == 8


def test_emqc_scorer_fit_transform_with_ml():
    from pycsamt.ai.processing.qc import EMQCScorer

    rng = np.random.default_rng(1)
    feat = np.column_stack(
        [
            rng.uniform(2.0, 20.0, 60),  # snr
            rng.uniform(0.0, 0.5, 60),  # swift_skew
            rng.uniform(-0.5, 0.5, 60),  # asym
            rng.uniform(20.0, 50.0, 60),  # phase_xy
            rng.uniform(-160.0, -120.0, 60),  # phase_yx
        ]
    )

    scorer = EMQCScorer(use_ml=True, n_estimators=10, random_state=0)
    assert scorer.fit(feat) is scorer
    assert scorer._model is not None

    scores = scorer.transform(feat)
    assert scores.shape == (60,)
    assert np.all((scores >= 0.0) & (scores <= 1.0))
    assert "ml+rules" in repr(scorer)


def test_emqc_scorer_score_table_and_site_transform(sites):
    from pycsamt.ai.processing.qc import EMQCScorer

    scorer = EMQCScorer(use_ml=False)
    scorer.fit(sites)

    tbl = scorer.score_table(sites)
    assert {"station", "freq", "score", "flag"}.issubset(tbl.columns)
    assert len(tbl) > 0
    assert set(np.unique(tbl["flag"])).issubset({0, 1})

    # transform() with a raw site collection (not ndarray/DataFrame) exercises
    # the `_to_feature_matrix` site-collection fallback branch.
    scores = scorer.transform(sites)
    assert scores.shape[0] == len(tbl)


def test_emqc_scorer_save_load_round_trip_with_ml_model(tmp_path):
    from pycsamt.ai.processing.qc import EMQCScorer

    rng = np.random.default_rng(2)
    feat = np.column_stack(
        [
            rng.uniform(2.0, 20.0, 40),
            rng.uniform(0.0, 0.5, 40),
            rng.uniform(-0.5, 0.5, 40),
            rng.uniform(20.0, 50.0, 40),
            rng.uniform(-160.0, -120.0, 40),
        ]
    )
    scorer = EMQCScorer(use_ml=True, n_estimators=10, random_state=0)
    scorer.fit(feat)
    scores = scorer.transform(feat)

    path = tmp_path / "qc_ml.npz"
    scorer.save(path)
    loaded = EMQCScorer.load(path)

    np.testing.assert_allclose(loaded.transform(feat), scores)


# ─────────────────────────────────────────────────────────────────────────────
# denoise.py
# ─────────────────────────────────────────────────────────────────────────────


def _denoise_X(n=24, n_comp=4, n_freqs=16, seed=5):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, n_comp, n_freqs)).astype(np.float32)


def test_emdenoiser_verbose_training_prints(capsys):
    from pycsamt.ai.processing.denoise import EMDenoiser

    X = _denoise_X(n=16, n_freqs=12)
    den = EMDenoiser(channels=(4, 8, 4))
    den.fit(X, epochs=10, verbose=True)
    out = capsys.readouterr().out
    assert "Epoch" in out


def test_emdenoiser_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.denoise import EMDenoiser

    X = _denoise_X(n=16, n_freqs=12)
    den = EMDenoiser(channels=(4, 8, 4))
    den.fit(X, epochs=2, verbose=False)
    assert den._backend_name == "torch"
    out_before = den.transform(X)

    path = tmp_path / "denoiser.npz"
    den.save(path)
    loaded = EMDenoiser.load(path)

    out_after = loaded.transform(X)
    np.testing.assert_allclose(out_before, out_after, rtol=1e-4, atol=1e-5)


def test_emdenoiser_numpy_fallback_scipy_missing(monkeypatch):
    import pycsamt.ai.processing.denoise as denoise_mod

    monkeypatch.setattr(denoise_mod, "active_backend", lambda: "none")
    monkeypatch.setitem(sys.modules, "scipy.ndimage", None)

    X = _denoise_X(n=8, n_freqs=10)
    den = denoise_mod.EMDenoiser()
    den.fit(X, epochs=1, verbose=False)
    assert den._use_numpy is True

    out = den.transform(X)
    assert out.shape == X.shape
    assert np.all(np.isfinite(out))


# ─────────────────────────────────────────────────────────────────────────────
# anomaly.py
# ─────────────────────────────────────────────────────────────────────────────


def _anomaly_X(n=40, n_feat=20, seed=6):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, n_feat)).astype(np.float32)


def test_anomaly_detector_verbose_training_prints(capsys):
    from pycsamt.ai.processing.anomaly import AnomalyDetector

    X = _anomaly_X()
    det = AnomalyDetector(latent_dim=4, channels=(8,))
    det.fit(X, epochs=10, verbose=True)
    out = capsys.readouterr().out
    assert "AnomalyDetector" in out


def test_anomaly_detector_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.anomaly import AnomalyDetector

    X = _anomaly_X()
    det = AnomalyDetector(latent_dim=4, channels=(8,))
    det.fit(X, epochs=2, verbose=False)
    assert det._backend_name == "torch"
    scores_before = det.transform(X)

    path = tmp_path / "anomaly_torch.npz"
    det.save(path)
    loaded = AnomalyDetector.load(path)

    scores_after = loaded.transform(X)
    np.testing.assert_allclose(scores_before, scores_after, rtol=1e-4, atol=1e-5)
    np.testing.assert_array_equal(loaded.flag_anomalies(X), det.flag_anomalies(X))


def test_anomaly_detector_pca_fallback_verbose_prints(monkeypatch, capsys):
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.anomaly as anomaly_mod

    monkeypatch.setattr(
        anomaly_mod.AnomalyDetector,
        "_fit_torch",
        lambda self, Xn, **kwargs: (_ for _ in ()).throw(ImportError("no torch")),
    )
    X = _anomaly_X(n=20, n_feat=10)
    det = anomaly_mod.AnomalyDetector(latent_dim=3)
    det.fit(X, epochs=1, verbose=True)
    out = capsys.readouterr().out
    assert "PCA fallback" in out
    assert det._use_pca is True


# ─────────────────────────────────────────────────────────────────────────────
# classify.py
# ─────────────────────────────────────────────────────────────────────────────


def _dim_Xy(n=80, seed=9):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, 5)).astype(np.float32)
    X[:, 0] = np.abs(X[:, 0]) * 4
    X[:, 1] = np.abs(X[:, 1]) * 0.3
    from pycsamt.ai.processing.classify import _rule_labels

    y = _rule_labels(X[:, 0], X[:, 1])
    strike = np.where(y == 1, rng.uniform(-90.0, 90.0, n), np.nan).astype(np.float32)
    return X, y, strike


def test_dim_classifier_strike_training_torch_verbose(capsys):
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, strike = _dim_Xy()
    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0)
    clf.fit(X, y, strike=strike, epochs=10, verbose=True)
    out = capsys.readouterr().out
    assert "DimClassifier" in out
    assert clf._backend_name == "torch"

    pred_strike = clf.predict_strike(X)
    assert pred_strike.shape == (len(X),)
    labels = clf.predict(X)
    assert np.all(np.isnan(pred_strike[labels != 1]))


def test_dim_classifier_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, strike = _dim_Xy(n=60)
    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0)
    clf.fit(X, y, strike=strike, epochs=2, verbose=False)
    proba_before = clf.transform(X)

    path = tmp_path / "dimclf_torch.npz"
    clf.save(path)
    loaded = DimensionalityClassifier.load(path)

    proba_after = loaded.transform(X)
    np.testing.assert_allclose(proba_before, proba_after, rtol=1e-4, atol=1e-5)
    np.testing.assert_array_equal(loaded.predict(X), clf.predict(X))


def test_dim_classifier_save_load_round_trip_rf(tmp_path, monkeypatch):
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.classify as classify_mod

    monkeypatch.setattr(
        classify_mod.DimensionalityClassifier,
        "_fit_torch",
        lambda self, Xn, y, strike, **kwargs: (_ for _ in ()).throw(
            ImportError("no torch")
        ),
    )
    X, y, _ = _dim_Xy(n=60)
    clf = classify_mod.DimensionalityClassifier()
    clf.fit(X, y, verbose=False)
    assert clf._use_rf is True
    proba_before = clf.transform(X)

    path = tmp_path / "dimclf_rf.npz"
    clf.save(path)
    loaded = classify_mod.DimensionalityClassifier.load(path)

    assert loaded._use_rf is True
    proba_after = loaded.transform(X)
    np.testing.assert_allclose(proba_before, proba_after)
    assert "rf" in repr(loaded)


def test_dim_classifier_rf_proba_maps_missing_class_columns():
    """RandomForestClassifier.predict_proba only returns a column per
    class actually seen in training -- with a class entirely absent
    from ``y``, ``_predict_proba`` must still return an
    (n_samples, n_classes) array with that column near-zero rather
    than raising a shape-mismatch error."""
    pytest.importorskip("sklearn")
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, _ = _dim_Xy(n=60)
    y = np.where(y == 0, 1, y)  # drop class 0 ("1D") entirely
    assert set(np.unique(y)) == {1, 2}

    clf = DimensionalityClassifier()
    Xn = (X - X.mean(0, keepdims=True)) / (X.std(0, keepdims=True) + 1e-8)
    clf._x_mean = X.mean(0, keepdims=True)
    clf._x_std = X.std(0, keepdims=True) + 1e-8
    clf._fit_rf(Xn, y, verbose=False)
    clf._use_rf = True
    clf._is_fitted = True

    proba = clf.transform(X)
    assert proba.shape == (len(X), clf.n_classes)
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-6)
    assert np.all(proba[:, 0] < 1e-6)  # missing class carries no mass


def test_dim_classifier_predict_table_from_sites(sites):
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, _ = _dim_Xy(n=60)
    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0)
    clf.fit(X, y, epochs=2, verbose=False)

    tbl = clf.predict_table(sites)
    assert {
        "station",
        "freq",
        "period",
        "dim",
        "dim_label",
        "strike",
        "confidence",
    }.issubset(tbl.columns)
    assert len(tbl) > 0
    assert set(tbl["dim_label"].unique()).issubset({"1D", "2D", "3D"})


# ─────────────────────────────────────────────────────────────────────────────
# processing/plot.py — plot_qc_feature_heatmap, plot_qc_summary
# ─────────────────────────────────────────────────────────────────────────────


def _qc_full_dataframe():
    stations = [f"S{i:02d}" for i in range(6)]
    freqs = np.array([1000.0, 100.0, 10.0, 1.0])
    rng = np.random.default_rng(3)
    rows = []
    for st in stations:
        for fr in freqs:
            rows.append(
                {
                    "station": st,
                    "freq": fr,
                    "snr": rng.uniform(2, 30),
                    "swift_skew": rng.uniform(0, 0.5),
                    "asym": rng.uniform(-0.5, 0.5),
                    "phase_xy": rng.uniform(20, 50),
                    "phase_yx": rng.uniform(-160, -120),
                    "score": rng.uniform(0, 1),
                }
            )
    return pd.DataFrame(rows)


def test_plot_qc_feature_heatmap_default_and_overrides():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.ai.processing.plot import plot_qc_feature_heatmap

    df = _qc_full_dataframe()

    fig = plot_qc_feature_heatmap(df, title="QC features")
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) >= 5
    plt.close(fig)

    fig2 = plot_qc_feature_heatmap(
        df,
        features=["snr", "asym"],
        cmaps={"snr": "viridis"},
        n_yticks=3,
    )
    assert len(fig2.axes) >= 2
    plt.close(fig2)

    # single-feature branch (n_feat == 1 -> axes wrapped in a list)
    fig3 = plot_qc_feature_heatmap(df, features=["snr"])
    plt.close(fig3)

    with pytest.raises(ValueError, match="No valid feature columns"):
        plot_qc_feature_heatmap(df, features=["not_a_column"])


def test_plot_qc_summary_dict_and_dataframe_inputs():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.ai.processing.plot import plot_qc_summary

    scores = {
        "L18": np.array([0.2, 0.6, 0.9, 0.4]),
        "L22": np.array([0.7, 0.8, 0.3]),
    }
    fig = plot_qc_summary(scores, suptitle="QC summary", spread_kind="box")
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) >= 3
    plt.close(fig)

    df = _qc_full_dataframe()
    fig2 = plot_qc_summary(df, spread_kind="strip")
    plt.close(fig2)


# ─────────────────────────────────────────────────────────────────────────────
# distortion.py
# ─────────────────────────────────────────────────────────────────────────────


def _distortion_Xy(n=60, seed=11):
    rng = np.random.default_rng(seed)
    X = np.zeros((n, 6), dtype=np.float32)
    X[:, 0] = rng.uniform(0, 5, n)  # beta_abs
    X[:, 1] = rng.uniform(0, 1, n)  # ellipt_abs
    X[:, 2] = rng.normal(0, 0.3, n)  # delta_log10_rho
    X[:, 3] = rng.uniform(-30, 30, n)  # twist_deg
    X[:, 4] = rng.uniform(-0.5, 0.5, n)  # shear
    X[:, 5] = rng.uniform(-0.1, 0.1, n)  # anisotropy
    from pycsamt.ai.processing.distortion import _rule_labels

    y = _rule_labels(X[:, 2], X[:, 3], X[:, 4])
    return X, y


def test_distortion_classifier_fit_transform_predict_torch():
    from pycsamt.ai.processing.distortion import DistortionTypeClassifier

    X, y = _distortion_Xy()
    clf = DistortionTypeClassifier(hidden=(16, 8), dropout=0.0)
    clf.fit(X, y, epochs=5, verbose=False)
    assert clf._backend_name == "torch"

    proba = clf.transform(X)
    assert proba.shape == (len(X), clf.n_classes)
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-5)

    labels = clf.predict(X)
    assert set(np.unique(labels)).issubset({0, 1, 2})
    assert "torch" in repr(clf)
    assert clf.history_["train_loss"]


def test_distortion_classifier_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.distortion import DistortionTypeClassifier

    X, y = _distortion_Xy()
    clf = DistortionTypeClassifier(hidden=(16, 8), dropout=0.0)
    clf.fit(X, y, epochs=5, verbose=False)
    proba_before = clf.transform(X)

    path = tmp_path / "distortion_torch.npz"
    clf.save(path)
    loaded = DistortionTypeClassifier.load(path)
    proba_after = loaded.transform(X)
    np.testing.assert_allclose(proba_before, proba_after)


def test_distortion_classifier_rf_proba_maps_missing_class_columns():
    """Same RandomForestClassifier.predict_proba column-mapping bug as
    DimensionalityClassifier's RF fallback: with a class entirely
    absent from ``y``, ``_predict_proba`` must still return an
    (n_samples, n_classes) array rather than raising a shape-mismatch
    error."""
    pytest.importorskip("sklearn")
    from pycsamt.ai.processing.distortion import DistortionTypeClassifier

    X, y = _distortion_Xy()
    y = np.where(y == 0, 1, y)  # drop class 0 ("clean") entirely
    assert set(np.unique(y)) == {1, 2}

    clf = DistortionTypeClassifier()
    clf._x_mean = X.mean(0, keepdims=True)
    clf._x_std = X.std(0, keepdims=True) + 1e-8
    Xn = (X - clf._x_mean) / clf._x_std
    clf._fit_rf(Xn, y, verbose=False)
    clf._use_rf = True
    clf._is_fitted = True

    proba = clf.transform(X)
    assert proba.shape == (len(X), clf.n_classes)
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-6)
    assert np.all(proba[:, 0] < 1e-6)  # missing class carries no mass


def test_distortion_classifier_predict_table_from_sites(sites):
    from pycsamt.ai.processing.distortion import DistortionTypeClassifier

    X, y = _distortion_Xy()
    clf = DistortionTypeClassifier(hidden=(16, 8), dropout=0.0)
    clf.fit(X, y, epochs=3, verbose=False)

    tbl = clf.predict_table(sites)
    if tbl.empty:
        pytest.skip("3edis dataset lacks a full GB + SS + phase overlap.")
    assert {"station", "regime", "regime_label", "confidence"}.issubset(
        tbl.columns
    )
    assert set(tbl["regime_label"].unique()).issubset(
        {"clean", "static_shift_only", "distorted"}
    )


# ─────────────────────────────────────────────────────────────────────────────
# uncertainty.py
# ─────────────────────────────────────────────────────────────────────────────


def _uncertainty_Xy_table(n=200, seed=13):
    """A synthetic features table shaped like
    build_uncertainty_features_table's output, with a genuine (not
    circular) relationship between the four model features and
    z_err_frac -- larger swift_skew/asym/phase spread -> larger
    error -- so the regressor has real signal to find."""
    rng = np.random.default_rng(seed)
    swift = rng.uniform(0, 1, n)
    asym = rng.uniform(-1, 1, n)
    phase_xy = rng.uniform(0, 90, n)
    phase_yx = rng.uniform(-180, -90, n)
    noise_level = 0.02 + 0.08 * swift + 0.01 * np.abs(asym)
    z_err_frac = np.clip(
        noise_level + rng.normal(0, 0.005, n), 0.005, None
    )
    # snr deliberately built as 1/z_err_frac, matching the real
    # circularity build_uncertainty_features_table produces, so the
    # excluded-column regression test below is a faithful check.
    snr = 1.0 / z_err_frac
    return pd.DataFrame(
        {
            "station": [f"S{i:03d}" for i in range(n)],
            "freq": rng.uniform(1, 1e4, n),
            "snr": snr,
            "swift_skew": swift,
            "asym": asym,
            "phase_xy": phase_xy,
            "phase_yx": phase_yx,
            "z_err_frac": z_err_frac,
        }
    )


def test_uncertainty_calibrator_excludes_snr_feature():
    """snr = amp/err and z_err_frac = err/amp share the same amp/err
    terms in the real feature table -- snr must never reach the
    regressor's input matrix."""
    from pycsamt.ai.processing.uncertainty import (
        _FEATURE_COLS,
        _TABLE_COLS,
        UncertaintyCalibrator,
    )

    assert "snr" not in _FEATURE_COLS
    assert "snr" in _TABLE_COLS
    assert UncertaintyCalibrator().n_features == len(_FEATURE_COLS) == 4

    df = _uncertainty_Xy_table(n=20)
    cal = UncertaintyCalibrator()
    X, _y = cal._coerce_Xy(df, None)
    assert X.shape[1] == 4


def test_uncertainty_calibrator_fit_transform_torch():
    from pycsamt.ai.processing.uncertainty import UncertaintyCalibrator

    df = _uncertainty_Xy_table()
    cal = UncertaintyCalibrator(hidden=(16, 8), dropout=0.0)
    cal.fit(df, epochs=40, seed=0, verbose=False)
    assert cal._backend_name == "torch"

    pred = cal.transform(df)
    assert pred.shape == (len(df),)
    assert np.all(pred > 0)  # fractional error must stay positive

    corr = np.corrcoef(
        np.log10(df["z_err_frac"].to_numpy()), np.log10(pred)
    )[0, 1]
    assert corr > 0.5, f"expected real predictive signal, got corr={corr}"
    assert "torch" in repr(cal)
    assert cal.history_["train_loss"]


def test_uncertainty_calibrator_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.uncertainty import UncertaintyCalibrator

    df = _uncertainty_Xy_table()
    cal = UncertaintyCalibrator(hidden=(16, 8), dropout=0.0)
    cal.fit(df, epochs=20, seed=0, verbose=False)
    pred_before = cal.transform(df)

    path = tmp_path / "uncertainty_torch.npz"
    cal.save(path)
    loaded = UncertaintyCalibrator.load(path)
    pred_after = loaded.transform(df)
    np.testing.assert_allclose(pred_before, pred_after)
    assert loaded._log_target == cal._log_target


def test_uncertainty_calibrator_save_load_round_trip_rf(tmp_path, monkeypatch):
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.uncertainty as unc_mod

    monkeypatch.setattr(
        unc_mod.UncertaintyCalibrator,
        "_fit_torch",
        lambda self, Xn, y, **kwargs: (_ for _ in ()).throw(
            ImportError("no torch")
        ),
    )
    df = _uncertainty_Xy_table()
    cal = unc_mod.UncertaintyCalibrator()
    cal.fit(df, epochs=10, verbose=False)
    assert cal._use_rf is True
    pred_before = cal.transform(df)

    path = tmp_path / "uncertainty_rf.npz"
    cal.save(path)
    loaded = unc_mod.UncertaintyCalibrator.load(path)
    assert loaded._use_rf is True
    pred_after = loaded.transform(df)
    np.testing.assert_allclose(pred_before, pred_after)
    assert "rf" in repr(loaded)


def test_uncertainty_calibrator_fit_requires_target():
    from pycsamt.ai.processing.uncertainty import UncertaintyCalibrator

    cal = UncertaintyCalibrator()
    with pytest.raises(ValueError, match="z_err_frac"):
        cal.fit(np.random.default_rng(0).standard_normal((10, 4)))

    df = _uncertainty_Xy_table(n=10).drop(columns=["z_err_frac"])
    with pytest.raises(ValueError, match="z_err_frac"):
        cal.fit(df)


def test_uncertainty_calibrator_predict_table_and_apply(sites):
    from pycsamt.ai.processing.uncertainty import (
        UncertaintyCalibrator,
        build_uncertainty_features_table,
    )

    feats = build_uncertainty_features_table(sites)
    if feats.empty:
        pytest.skip("3edis dataset has no usable z_err.")

    cal = UncertaintyCalibrator(hidden=(16, 8), dropout=0.0)
    cal.fit(feats, epochs=15, seed=0, verbose=False)

    table = cal.predict_table(sites)
    assert {"station", "freq", "z_err_frac", "z_err_frac_calibrated"}.issubset(
        table.columns
    )
    assert (table["z_err_frac_calibrated"] > 0).all()

    calibrated = cal.apply(sites, inplace=False)
    assert len(list(calibrated)) == len(list(sites))


if __name__ == "__main__":
    pytest.main([__file__])

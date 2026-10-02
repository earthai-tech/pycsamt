# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional coverage for pycsamt.ai.processing.distortion.

test_ai_processing_coverage.py already exercises
``DistortionTypeClassifier``'s torch fit/transform/predict/save/load round
trip, the random-forest missing-class-column mapping, and
``predict_table`` against the real 3edis sites. This file fills the
remaining gap: ``build_distortion_features_table``'s early-return
branches (each upstream table being empty/missing columns) and its
``ImportError`` branch, the DataFrame-input path through ``fit``/
``from_features_table``, the RF backend reached through ``fit`` itself
(not only via the private ``_fit_rf`` helper), the RF save/load round
trip, ``transform``/``predict`` called before ``fit``, and ``__repr__``.

Per the task brief: ``DistortionTypeClassifier``'s self-training labels
are known to flip run-to-run on tiny real surveys (a diagnosed, accepted
instability, not a bug to fix here) -- assertions below check structural
properties (shapes, valid probability simplices, membership in the
label set) rather than pinning a specific predicted class.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.ai.processing.distortion import (
    DistortionTypeClassifier,
    _DISTORTION_LABELS,
    _FEATURE_COLS,
    build_distortion_features_table,
)


def _distortion_Xy(n=60, seed=11):
    rng = np.random.default_rng(seed)
    X = np.zeros((n, 6), dtype=np.float32)
    X[:, 0] = rng.uniform(0, 5, n)
    X[:, 1] = rng.uniform(0, 1, n)
    X[:, 2] = rng.normal(0, 0.3, n)
    X[:, 3] = rng.uniform(-30, 30, n)
    X[:, 4] = rng.uniform(-0.5, 0.5, n)
    X[:, 5] = rng.uniform(-0.1, 0.1, n)
    from pycsamt.ai.processing.distortion import _rule_labels

    y = _rule_labels(X[:, 2], X[:, 3], X[:, 4])
    return X, y


def _distortion_df(n=40, seed=3):
    X, y = _distortion_Xy(n=n, seed=seed)
    df = pd.DataFrame(X, columns=_FEATURE_COLS)
    df.insert(0, "station", [f"S{i:03d}" for i in range(n)])
    return df


# ---------------------------------------------------------------------------
# build_distortion_features_table — early-return branches
# ---------------------------------------------------------------------------


def test_build_features_table_empty_phase_df(monkeypatch):
    import pycsamt.emtools.dimensionality as dim_mod

    monkeypatch.setattr(
        dim_mod, "phase_features_table", lambda *a, **kw: pd.DataFrame()
    )
    df = build_distortion_features_table([object()])
    assert list(df.columns) == ["station", *_FEATURE_COLS]
    assert df.empty


def test_build_features_table_empty_gb_df(monkeypatch):
    import pycsamt.emtools.dimensionality as dim_mod
    import pycsamt.emtools.gb as gb_mod

    monkeypatch.setattr(
        dim_mod,
        "phase_features_table",
        lambda *a, **kw: pd.DataFrame(
            {"station": ["A"], "beta_abs": [1.0], "ellipt_abs": [0.5]}
        ),
    )
    monkeypatch.setattr(gb_mod, "groom_bailey_table", lambda *a, **kw: pd.DataFrame())
    df = build_distortion_features_table([object()])
    assert df.empty


def test_build_features_table_gb_missing_status_column(monkeypatch):
    import pycsamt.emtools.dimensionality as dim_mod
    import pycsamt.emtools.gb as gb_mod

    monkeypatch.setattr(
        dim_mod,
        "phase_features_table",
        lambda *a, **kw: pd.DataFrame(
            {"station": ["A"], "beta_abs": [1.0], "ellipt_abs": [0.5]}
        ),
    )
    monkeypatch.setattr(
        gb_mod,
        "groom_bailey_table",
        lambda *a, **kw: pd.DataFrame({"station": ["A"], "twist_deg": [1.0]}),
    )
    df = build_distortion_features_table([object()])
    assert df.empty


def test_build_features_table_gb_all_status_not_ok(monkeypatch):
    import pycsamt.emtools.dimensionality as dim_mod
    import pycsamt.emtools.gb as gb_mod

    monkeypatch.setattr(
        dim_mod,
        "phase_features_table",
        lambda *a, **kw: pd.DataFrame(
            {"station": ["A"], "beta_abs": [1.0], "ellipt_abs": [0.5]}
        ),
    )
    monkeypatch.setattr(
        gb_mod,
        "groom_bailey_table",
        lambda *a, **kw: pd.DataFrame(
            {
                "station": ["A"],
                "twist_deg": [1.0],
                "shear": [0.1],
                "anisotropy": [0.01],
                "status": ["failed"],
            }
        ),
    )
    df = build_distortion_features_table([object()])
    assert df.empty


def test_build_features_table_empty_ss_df(monkeypatch):
    import pycsamt.emtools.dimensionality as dim_mod
    import pycsamt.emtools.gb as gb_mod
    import pycsamt.emtools.ss as ss_mod

    monkeypatch.setattr(
        dim_mod,
        "phase_features_table",
        lambda *a, **kw: pd.DataFrame(
            {"station": ["A"], "beta_abs": [1.0], "ellipt_abs": [0.5]}
        ),
    )
    monkeypatch.setattr(
        gb_mod,
        "groom_bailey_table",
        lambda *a, **kw: pd.DataFrame(
            {
                "station": ["A"],
                "twist_deg": [1.0],
                "shear": [0.1],
                "anisotropy": [0.01],
                "status": ["ok"],
            }
        ),
    )
    monkeypatch.setattr(ss_mod, "estimate_ss_ama", lambda *a, **kw: pd.DataFrame())
    df = build_distortion_features_table([object()])
    assert df.empty


def test_build_features_table_full_merge_happy_path(monkeypatch):
    import pycsamt.emtools.dimensionality as dim_mod
    import pycsamt.emtools.gb as gb_mod
    import pycsamt.emtools.ss as ss_mod

    monkeypatch.setattr(
        dim_mod,
        "phase_features_table",
        lambda *a, **kw: pd.DataFrame(
            {
                "station": ["A", "A", "B", "B"],
                "beta_abs": [1.0, 2.0, 3.0, 4.0],
                "ellipt_abs": [0.1, 0.2, 0.3, 0.4],
            }
        ),
    )
    monkeypatch.setattr(
        gb_mod,
        "groom_bailey_table",
        lambda *a, **kw: pd.DataFrame(
            {
                "station": ["A", "B"],
                "twist_deg": [1.0, -2.0],
                "shear": [0.1, -0.1],
                "anisotropy": [0.01, 0.02],
                "status": ["ok", "ok"],
            }
        ),
    )
    monkeypatch.setattr(
        ss_mod,
        "estimate_ss_ama",
        lambda *a, **kw: pd.DataFrame(
            {"station": ["A", "B"], "delta_log10_rho": [0.05, -0.2]}
        ),
    )
    df = build_distortion_features_table([object()])
    assert list(df.columns) == ["station", *_FEATURE_COLS]
    assert set(df["station"]) == {"A", "B"}
    # median beta_abs per station
    a_row = df[df["station"] == "A"].iloc[0]
    assert a_row["beta_abs"] == pytest.approx(1.5)


def test_build_features_table_import_error(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def _fake_import(name, *a, **kw):
        if name == "pycsamt.emtools.dimensionality":
            raise ImportError("blocked")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", _fake_import)
    with pytest.raises(ImportError, match="emtools is required"):
        build_distortion_features_table([object()])


# ---------------------------------------------------------------------------
# fit()/transform()/predict() guard rails
# ---------------------------------------------------------------------------


def test_transform_before_fit_raises():
    clf = DistortionTypeClassifier()
    with pytest.raises(RuntimeError, match="Call fit"):
        clf.transform(_distortion_Xy()[0])


def test_repr_unfitted():
    clf = DistortionTypeClassifier()
    assert "unfitted" in repr(clf)
    assert "torch" in repr(clf)  # default backend guess before any fit


# ---------------------------------------------------------------------------
# DataFrame input path (fit / from_features_table / predict_table empty)
# ---------------------------------------------------------------------------


def test_from_features_table_classmethod_trains_a_classifier():
    df = _distortion_df()
    clf = DistortionTypeClassifier.from_features_table(
        df, epochs=3, verbose=False
    )
    assert clf._is_fitted
    proba = clf.transform(df)
    assert proba.shape == (len(df), clf.n_classes)
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-5)
    labels = clf.predict(df)
    assert set(np.unique(labels)).issubset(set(range(clf.n_classes)))


def test_fit_with_dataframe_x_and_explicit_y_overrides_rule_labels():
    df = _distortion_df(n=30)
    X, _ = _distortion_Xy(n=30, seed=3)
    forced_y = np.zeros(30, dtype=int)  # force every label to "clean"
    clf = DistortionTypeClassifier(hidden=(8, 4))
    clf.fit(df, forced_y, epochs=3, verbose=False)
    assert clf._is_fitted


def test_from_features_table_with_label_col():
    df = _distortion_df(n=25)
    df["my_label"] = 2  # force every label to "distorted"
    clf = DistortionTypeClassifier.from_features_table(
        df, label_col="my_label", epochs=3, verbose=False
    )
    assert clf._is_fitted


def test_predict_table_returns_empty_df_when_features_are_empty(monkeypatch):
    import pycsamt.ai.processing.distortion as distortion_mod

    monkeypatch.setattr(
        distortion_mod,
        "build_distortion_features_table",
        lambda *a, **kw: pd.DataFrame(columns=["station", *_FEATURE_COLS]),
    )
    clf = DistortionTypeClassifier(hidden=(8, 4))
    X, y = _distortion_Xy(n=20)
    clf.fit(X, y, epochs=2, verbose=False)
    out = clf.predict_table([object()])
    assert out.empty


# ---------------------------------------------------------------------------
# RandomForest backend reached through fit() itself
# ---------------------------------------------------------------------------


def test_fit_falls_back_to_rf_when_no_dl_backend(monkeypatch):
    # fit() only branches on active_backend() to pick tensorflow vs. torch
    # (the same convention as classify.py's DimensionalityClassifier) --
    # it does not special-case a "none" active_backend() the way
    # EMDenoiser.fit() does. With real PyTorch installed in this test
    # environment, the RF fallback is only actually reached when building
    # the torch network itself fails, so that failure is forced directly.
    import pycsamt.ai.processing.distortion as distortion_mod

    def _boom(*a, **kw):
        raise ImportError("PyTorch is required for DistortionTypeClassifier")

    monkeypatch.setattr(distortion_mod, "_build_distortion_mlp_torch", _boom)
    X, y = _distortion_Xy(n=40)
    clf = DistortionTypeClassifier()
    clf.fit(X, y, verbose=False)
    assert clf._use_rf is True
    assert clf._backend_name == "rf"
    assert "rf" in repr(clf)

    proba = clf.transform(X)
    assert proba.shape == (len(X), clf.n_classes)
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-6)


def test_fit_rf_import_error_when_sklearn_unavailable(monkeypatch):
    import builtins

    import pycsamt.ai.processing.distortion as distortion_mod

    def _boom(*a, **kw):
        raise ImportError("PyTorch is required for DistortionTypeClassifier")

    monkeypatch.setattr(distortion_mod, "_build_distortion_mlp_torch", _boom)

    real_import = builtins.__import__

    def _fake_import(name, *a, **kw):
        if name == "sklearn.ensemble":
            raise ImportError("blocked")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", _fake_import)
    X, y = _distortion_Xy(n=20)
    clf = DistortionTypeClassifier()
    with pytest.raises(ImportError, match="PyTorch, TensorFlow, or scikit-learn"):
        clf.fit(X, y, verbose=False)


def test_rf_save_load_round_trip(tmp_path):
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.distortion as distortion_mod

    X, y = _distortion_Xy(n=40)
    clf = DistortionTypeClassifier()
    clf._x_mean = X.mean(0, keepdims=True)
    clf._x_std = X.std(0, keepdims=True) + 1e-8
    Xn = (X - clf._x_mean) / clf._x_std
    clf._fit_rf(Xn, y, verbose=True)
    clf._use_rf = True
    clf._backend_name = "rf"
    clf._is_fitted = True
    proba_before = clf.transform(X)

    path = tmp_path / "distortion_rf.npz"
    clf.save(path)
    loaded = distortion_mod.DistortionTypeClassifier.load(path)
    assert loaded._use_rf is True
    proba_after = loaded.transform(X)
    np.testing.assert_allclose(proba_before, proba_after)


def test_predict_proba_rf_all_nan_returns_uniform():
    X, y = _distortion_Xy(n=10)
    clf = DistortionTypeClassifier()
    clf._x_mean = X.mean(0, keepdims=True)
    clf._x_std = X.std(0, keepdims=True) + 1e-8
    Xn = (X - clf._x_mean) / clf._x_std
    clf._fit_rf(Xn, y, verbose=False)
    clf._use_rf = True
    clf._is_fitted = True

    X_nan = np.full_like(X, np.nan)
    proba = clf.transform(X_nan)
    assert np.allclose(proba, 1.0 / clf.n_classes)


# ---------------------------------------------------------------------------
# _coerce_Xy — short ndarray column count falls back to all-zero labels
# ---------------------------------------------------------------------------


def test_coerce_xy_short_columns_default_zero_labels():
    clf = DistortionTypeClassifier(hidden=(8, 4))
    X = np.random.default_rng(0).normal(size=(15, 3)).astype(np.float32)
    clf.fit(X, epochs=2, verbose=False)
    assert clf._is_fitted


def test_distortion_labels_constant_matches_classes():
    assert len(_DISTORTION_LABELS) == 3

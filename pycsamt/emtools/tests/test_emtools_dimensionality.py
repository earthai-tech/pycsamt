"""Tests for pycsamt.emtools.dimensionality"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from pycsamt.api import reset_api_view
from pycsamt.emtools.dimensionality import (
    classify_dimensionality,
    encode_dimensionality,
    mask_by_dimensionality,
    phase_features_table,
)


@pytest.fixture(autouse=True)
def _no_api_view():
    reset_api_view()
    yield
    reset_api_view()


# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """Pure 1-D/2-D isotropic Z (zero diagonal)."""
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _3d_z(freqs: np.ndarray, skew_frac: float = 0.5) -> np.ndarray:
    """Z with non-zero diagonal (3-D character)."""
    z = _iso_z(freqs)
    amp = np.abs(z[:, 0, 1])
    z[:, 0, 0] = skew_frac * amp * (0.6 + 0.8j)
    z[:, 1, 1] = -skew_frac * amp * (0.5 + 0.7j)
    return z


def _iso_site(name: str, n: int = 10) -> _FakeSite:
    fr = _freqs(n)
    return _FakeSite(name, _iso_z(fr), fr)


def _3d_site(name: str, n: int = 10) -> _FakeSite:
    fr = _freqs(n)
    return _FakeSite(name, _3d_z(fr), fr)


# ─────────────────────────────────────────────────────────────────────────────
# phase_features_table
# ─────────────────────────────────────────────────────────────────────────────


class TestPhaseFeatureTable:
    def test_returns_dataframe(self):
        import pandas as pd

        sites = [_iso_site("S00")]
        df = phase_features_table(sites, api=False)
        assert isinstance(df, pd.DataFrame)

    def test_columns_present(self):
        sites = [_iso_site("S00")]
        df = phase_features_table(sites, api=False)
        for col in ("station", "freq", "period", "beta_abs", "ellipt_abs"):
            assert col in df.columns

    def test_row_count(self):
        n = 8
        sites = [_iso_site("S00", n=n)]
        df = phase_features_table(sites, api=False)
        assert len(df) == n

    def test_multiple_sites(self):
        sites = [_iso_site(f"S{i:02d}") for i in range(4)]
        df = phase_features_table(sites, api=False)
        assert df["station"].nunique() == 4

    def test_empty_input_returns_empty(self):
        import pandas as pd

        from pycsamt.api.view.frame import APIFrame

        df = phase_features_table([])
        assert isinstance(df, (pd.DataFrame, APIFrame))

    def test_period_positive(self):
        sites = [_iso_site("S00")]
        df = phase_features_table(sites, api=False)
        assert (df["period"].dropna() > 0).all()

    def test_freq_positive(self):
        sites = [_iso_site("S00")]
        df = phase_features_table(sites, api=False)
        assert (df["freq"].dropna() > 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# classify_dimensionality
# ─────────────────────────────────────────────────────────────────────────────


class TestClassifyDimensionality:
    def test_returns_dataframe(self):
        import pandas as pd

        sites = [_iso_site("S00")]
        df = classify_dimensionality(sites, api=False)
        assert isinstance(df, pd.DataFrame)

    def test_columns_present(self):
        sites = [_iso_site("S00")]
        df = classify_dimensionality(sites, api=False)
        assert "station" in df.columns

    def test_empty_input(self):
        import pandas as pd

        from pycsamt.api.view.frame import APIFrame

        df = classify_dimensionality([])
        assert isinstance(df, (pd.DataFrame, APIFrame))

    def test_multiple_sites(self):
        sites = [_iso_site(f"S{i:02d}") for i in range(3)]
        df = classify_dimensionality(sites, api=False)
        assert df["station"].nunique() == 3

    def test_high_skew_threshold(self):
        """Very high skew_th → even 3-D stations pass classification."""
        sites = [_iso_site("S00"), _3d_site("S01")]
        df = classify_dimensionality(sites, skew_th=100.0, api=False)
        assert len(df) > 0

    def test_class_column_when_present(self):
        sites = [_iso_site("S00"), _3d_site("S01")]
        df = classify_dimensionality(sites, api=False)
        if "class" in df.columns:
            assert df["class"].notna().any()


# ─────────────────────────────────────────────────────────────────────────────
# mask_by_dimensionality
# ─────────────────────────────────────────────────────────────────────────────


class TestMaskByDimensionality:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        sites = [_iso_site("S00"), _iso_site("S01")]
        result = mask_by_dimensionality(sites)
        assert isinstance(result, Sites)

    def test_site_count_preserved(self):
        sites = [_iso_site(f"S{i:02d}") for i in range(3)]
        result = mask_by_dimensionality(sites)
        assert sum(1 for _ in result) == 3

    def test_keep_all_classes(self):
        """keep=(0, 1) is the default — all dimensionality classes pass."""
        from pycsamt.site.base import Sites

        sites = [_3d_site("S00")]
        result = mask_by_dimensionality(sites, keep=(0, 1))
        assert isinstance(result, Sites)

    def test_empty_input(self):
        from pycsamt.site.base import Sites

        result = mask_by_dimensionality([])
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# encode_dimensionality
# ─────────────────────────────────────────────────────────────────────────────


class TestEncodeDimensionality:
    def test_empty_model_returns_dataframe(self):
        import pandas as pd

        from pycsamt.api.view.frame import APIFrame

        sites = [_iso_site("S00")]
        df = encode_dimensionality(sites, {})
        assert isinstance(df, (pd.DataFrame, APIFrame))

    def test_empty_model_empty_df(self):
        sites = [_iso_site("S00")]
        df = encode_dimensionality(sites, {}, api=False)
        assert df.empty

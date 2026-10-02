"""Tests for pycsamt.emtools.dimensionality"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pycsamt.api import reset_api_view
from pycsamt.emtools.dimensionality import (
    classify_dimensionality,
    encode_dimensionality,
    learn_dim_dictionary,
    mask_by_dictionary,
    mask_by_dimensionality,
    phase_features_table,
    plot_atom_psection,
    plot_dim_confidence_grid,
    plot_dim_map,
    plot_dim_occupancy_area,
    pre2d_inversion_assessment,
    project_to_2d,
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
    def __init__(self, station, z, freq, coords=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if coords is not None:
            self.coords = coords

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

    def test_trained_model_encodes(self):
        sites = [_iso_site(f"S{i:02d}", n=12) for i in range(3)] + [
            _3d_site("S03", n=12)
        ]
        model = learn_dim_dictionary(sites, n_atoms=4, n_iter=3, code_iter=5)
        assert model["D"] is not None
        df = encode_dimensionality(sites, model, code_iter=5, api=False)
        assert "dim_pred" in df.columns
        assert len(df) > 0

    def test_empty_sites_encoding(self):
        model = dict(D=None, A=None, mu=None, sd=None, feat=[])
        df = encode_dimensionality([], model, api=False)
        assert df.empty


# ─────────────────────────────────────────────────────────────────────────────
# learn_dim_dictionary
# ─────────────────────────────────────────────────────────────────────────────


class TestLearnDimDictionary:
    def test_empty_input_returns_none_model(self):
        model = learn_dim_dictionary([])
        assert model["D"] is None
        assert model["meta"]["samples"] == 0

    def test_nonempty_learns_dictionary(self):
        sites = [_iso_site(f"S{i:02d}", n=10) for i in range(2)]
        model = learn_dim_dictionary(sites, n_atoms=3, n_iter=2, code_iter=5)
        assert model["D"].shape[0] == 4  # 4 features
        assert model["A"].shape[1] == model["meta"]["samples"]


# ─────────────────────────────────────────────────────────────────────────────
# mask_by_dictionary
# ─────────────────────────────────────────────────────────────────────────────


class TestMaskByDictionary:
    def test_empty_model_returns_sites_unchanged(self):
        from pycsamt.site.base import Sites

        sites = [_iso_site("S00")]
        result = mask_by_dictionary(sites, {})
        assert isinstance(result, Sites)

    def test_trained_model_masks(self):
        from pycsamt.site.base import Sites

        sites = [_iso_site(f"S{i:02d}", n=10) for i in range(2)]
        model = learn_dim_dictionary(sites, n_atoms=3, n_iter=2, code_iter=5)
        result = mask_by_dictionary(sites, model, code_iter=5)
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# pre2d_inversion_assessment
# ─────────────────────────────────────────────────────────────────────────────


class TestPre2dInversionAssessment:
    def test_returns_dataframe_with_recommendation(self):
        sites = [_iso_site(f"S{i:02d}", n=12) for i in range(3)]
        df = pre2d_inversion_assessment(sites, api=False)
        assert "recommendation" in df.columns
        assert len(df) == 3

    def test_band_filters_rows(self):
        sites = [_iso_site("S00", n=20)]
        df = pre2d_inversion_assessment(sites, band=(1e-3, 1e-1), api=False)
        assert "period_min_s" in df.columns

    def test_rotation_and_groom_bailey_flags(self):
        sites = [_iso_site("S00", n=12)]
        df = pre2d_inversion_assessment(
            sites,
            rotation_applied=True,
            rotation_method="phase_tensor",
            groom_bailey_attempted=True,
            groom_bailey_applied=True,
            groom_bailey_reason="custom reason",
            api=False,
        )
        row = df.iloc[0]
        assert row["rotated_to_strike"] is True or bool(row["rotated_to_strike"])
        assert row["groom_bailey_reason"] == "custom reason"

    def test_3d_station_recommends_review(self):
        sites = [_3d_site("S00", n=20)]
        df = pre2d_inversion_assessment(sites, skew_th=0.001, api=False)
        assert df.iloc[0]["recommendation"] in (
            "review_3d_effects_before_2d",
            "unstable_strike_review_band",
            "acceptable_for_2d_with_documented_rotation",
        )


# ─────────────────────────────────────────────────────────────────────────────
# project_to_2d
# ─────────────────────────────────────────────────────────────────────────────


class TestProjectTo2d:
    def test_default_rotate_to_strike(self):
        from pycsamt.site.base import Sites

        sites = [_iso_site(f"S{i:02d}", n=10) for i in range(3)]
        result = project_to_2d(sites)
        assert isinstance(result, Sites)

    def test_explicit_strike_angle(self):
        from pycsamt.site.base import Sites

        sites = [_iso_site("S00", n=10)]
        result = project_to_2d(sites, strike=30.0)
        assert isinstance(result, Sites)

    def test_no_antisym(self):
        from pycsamt.site.base import Sites

        sites = [_iso_site("S00", n=10)]
        result = project_to_2d(sites, strike=10.0, antisym=False)
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# plot_atom_psection / plot_dim_confidence_grid / plot_dim_occupancy_area /
# plot_dim_map
# ─────────────────────────────────────────────────────────────────────────────


def _trained(n_sites=3, n=12):
    sites = [_iso_site(f"S{i:02d}", n=n) for i in range(n_sites)]
    model = learn_dim_dictionary(sites, n_atoms=3, n_iter=2, code_iter=5)
    return sites, model


class TestPlotAtomPsection:
    def test_returns_axes(self):
        sites, model = _trained()
        ax = plot_atom_psection(sites, model)
        plt.close("all")
        assert ax is not None

    def test_no_codes_when_model_empty(self):
        sites = [_iso_site("S00")]
        ax = plot_atom_psection(sites, {})
        texts = [t.get_text() for t in ax.texts]
        plt.close("all")
        assert "no codes" in texts

    def test_energy_variants(self):
        sites, model = _trained()
        for energy in ("l1", "max", "l2"):
            ax = plot_atom_psection(sites, model, energy=energy)
            plt.close("all")
            assert ax is not None


class TestPlotDimConfidenceGrid:
    def test_returns_axes(self):
        sites = [_iso_site(f"S{i:02d}", n=12) for i in range(3)]
        ax = plot_dim_confidence_grid(sites)
        plt.close("all")
        assert ax is not None

    def test_empty_sites_no_data(self):
        ax = plot_dim_confidence_grid([])
        texts = [t.get_text() for t in ax.texts]
        plt.close("all")
        assert "no data" in texts


class TestPlotDimOccupancyArea:
    def test_returns_axes(self):
        sites = [_iso_site(f"S{i:02d}", n=12) for i in range(3)]
        ax = plot_dim_occupancy_area(sites)
        plt.close("all")
        assert ax is not None

    def test_empty_sites_no_data(self):
        ax = plot_dim_occupancy_area([])
        texts = [t.get_text() for t in ax.texts]
        plt.close("all")
        assert "no data" in texts


class TestPlotDimMap:
    def test_no_coords_message(self):
        sites = [_iso_site(f"S{i:02d}", n=12) for i in range(2)]
        ax = plot_dim_map(sites)
        texts = [t.get_text() for t in ax.texts]
        plt.close("all")
        assert "no coords" in texts

    def test_with_coords_plots_scatter(self):
        fr = _freqs(12)
        sites = [
            _FakeSite("S00", _iso_z(fr), fr, coords=(10.0, 20.0, 0.0)),
            _FakeSite("S01", _3d_z(fr), fr, coords=(10.1, 20.1, 0.0)),
        ]
        ax = plot_dim_map(sites, period=1.0)
        plt.close("all")
        assert ax is not None

    def test_empty_sites_no_data(self):
        ax = plot_dim_map([])
        texts = [t.get_text() for t in ax.texts]
        plt.close("all")
        assert "no data" in texts

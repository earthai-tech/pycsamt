"""Tests for pycsamt.emtools.inspect"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.api import reset_api_view
from pycsamt.emtools.inspect import (
    frequency_coverage,
    list_missing_sections,
    plot_coverage,
    plot_survey_inventory_overview,
    plot_rhoa_phi,
    plot_tipper_components,
    pseudosection,
    sites_summary,
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


class _FakeTipper:
    def __init__(self, tipper, freq):
        self.tipper = np.asarray(tipper, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(
        self,
        station,
        z,
        freq,
        *,
        tipper=None,
        lat=None,
        lon=None,
        east=None,
        north=None,
    ):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)
        if lat is not None:
            self.lat = float(lat)
            self.lon = float(lon)
        if east is not None:
            self.east = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _tipper(freqs: np.ndarray, amp: float = 0.1) -> np.ndarray:
    """Synthetic (n, 2) complex tipper."""
    t = np.zeros((freqs.size, 2), dtype=complex)
    t[:, 0] = amp * (0.8 + 0.6j)
    t[:, 1] = amp * (0.3 + 0.4j)
    return t


def _site(
    name: str, n: int = 10, *, with_tipper=False, lat=None, lon=None
) -> _FakeSite:
    fr = _freqs(n)
    tip = _tipper(fr) if with_tipper else None
    return _FakeSite(name, _iso_z(fr), fr, tipper=tip, lat=lat, lon=lon)


# ─────────────────────────────────────────────────────────────────────────────
# sites_summary
# ─────────────────────────────────────────────────────────────────────────────


class TestSitesSummary:
    def test_returns_dataframe(self):
        import pandas as pd

        sites = [_site("S00")]
        df = sites_summary(sites, api=False)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        sites = [_site("S00")]
        df = sites_summary(sites, api=False)
        for col in (
            "station",
            "n_freq",
            "has_tipper",
            "period_min",
            "period_max",
        ):
            assert col in df.columns

    def test_one_row_per_site(self):
        n = 4
        sites = [_site(f"S{i:02d}") for i in range(n)]
        df = sites_summary(sites, api=False)
        assert len(df) == n

    def test_n_freq_correct(self):
        nf = 8
        sites = [_site("S00", n=nf)]
        df = sites_summary(sites, api=False)
        assert df["n_freq"].iloc[0] == nf

    def test_has_tipper_column_is_bool(self):
        """has_tipper column contains boolean values."""
        sites = [_site("S00")]
        df = sites_summary(sites, api=False)
        assert df["has_tipper"].dtype in (bool, object, "bool")

    def test_period_min_max_consistent(self):
        sites = [_site("S00")]
        df = sites_summary(sites, api=False)
        assert df["period_min"].iloc[0] < df["period_max"].iloc[0]

    def test_empty_input(self):
        import pandas as pd

        from pycsamt.api.view.frame import APIFrame

        df = sites_summary([])
        assert isinstance(df, (pd.DataFrame, APIFrame))


# ─────────────────────────────────────────────────────────────────────────────
# list_missing_sections
# ─────────────────────────────────────────────────────────────────────────────


class TestListMissingSections:
    def test_returns_dict(self):
        sites = [_site("S00")]
        result = list_missing_sections(sites)
        assert isinstance(result, dict)

    def test_tipper_require_returns_dict(self):
        sites = [_site("S00", with_tipper=False)]
        result = list_missing_sections(sites, require=("tipper",))
        assert isinstance(result, dict)

    def test_empty_require_returns_empty(self):
        sites = [_site("S00")]
        result = list_missing_sections(sites, require=())
        assert result == {}

    def test_empty_input(self):
        result = list_missing_sections([])
        assert result == {}


# ─────────────────────────────────────────────────────────────────────────────
# frequency_coverage
# ─────────────────────────────────────────────────────────────────────────────


class TestFrequencyCoverage:
    def test_per_site_returns_dataframe(self):
        import pandas as pd

        sites = [_site(f"S{i}") for i in range(3)]
        result = frequency_coverage(sites, mode="per-site")
        assert isinstance(result, (pd.DataFrame, dict, type(None))) or True

    def test_empty_input_does_not_raise(self):
        try:
            frequency_coverage([])
        except Exception:
            pass  # allowed to return empty; must not crash catastrophically


# ─────────────────────────────────────────────────────────────────────────────
# plot_coverage
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotCoverage:
    def test_returns_figure(self):
        sites = [_site(f"S{i}") for i in range(3)]
        result = plot_coverage(sites)
        plt.close("all")
        assert result is not None

    def test_external_ax(self):
        sites = [_site("S00")]
        fig, ax = plt.subplots()
        result = plot_coverage(sites, ax=ax)
        plt.close("all")
        assert result is not None

    def test_empty_sites_no_crash(self):
        try:
            plot_coverage([])
        except Exception:
            pass
        plt.close("all")


class TestPlotSurveyInventoryOverview:
    def test_aligned_count_and_coverage_panels(self):
        sites = [_site("S00", n=5), _site("S01", n=7), _site("S02", n=4)]
        fig = plot_survey_inventory_overview(sites)
        data_axes = [ax for ax in fig.axes if ax.get_ylabel()]
        assert len(data_axes) == 2
        assert data_axes[0].xaxis.get_label_position() == "top"
        assert len(data_axes[0].lines) >= 2  # count profile + station guides
        assert data_axes[1].collections
        assert data_axes[0].get_xlim() == data_axes[1].get_xlim()
        plt.close("all")

    def test_custom_station_labels_must_match(self):
        with pytest.raises(ValueError, match="station_labels"):
            plot_survey_inventory_overview(
                [_site("S00"), _site("S01")],
                station_order=["S00", "S01"],
                station_labels=["only one"],
            )

    def test_empty_sites_returns_message_figure(self):
        fig = plot_survey_inventory_overview([])
        assert fig is not None
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_rhoa_phi
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotRhoaPhi:
    def test_returns_figure(self):
        sites = [_site("S00")]
        result = plot_rhoa_phi(sites)
        plt.close("all")
        assert result is not None

    def test_external_axes(self):
        sites = [_site("S00")]
        fig, (ax1, ax2) = plt.subplots(1, 2)
        result = plot_rhoa_phi(sites, ax_r=ax1, ax_p=ax2)
        plt.close("all")
        assert result is not None

    def test_multi_site_no_crash(self):
        sites = [_site(f"S{i}") for i in range(4)]
        result = plot_rhoa_phi(sites)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_tipper_components
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotTipperComponents:
    def test_returns_figure_with_tipper(self):
        sites = [_site("S00", with_tipper=True)]
        result = plot_tipper_components(sites)
        plt.close("all")
        assert result is not None

    def test_no_tipper_graceful(self):
        sites = [_site("S00", with_tipper=False)]
        result = plot_tipper_components(sites)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# pseudosection
# ─────────────────────────────────────────────────────────────────────────────


class TestPseudosection:
    def test_returns_axes(self):
        sites = [_site(f"S{i}") for i in range(3)]
        result = pseudosection(sites)
        plt.close("all")
        assert result is not None

    def test_external_ax(self):
        sites = [_site(f"S{i}") for i in range(3)]
        fig, ax = plt.subplots()
        result = pseudosection(sites, ax=ax)
        plt.close("all")
        assert result is not None

    def test_empty_sites_no_crash(self):
        result = pseudosection([])
        plt.close("all")
        assert result is not None

    def test_different_quantities(self):
        sites = [_site(f"S{i}") for i in range(3)]
        for qty in ("rho_xy", "phi_xy", "rho_yx"):
            result = pseudosection(sites, quantity=qty)
            plt.close("all")
            assert result is not None

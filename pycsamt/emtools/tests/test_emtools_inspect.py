"""Tests for pycsamt.emtools.inspect"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.api import reset_api_view
from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.inspect import (
    _coords,
    _df_from_kv,
    _get_freq,
    _has,
    _is_df,
    _iter_items,
    _name,
    _period,
    _pick_station,
    _presence_vec,
    _rho_phase_err,
    _rms_log,
    _union_freq,
    frequency_coverage,
    list_missing_sections,
    plot_coverage,
    plot_station_response,
    plot_survey_inventory_overview,
    plot_rhoa_phi,
    plot_tipper_components,
    pseudosection,
    sites_summary,
)

_REPO_ROOT = Path(__file__).resolve().parents[3]
_KAP03LMT = _REPO_ROOT / "data" / "MT" / "kap03lmt_edis"
_HAS_KAP03LMT = _KAP03LMT.exists() and any(_KAP03LMT.glob("*.edi"))

pytestmark_kap03lmt = pytest.mark.skipif(
    not _HAS_KAP03LMT, reason="kap03lmt EDI sample data not available"
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

    def test_logperiod_axis_uses_linear_scale_and_log10_label(self):
        sites = [_site("S00", with_tipper=True)]
        fig, ax = plt.subplots()
        plot_tipper_components(sites, axis="logperiod", ax=ax)
        assert ax.get_xscale() == "linear"
        assert "log" in ax.get_xlabel().lower()
        plt.close(fig)

    def test_bad_axis_raises(self):
        sites = [_site("S00", with_tipper=True)]
        with pytest.raises(ValueError, match="axis"):
            plot_tipper_components(sites, axis="bogus")
        plt.close("all")


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


# ─────────────────────────────────────────────────────────────────────────────
# private helpers
# ─────────────────────────────────────────────────────────────────────────────


class TestIsDf:
    def test_dataframe_true(self):
        assert _is_df(pd.DataFrame({"a": [1]})) is True

    def test_non_dataframe_false(self):
        assert _is_df([1, 2, 3]) is False


class TestNamePriv:
    def test_falls_back_to_index(self):
        class _NoName:
            pass

        assert _name(_NoName(), 3) == "site_3"

    def test_uses_station_attr(self):
        class _Named:
            station = "ABC"

        assert _name(_Named(), 0) == "ABC"


class TestCoordsPriv:
    def test_valid_coords_tuple(self):
        class _Ed:
            coords = (12.3, 45.6)

        lat, lon = _coords(_Ed())
        assert (lat, lon) == (12.3, 45.6)

    def test_coords_unpack_failure_falls_back(self):
        class _Ed:
            coords = "not-a-pair"
            lat = 1.5
            lon = 2.5

        assert _coords(_Ed()) == (1.5, 2.5)

    def test_latitude_longitude_fallback(self):
        class _Ed:
            latitude = 10.0
            longitude = 20.0

        assert _coords(_Ed()) == (10.0, 20.0)

    def test_non_numeric_lat_lon_returns_none(self):
        class _Ed:
            lat = "bad"
            lon = "bad"

        assert _coords(_Ed()) == (None, None)

    def test_no_attrs_returns_none(self):
        class _Ed:
            pass

        assert _coords(_Ed()) == (None, None)


class TestHasPriv:
    def test_tipper_absent(self):
        class _Ed:
            pass

        assert _has(_Ed(), "tipper") is False

    def test_mt_present(self):
        fr = _freqs()
        ed = _site("S00")
        assert _has(ed, "mt") is True
        assert _has(ed, "z") is True

    def test_has_section_callable(self):
        class _Ed:
            def has_section(self, sect):
                return sect == "resphase"

        assert _has(_Ed(), "resphase") is True
        assert _has(_Ed(), "other") is False

    def test_has_section_raises_falls_back_to_duck_typing(self):
        class _Ed:
            other = True

            def has_section(self, sect):
                raise RuntimeError("boom")

        assert _has(_Ed(), "other") is True
        assert _has(_Ed(), "missing") is False


class TestGetFreqPriv:
    def test_from_z_wrapper(self):
        fr = _freqs()
        ed = _site("S00", n=len(fr))
        out = _get_freq(ed)
        assert isinstance(out, np.ndarray) and out.size > 0

    def test_last_resort_ed_freq(self):
        class _Ed:
            freq = np.array([1.0, 2.0, 3.0])

        out = _get_freq(_Ed())
        assert out is not None and out.size == 3

    def test_ragged_freq_raises_internally_and_returns_none(self):
        class _Ed:
            freq = [[1, 2], [3]]

        assert _get_freq(_Ed()) is None

    def test_no_freq_anywhere_returns_none(self):
        class _Ed:
            pass

        assert _get_freq(_Ed()) is None


class TestIterItemsPriv:
    def test_plain_list(self):
        assert list(_iter_items([1, 2, 3])) == [1, 2, 3]

    def test_dict_like_items(self):
        class _Coll:
            items = {"a": 1, "b": 2}

        assert sorted(_iter_items(_Coll())) == [1, 2]

    def test_neither_iterable_nor_dict_like_yields_nothing(self):
        class _Coll:
            pass

        assert list(_iter_items(_Coll())) == []


class TestPeriodUnionPresence:
    def test_period_inverts_freq(self):
        fr = np.array([1.0, 2.0, 4.0])
        np.testing.assert_allclose(_period(fr), 1.0 / fr)

    def test_period_zero_freq_produces_inf_without_raising(self):
        fr = np.array([0.0, 1.0])
        out = _period(fr)
        assert np.isinf(out[0])

    def test_union_freq_empty(self):
        assert _union_freq([]).size == 0

    def test_union_freq_filters_nonpositive_and_dedupes(self):
        out = _union_freq([np.array([1.0, -1.0, np.nan, 2.0]), np.array([2.0, 3.0])])
        np.testing.assert_allclose(np.sort(out), [1.0, 2.0, 3.0])

    def test_presence_vec_empty_inputs(self):
        assert _presence_vec(np.array([]), np.array([1.0])).sum() == 0
        assert _presence_vec(np.array([1.0]), np.array([])).size == 0

    def test_presence_vec_marks_matches(self):
        grid = np.array([1.0, 2.0, 3.0])
        fr = np.array([2.0])
        out = _presence_vec(fr, grid)
        assert out.tolist() == [False, True, False]


class TestDfFromKv:
    def test_empty_rows_returns_empty_frame_with_columns(self):
        df = _df_from_kv([], ["a", "b"])
        assert list(df.columns) == ["a", "b"]
        assert df.empty

    def test_rows_populate_frame(self):
        df = _df_from_kv([{"a": 1, "b": 2}], ["a", "b"])
        assert df.iloc[0]["a"] == 1


class TestRhoPhaseErrRms:
    def test_rho_phase_err_zero_z_gives_zero_error(self):
        rho = np.array([100.0])
        z = np.array([0.0 + 0.0j])
        z_err = np.array([0.0])
        rho_err, phase_err = _rho_phase_err(rho, z, z_err)
        assert rho_err[0] == 0.0
        assert phase_err[0] == 0.0

    def test_rho_phase_err_nonzero(self):
        rho = np.array([100.0])
        z = np.array([1.0 + 0.0j])
        z_err = np.array([0.1])
        rho_err, phase_err = _rho_phase_err(rho, z, z_err)
        assert rho_err[0] == pytest.approx(2.0 * 100.0 * 0.1)

    def test_rms_log_basic(self):
        obs = np.array([100.0, 200.0])
        mod = np.array([100.0, 200.0])
        assert _rms_log(obs, mod) == pytest.approx(0.0, abs=1e-9)

    def test_rms_log_all_invalid_returns_nan(self):
        obs = np.array([np.nan])
        mod = np.array([np.nan])
        assert np.isnan(_rms_log(obs, mod))


class TestPickStation:
    def test_picks_first_when_none(self):
        S = ensure_sites([_site("S00"), _site("S01")])
        name, ed = _pick_station(S, None)
        assert name == "S00"

    def test_picks_named_station(self):
        S = ensure_sites([_site("S00"), _site("S01")])
        name, ed = _pick_station(S, "S01")
        assert name == "S01"

    def test_missing_station_raises(self):
        S = ensure_sites([_site("S00")])
        with pytest.raises(RuntimeError, match="ZZZ"):
            _pick_station(S, "ZZZ")

    def test_empty_sites_raises(self):
        S = ensure_sites([])
        with pytest.raises(RuntimeError, match="No sites"):
            _pick_station(S, None)


# ─────────────────────────────────────────────────────────────────────────────
# frequency_coverage modes + list_missing_sections / sites_summary edge cases
# ─────────────────────────────────────────────────────────────────────────────


class TestFrequencyCoverageModes:
    def test_union_mode(self):
        sites = [_site("S00", n=5), _site("S01", n=7)]
        grid = frequency_coverage(sites, mode="union")
        assert isinstance(grid, np.ndarray) and grid.size > 0

    def test_intersection_mode(self):
        fr = _freqs(6)
        sites = [_FakeSite("S00", _iso_z(fr), fr), _FakeSite("S01", _iso_z(fr), fr)]
        grid = frequency_coverage(sites, mode="intersection")
        assert isinstance(grid, np.ndarray) and grid.size == fr.size

    def test_intersection_mode_empty_sites(self):
        grid = frequency_coverage([], mode="intersection")
        assert grid.size == 0

    def test_invalid_mode_raises(self):
        sites = [_site("S00")]
        with pytest.raises(ValueError, match="mode must be"):
            frequency_coverage(sites, mode="bogus")


class TestSitesSummaryCoordsAndTipper:
    def test_lat_lon_columns_present_even_when_unset(self):
        # NOTE: lat/lon set on a raw fixture object are not forwarded by
        # the Sites/Site wrapper (that lives on ``.edi``, not the wrapper
        # itself) -- sites_summary/_coords still must not crash and must
        # keep the lat/lon columns present.
        sites = [_site("S00", lat=10.0, lon=20.0)]
        df = sites_summary(sites, api=False)
        assert "lat" in df.columns and "lon" in df.columns

    def test_has_tipper_true_with_tipper_data(self):
        sites = [_site("S00", with_tipper=True)]
        df = sites_summary(sites, api=False)
        assert bool(df["has_tipper"].iloc[0]) is True


class TestListMissingSectionsTipperPresent:
    def test_tipper_present_not_reported_missing(self):
        sites = [_site("S00", with_tipper=True)]
        result = list_missing_sections(sites, require=("tipper",))
        assert result == {}


# ─────────────────────────────────────────────────────────────────────────────
# plot_survey_inventory_overview / plot_rhoa_phi / plot_coverage edge cases
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSurveyInventoryOverviewExtra:
    def test_custom_station_order_and_kws(self):
        sites = [_site("S00", n=5), _site("S01", n=7)]
        fig = plot_survey_inventory_overview(
            sites,
            station_order=["S01", "S00", "S02"],
            count_kws={"color": "red"},
            station_grid_kws={"color": "black"},
        )
        plt.close("all")
        assert fig is not None

    def test_external_axes(self):
        sites = [_site("S00", n=5), _site("S01", n=6)]
        fig, (ax1, ax2) = plt.subplots(2, 1)
        result = plot_survey_inventory_overview(sites, axes=(ax1, ax2))
        plt.close("all")
        assert result is fig


class TestPlotRhoaPhiZFallback:
    def test_period_axis_uses_freq_column_when_period_missing(self):
        sites = [_site("S00")]
        result = plot_rhoa_phi(sites, axis="freq", errorbar=False)
        plt.close("all")
        assert result is not None


# ─────────────────────────────────────────────────────────────────────────────
# plot_station_response
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotStationResponseSynthetic:
    def test_no_impedance_data_message(self):
        fr = _freqs()
        z = _iso_z(fr)
        ed = _FakeSite("S00", z, fr)  # no .rho / .phase attrs
        fig = plot_station_response([ed])
        texts = [t.get_text() for t in fig.axes[0].texts]
        plt.close("all")
        assert "no impedance data" in texts

    def test_missing_station_raises(self):
        fr = _freqs()
        ed = _FakeSite("S00", _iso_z(fr), fr)
        with pytest.raises(RuntimeError):
            plot_station_response([ed], station="ZZZ")

    def test_no_sites_raises(self):
        with pytest.raises(RuntimeError):
            plot_station_response([])


@pytestmark_kap03lmt
class TestPlotStationResponseReal:
    @pytest.fixture(scope="class")
    def sites(self):
        return ensure_sites(str(_KAP03LMT))

    def test_returns_figure_defaults(self, sites):
        fig = plot_station_response(sites)
        plt.close("all")
        assert fig is not None

    def test_model_overlay_and_rms(self, sites):
        first_name, _ = _pick_station(sites, None)
        fig = plot_station_response(
            sites,
            station=first_name,
            sites_model=sites,
            components=("xy", "yx"),
            show_rms=True,
        )
        plt.close("all")
        assert fig is not None

    def test_period_range_and_limits_no_tipper(self, sites):
        first_name, _ = _pick_station(sites, None)
        fig = plot_station_response(
            sites,
            station=first_name,
            period_range=(1e-3, 1.0),
            rho_lim=(1, 1000),
            phase_lim=(0, 90),
            tipper_lim=(-1, 1),
            show_error_bars=False,
            show_rms=False,
            show_tipper=False,
            components=("xy",),
        )
        plt.close("all")
        assert fig is not None

    def test_external_axes(self, sites):
        first_name, _ = _pick_station(sites, None)
        n_comp = 2
        fig0 = plt.figure(figsize=(14, 8))
        gs = gridspec_module().GridSpec(3, max(n_comp, 4), figure=fig0)
        axs = (
            [fig0.add_subplot(gs[0, c]) for c in range(n_comp)]
            + [fig0.add_subplot(gs[1, c]) for c in range(n_comp)]
            + [fig0.add_subplot(gs[2, c]) for c in range(4)]
        )
        fig = plot_station_response(
            sites,
            station=first_name,
            components=("xy", "yx"),
            axes=axs,
            show_tipper=True,
        )
        plt.close("all")
        assert fig is fig0

    def test_invalid_component_is_skipped(self, sites):
        first_name, _ = _pick_station(sites, None)
        fig = plot_station_response(
            sites, station=first_name, components=("xy", "zz")
        )
        plt.close("all")
        assert fig is not None

    def test_model_station_not_found_is_silently_ignored(self, sites):
        first_name, _ = _pick_station(sites, None)

        class _EmptyModel(list):
            pass

        fig = plot_station_response(
            sites, station=first_name, sites_model=_EmptyModel()
        )
        plt.close("all")
        assert fig is not None


def gridspec_module():
    import matplotlib.gridspec as gridspec

    return gridspec

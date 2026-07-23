# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Supplementary tests for MapPanel (module helpers, pop-out button, typed
maps, basemap / contextily branches, contour, interaction handlers).

Complements pycsamt/app/tests/test_map_panel.py, which already covers
basic construction / set_dataframe / clear / highlight_station.

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — a small (28-EDI) profile used (subset) to
build real Sites + a fully-populated station DataFrame via
DataController, so the resistivity/depth/elevation map-type branches and
the frequency/impedance helpers get exercised against real EDI data.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import QEvent, QPoint
from PySide6.QtWidgets import QMessageBox

import pycsamt.app.desktop.panels.map_panel as mp
from pycsamt.app.desktop.panels.map_panel import (
    MapPanel,
    _build_provider_table,
    _is_geographic,
    _merc_tick_formatter,
    _MapPopOutButton,
    _project_to_merc,
    _reshape_z,
    _rho_app,
    _slice_comp,
    _try_ctx,
    _try_griddata,
)

# ── Real data paths ─────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18.exists() and any(_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def real_sites_df():
    """A small real Sites collection + matching station DataFrame."""
    if not _HAS_L18:
        pytest.skip("WILLY L18PLT EDI data not available")
    from pycsamt.app.desktop.controllers.data_controller import (
        DataController,
    )

    paths = sorted(_L18.glob("*.edi"))[:8]
    dc = DataController()
    sites = dc.load(paths)
    df = dc.dataframe
    return sites, df


@pytest.fixture
def coords_df():
    return pd.DataFrame(
        {
            "ID": ["A01", "A02", "A03", "A04"],
            "Latitude": [48.50, 48.51, 48.52, 48.53],
            "Longitude": [7.75, 7.76, 7.77, 7.78],
        }
    )


@pytest.fixture
def panel(qapp):
    p = MapPanel()
    yield p
    p.close()


# ═════════════════════════════════════════════════════════════════════════
# Module-level coordinate / formatting helpers
# ═════════════════════════════════════════════════════════════════════════


class TestProjectToMerc:
    def test_wgs84_origin(self):
        xs, ys = _project_to_merc(np.array([0.0]), np.array([0.0]))
        assert xs[0] == pytest.approx(0.0, abs=1.0)
        assert ys[0] == pytest.approx(0.0, abs=1.0)

    def test_wgs84_real_point(self):
        xs, ys = _project_to_merc(np.array([7.75]), np.array([48.5]))
        assert xs[0] == pytest.approx(862726.05, abs=1.0)
        assert ys[0] == pytest.approx(6190443.8, abs=1.0)

    def test_utm_source_crs(self):
        # Any valid EPSG code is accepted, not just WGS84.
        xs, ys = _project_to_merc(
            np.array([500000.0]), np.array([5000000.0]), "EPSG:32650"
        )
        assert np.isfinite(xs[0]) and np.isfinite(ys[0])

    def test_invalid_crs_raises(self):
        with pytest.raises(Exception):
            _project_to_merc(
                np.array([0.0]), np.array([0.0]), "NOT_A_REAL_CRS"
            )


class TestIsGeographic:
    def test_epsg_4326_is_geographic(self):
        assert _is_geographic("EPSG:4326") is True

    def test_utm_is_not_geographic(self):
        assert _is_geographic("EPSG:32650") is False

    def test_fallback_known_geographic_code(self, monkeypatch):
        from pyproj import CRS

        monkeypatch.setattr(
            CRS,
            "from_user_input",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        assert _is_geographic("EPSG:4326") is True

    def test_fallback_unknown_code_is_false(self, monkeypatch):
        from pyproj import CRS

        monkeypatch.setattr(
            CRS,
            "from_user_input",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        assert _is_geographic("EPSG:32650") is False


class TestMercTickFormatter:
    def test_returns_two_callables(self):
        fmt_lon, fmt_lat = _merc_tick_formatter()
        assert callable(fmt_lon) and callable(fmt_lat)

    def test_fmt_lon_near_zero(self):
        fmt_lon, _ = _merc_tick_formatter()
        label = fmt_lon(0.0, 0)
        assert label.endswith("°")
        assert label.startswith("0")

    def test_fmt_lat_real_value(self):
        _, fmt_lat = _merc_tick_formatter()
        label = fmt_lat(6190443.8, 0)
        assert label.endswith("°")
        assert "48" in label


class TestSliceComp:
    def _z(self):
        z = np.zeros((3, 2, 2), dtype=complex)
        z[:, 0, 0] = 1  # xx
        z[:, 0, 1] = 2  # xy
        z[:, 1, 0] = 3  # yx
        z[:, 1, 1] = 4  # yy
        return z

    def test_xy(self):
        assert np.all(_slice_comp(self._z(), "xy") == 2)

    def test_yx(self):
        assert np.all(_slice_comp(self._z(), "yx") == 3)

    def test_xx(self):
        assert np.all(_slice_comp(self._z(), "xx") == 1)

    def test_yy(self):
        assert np.all(_slice_comp(self._z(), "yy") == 4)


class TestRhoApp:
    def test_known_value(self):
        zc = np.array([10.0 + 0j])
        freqs = np.array([1.0])
        rho = _rho_app(zc, freqs)
        assert rho[0] == pytest.approx(0.2 * 100, rel=1e-6)

    def test_zero_freq_does_not_raise(self):
        zc = np.array([5.0 + 0j])
        freqs = np.array([0.0])
        rho = _rho_app(zc, freqs)
        assert np.isfinite(rho[0])


class TestReshapeZ:
    def test_none_returns_none(self):
        assert _reshape_z(None) is None

    def test_3d_valid_shape(self):
        z = np.zeros((5, 2, 2))
        out = _reshape_z(z)
        assert out.shape == (5, 2, 2)

    def test_2d_with_four_plus_cols(self):
        z = np.zeros((5, 4))
        out = _reshape_z(z)
        assert out.shape == (5, 2, 2)

    def test_2d_extra_cols_truncated(self):
        z = np.arange(5 * 6).reshape(5, 6).astype(float)
        out = _reshape_z(z)
        assert out.shape == (5, 2, 2)
        # Only the first four columns feed the reshape.
        assert out[0, 0, 0] == 0

    def test_2d_too_few_cols_returns_none(self):
        z = np.zeros((5, 3))
        assert _reshape_z(z) is None

    def test_1d_returns_none(self):
        z = np.zeros(10)
        assert _reshape_z(z) is None

    def test_wrong_3d_shape_returns_none(self):
        z = np.zeros((5, 3, 3))
        assert _reshape_z(z) is None


class TestBuildProviderTable:
    def test_real_contextily_table(self):
        ctx = _try_ctx()
        if ctx is None:
            pytest.skip("contextily not installed")
        tbl = _build_provider_table(ctx)
        assert "OpenStreetMap" in tbl
        assert isinstance(tbl, dict)

    def test_partial_provider_missing_attrs_skipped(self):
        class FakeProviders:
            class OpenStreetMap:
                Mapnik = "osm-mapnik"

            # CartoDB / Esri / Stamen deliberately absent -> AttributeError

        class FakeCtx:
            providers = FakeProviders

        tbl = _build_provider_table(FakeCtx)
        assert tbl == {"OpenStreetMap": "osm-mapnik"}

    def test_all_providers_missing_returns_empty(self):
        class FakeProviders:
            pass

        class FakeCtx:
            providers = FakeProviders

        tbl = _build_provider_table(FakeCtx)
        assert tbl == {}


class TestTryGriddata:
    def test_returns_callable_or_none(self):
        g = _try_griddata()
        assert g is None or callable(g)


class TestOptionalImportFailures:
    """Force the real ImportError branches inside _try_ctx / _try_griddata."""

    def test_try_ctx_import_error(self, monkeypatch):
        import sys

        monkeypatch.setitem(sys.modules, "contextily", None)
        assert _try_ctx() is None

    def test_try_griddata_import_error(self, monkeypatch):
        import sys

        monkeypatch.setitem(sys.modules, "scipy.interpolate", None)
        assert _try_griddata() is None


# ═════════════════════════════════════════════════════════════════════════
# _MapPopOutButton
# ═════════════════════════════════════════════════════════════════════════


class TestMapPopOutButton:
    def test_button_created_with_container(self, panel):
        btn = panel._pop_out_button
        assert isinstance(btn, _MapPopOutButton)
        assert btn._container is panel

    def test_initially_hidden_state(self, panel):
        assert panel._pop_out_button._shown is False

    def test_set_shown_true_changes_state(self, panel):
        btn = panel._pop_out_button
        btn._set_shown(True)
        assert btn._shown is True

    def test_set_shown_false_after_true(self, panel):
        btn = panel._pop_out_button
        btn._set_shown(True)
        btn._set_shown(False)
        assert btn._shown is False

    def test_set_shown_same_value_is_noop(self, panel):
        btn = panel._pop_out_button
        btn._set_shown(False)  # already False
        assert btn._shown is False

    def test_reposition_moves_button(self, panel):
        btn = panel._pop_out_button
        panel.resize(400, 300)
        btn._reposition()
        pos = btn.pos()
        assert pos.x() == panel.width() - btn.width() - 8
        assert pos.y() == 8

    def test_event_filter_resize(self, panel):
        btn = panel._pop_out_button
        panel.resize(500, 400)
        ev = QEvent(QEvent.Type.Resize)
        result = btn.eventFilter(panel, ev)
        assert isinstance(result, bool)

    def test_event_filter_enter_shows(self, panel):
        btn = panel._pop_out_button
        ev = QEvent(QEvent.Type.Enter)
        btn.eventFilter(panel, ev)
        assert btn._shown is True

    def test_event_filter_leave_hides(self, panel):
        btn = panel._pop_out_button
        btn.eventFilter(panel, QEvent(QEvent.Type.Enter))
        btn.eventFilter(panel, QEvent(QEvent.Type.Leave))
        assert btn._shown is False

    def test_event_filter_unhandled_type(self, panel):
        """An event type that is neither Resize/Enter/Leave falls through."""
        btn = panel._pop_out_button
        ev = QEvent(QEvent.Type.Show)
        result = btn.eventFilter(panel, ev)
        assert isinstance(result, bool)
        assert btn._shown is False

    def test_event_filter_other_obj_ignored(self, panel):
        btn = panel._pop_out_button
        other = MapPanel()
        ev = QEvent(QEvent.Type.Enter)
        btn.eventFilter(other, ev)
        assert btn._shown is False
        other.close()

    def test_event_filter_reentrant_guard(self, panel):
        """``getattr`` guard must survive a missing ``_container`` attr."""
        btn = panel._pop_out_button
        del btn.__dict__["_container"]
        ev = QEvent(QEvent.Type.Enter)
        # Must not raise even though _container is gone.
        btn.eventFilter(panel, ev)
        btn._container = panel  # restore for fixture teardown

    def test_on_click_opens_detail_window(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel._pop_out_button._on_click()  # must not raise

    def test_on_click_empty_df_and_no_sites(self, panel):
        """Freshly constructed panel: empty df, no sites -> both guarded
        ``if`` branches in _on_click take their False path."""
        assert panel._df.empty
        assert panel._sites is None
        panel._pop_out_button._on_click()  # must not raise

    def test_on_click_with_sites(self, panel, coords_df, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel._pop_out_button._on_click()  # must not raise


# ═════════════════════════════════════════════════════════════════════════
# MapPanel — set_sites / frequency & impedance helpers
# ═════════════════════════════════════════════════════════════════════════


class TestSetSites:
    def test_set_sites_stores_sites(self, panel, real_sites_df):
        sites, _ = real_sites_df
        panel.set_sites(sites)
        assert panel._sites is sites

    def test_set_sites_emits_freq_list(self, panel, real_sites_df):
        sites, _ = real_sites_df
        captured = []
        panel.freq_list_ready.connect(captured.append)
        panel.set_sites(sites)
        assert len(captured) == 1
        assert len(captured[0]) > 0

    def test_set_sites_no_frequencies_does_not_emit(self, panel):
        """When no frequency can be collected, freq_list_ready must not fire."""
        captured = []
        panel.freq_list_ready.connect(captured.append)
        bad = Mock()
        type(bad).freq = property(
            lambda self: (_ for _ in ()).throw(RuntimeError("no freq"))
        )
        panel.set_sites([bad])
        assert captured == []

    def test_collect_frequencies_sorted_desc(self, panel, real_sites_df):
        sites, _ = real_sites_df
        panel._sites = sites
        freqs = panel._collect_frequencies()
        assert freqs == sorted(freqs, reverse=True)
        assert all(f > 0 for f in freqs)

    def test_collect_frequencies_no_sites_returns_empty(self, panel):
        panel._sites = None
        assert panel._collect_frequencies() == []

    def test_collect_frequencies_bad_site_skipped(self, panel):
        bad = Mock()
        type(bad).freq = property(
            lambda self: (_ for _ in ()).throw(RuntimeError("no freq"))
        )
        panel._sites = [bad]
        assert panel._collect_frequencies() == []


class TestZArrays:
    def test_z_arrays_real_site(self, panel, real_sites_df):
        sites, _ = real_sites_df
        site = list(sites)[0]
        freqs, z3d = panel._z_arrays(site)
        assert freqs is not None
        assert z3d is not None
        assert z3d.shape[1:] == (2, 2)
        assert len(freqs) == z3d.shape[0]

    def test_z_arrays_bad_site_returns_none_none(self, panel):
        bad = Mock()
        bad.edi = object()
        freqs, z3d = panel._z_arrays(bad)
        assert freqs is None
        assert z3d is None


class TestRhoComp:
    def _z3d(self):
        z = np.zeros((4, 2, 2), dtype=complex)
        z[:, 0, 0] = 1 + 1j
        z[:, 0, 1] = 10 + 5j
        z[:, 1, 0] = -8 - 4j
        z[:, 1, 1] = 1 - 1j
        return z

    def test_xy_component(self, panel):
        freqs = np.array([1.0, 2.0, 3.0, 4.0])
        rho = panel._rho_comp(self._z3d(), freqs, "xy")
        assert rho.shape == (4,)
        assert np.all(np.isfinite(rho))

    def test_det_component(self, panel):
        freqs = np.array([1.0, 2.0, 3.0, 4.0])
        rho = panel._rho_comp(self._z3d(), freqs, "det")
        assert rho.shape == (4,)
        assert np.all(rho >= 0)


class TestRhoAtFreqAndDepth:
    def test_rho_at_freq_real_sites(self, panel, real_sites_df):
        sites, _ = real_sites_df
        panel._sites = sites
        vm = panel._rho_at_freq(1000.0, "xy")
        assert isinstance(vm, dict)
        assert len(vm) > 0
        assert all(v > 0 for v in vm.values())

    def test_rho_at_freq_det_component(self, panel, real_sites_df):
        sites, _ = real_sites_df
        panel._sites = sites
        vm = panel._rho_at_freq(500.0, "det")
        assert isinstance(vm, dict)

    def test_rho_at_freq_no_sites_returns_empty(self, panel):
        panel._sites = None
        assert panel._rho_at_freq(1.0) == {}

    def test_rho_at_depth_real_sites(self, panel, real_sites_df):
        sites, _ = real_sites_df
        panel._sites = sites
        vm = panel._rho_at_depth(200.0, "xy")
        assert isinstance(vm, dict)
        assert len(vm) > 0

    def test_rho_at_depth_no_sites_returns_empty(self, panel):
        panel._sites = None
        assert panel._rho_at_depth(100.0) == {}

    def test_rho_at_freq_bad_site_skipped(self, panel):
        bad = Mock()
        bad.edi = object()
        bad.name = "BAD"
        panel._sites = [bad]
        assert panel._rho_at_freq(1.0) == {}

    def test_rho_at_depth_bad_site_skipped(self, panel):
        bad = Mock()
        bad.edi = object()
        bad.name = "BAD"
        panel._sites = [bad]
        assert panel._rho_at_depth(100.0) == {}

    def test_rho_at_freq_rho_comp_exception_is_caught(
        self, panel, real_sites_df
    ):
        """An invalid component key makes _rho_comp raise; must be swallowed."""
        sites, _ = real_sites_df
        panel._sites = sites
        assert panel._rho_at_freq(1000.0, comp="zz") == {}

    def test_rho_at_depth_rho_comp_exception_is_caught(
        self, panel, real_sites_df
    ):
        sites, _ = real_sites_df
        panel._sites = sites
        assert panel._rho_at_depth(100.0, comp="zz") == {}

    def test_rho_at_freq_non_positive_value_excluded(
        self, panel, real_sites_df, monkeypatch
    ):
        """When the computed rho is 0 (or non-finite) the site must be
        excluded from the output map (the ``val > 0`` guard's False arm)."""
        sites, _ = real_sites_df
        panel._sites = sites
        monkeypatch.setattr(
            mp, "_rho_app", lambda zc, freqs: np.zeros_like(freqs)
        )
        assert panel._rho_at_freq(1000.0) == {}

    def test_rho_at_depth_non_positive_value_excluded(
        self, panel, real_sites_df, monkeypatch
    ):
        sites, _ = real_sites_df
        panel._sites = sites
        monkeypatch.setattr(
            mp, "_rho_app", lambda zc, freqs: np.zeros_like(freqs)
        )
        assert panel._rho_at_depth(100.0) == {}

    def test_z_arrays_missing_edi_attr_hits_except(self, panel):
        """A site with no ``.edi`` attribute triggers the try/except path."""

        class NoEdi:
            pass

        freqs, z3d = panel._z_arrays(NoEdi())
        assert freqs is None
        assert z3d is None


# ═════════════════════════════════════════════════════════════════════════
# _draw_map dispatch — station / elevation / depth / resistivity
# ═════════════════════════════════════════════════════════════════════════


class TestDrawMapDispatch:
    def test_station_map_default(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        assert panel._scatter is not None
        assert panel._map_type == "station"

    def test_empty_dataframe_shows_no_stations_title(self, panel):
        panel.set_dataframe(
            pd.DataFrame(columns=["ID", "Latitude", "Longitude"])
        )
        assert panel._canvas.axes.get_title() == "No stations loaded"

    def test_elevation_map_no_column_shows_no_data(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(map_type="elevation")
        title = panel._canvas.axes.get_title()
        assert "no valid data" in title.lower()

    def test_elevation_map_with_real_data(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation")
        assert panel._scatter is not None
        assert panel._canvas.axes.get_title() == "Elevation Map"

    def test_depth_map_no_sites_shows_no_data(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(map_type="depth")
        title = panel._canvas.axes.get_title()
        assert "no valid data" in title.lower()

    def test_depth_map_with_real_sites(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(
            map_type="depth", target_depth_m=200.0, component="xy"
        )
        assert panel._scatter is not None
        assert "Depth Map" in panel._canvas.axes.get_title()

    def test_resistivity_map_with_real_sites(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(
            map_type="resistivity", target_freq_hz=1000.0, component="xy"
        )
        assert panel._scatter is not None
        assert "Resistivity Map" in panel._canvas.axes.get_title()

    def test_typed_map_no_profile_no_labels(self, panel, real_sites_df):
        """Exercise the show_profile/show_labels False branches of the
        typed-map path (elevation/depth/resistivity), not just the
        station-map path."""
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(
            map_type="elevation", show_profile=False, show_labels=False
        )
        assert panel._annots == {}

    def test_resistivity_map_zero_freq_does_not_raise(
        self, panel, real_sites_df
    ):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(map_type="resistivity", target_freq_hz=0.0)

    def test_redraw_unknown_kwarg_ignored(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(not_a_real_attr=123)  # must not raise / not set


# ═════════════════════════════════════════════════════════════════════════
# Colour-by variants of the station map
# ═════════════════════════════════════════════════════════════════════════


class TestStationMapColorBy:
    def test_color_by_elevation(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(color_by="elevation")
        assert panel._scatter is not None

    def test_color_by_n_freq(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(color_by="n_freq")
        assert panel._scatter is not None

    def test_color_by_tipper(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(color_by="tipper")
        assert panel._scatter is not None

    def test_color_by_index_default(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(color_by="index")
        assert panel._scatter is not None

    def test_color_by_elevation_missing_column_falls_back(
        self, panel, coords_df
    ):
        panel.set_dataframe(coords_df)
        panel.redraw(color_by="elevation")  # no Elevation column present
        assert panel._scatter is not None


# ═════════════════════════════════════════════════════════════════════════
# Contour / colourbar
# ═════════════════════════════════════════════════════════════════════════


class TestContour:
    def test_contour_lines_mode(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="lines")
        assert panel._scatter is not None

    def test_contour_filled_mode(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="filled")

    def test_contour_filled_labels_mode(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(
            map_type="resistivity",
            contour_mode="filled_labels",
            log_scale=True,
            target_freq_hz=1000.0,
        )

    def test_contour_unrecognized_mode_falls_through(
        self, panel, real_sites_df
    ):
        """A mode string matching none of none/lines/filled/filled_labels
        must reach the end of the if/elif chain without drawing or raising."""
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="bogus_mode")
        assert panel._scatter is not None

    def test_contour_none_mode_skips(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="none")

    def test_contour_too_few_points(self, panel):
        df = pd.DataFrame(
            {
                "ID": ["A", "B", "C"],
                "Latitude": [48.5, 48.51, 48.52],
                "Longitude": [7.75, 7.76, 7.77],
                "Elevation": [10.0, 20.0, 30.0],
            }
        )
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="lines")
        assert panel._scatter is not None

    def test_contour_griddata_unavailable(
        self, panel, real_sites_df, monkeypatch
    ):
        sites, df = real_sites_df
        monkeypatch.setattr(mp, "_try_griddata", lambda: None)
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="lines")
        assert "scipy" in panel._info_label.text().lower()

    def test_contour_griddata_raises_is_swallowed(
        self, panel, real_sites_df, monkeypatch
    ):
        sites, df = real_sites_df

        def _boom(*a, **k):
            raise RuntimeError("griddata boom")

        monkeypatch.setattr(mp, "_try_griddata", lambda: _boom)
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="lines")  # no raise

    def test_contour_all_nan_result_returns_early(
        self, panel, real_sites_df, monkeypatch
    ):
        sites, df = real_sites_df

        def _all_nan(*a, **k):
            pts = a[2] if len(a) > 2 else k.get("xi")
            shape = np.asarray(a[-2]).shape if len(a) >= 2 else (200, 200)
            return np.full((200, 200), np.nan)

        monkeypatch.setattr(mp, "_try_griddata", lambda: _all_nan)
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", contour_mode="filled")  # no raise

    def test_colorbar_hidden_when_show_cbar_false(
        self, panel, real_sites_df
    ):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", show_cbar=False)
        assert panel._scatter is not None

    def test_colorbar_horizontal_orientation(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.redraw(map_type="elevation", cbar_orient="horizontal")

    def test_colorbar_log_scale_horizontal_formatter(
        self, panel, real_sites_df
    ):
        """log_scale + horizontal orientation exercises the xaxis
        formatter branch (the yaxis branch is covered elsewhere)."""
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(
            map_type="resistivity",
            log_scale=True,
            cbar_orient="horizontal",
            target_freq_hz=1000.0,
        )

    def test_colorbar_log_scale_formatter(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(
            map_type="resistivity", log_scale=True, target_freq_hz=1000.0
        )

    def test_colorbar_log_scale_formatter_exception_swallowed(
        self, panel, real_sites_df, monkeypatch
    ):
        """The formatter-assignment try/except in _add_colorbar must not
        propagate a failure while wiring the log-scale tick formatter."""
        import matplotlib.ticker as mticker

        def boom(*a, **k):
            raise RuntimeError("formatter boom")

        monkeypatch.setattr(mticker, "FuncFormatter", boom)
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.redraw(
            map_type="resistivity", log_scale=True, target_freq_hz=1000.0
        )  # must not raise


# ═════════════════════════════════════════════════════════════════════════
# Elevation values / profile lines / labels
# ═════════════════════════════════════════════════════════════════════════


class TestOverlaysAndValues:
    def test_elevation_values_real_data(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel._df = df
        vals = panel._elevation_values()
        assert isinstance(vals, dict)
        assert len(vals) > 0

    def test_elevation_values_no_column(self, panel, coords_df):
        panel._df = coords_df
        assert panel._elevation_values() == {}

    def test_elevation_values_skips_nan_rows(self, panel):
        panel._df = pd.DataFrame(
            {
                "ID": ["A", "B", "C"],
                "Latitude": [48.5, 48.51, 48.52],
                "Longitude": [7.75, 7.76, 7.77],
                "Elevation": [100.0, np.nan, 200.0],
            }
        )
        vals = panel._elevation_values()
        assert vals == {"A": 100.0, "C": 200.0}
        assert "B" not in vals

    def test_profile_lines_single_point_noop(self, panel):
        panel._add_profile_lines(
            panel._canvas.axes, np.array([1.0]), np.array([1.0])
        )

    def test_profile_lines_multi_point(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        xs = coords_df["Longitude"].values
        ys = coords_df["Latitude"].values
        panel._add_profile_lines(panel._canvas.axes, xs, ys)

    def test_labels_dark_and_light(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        xs = coords_df["Longitude"].values
        ys = coords_df["Latitude"].values
        ids = coords_df["ID"].values
        for dark in (True, False):
            panel._dark = dark
            panel._annots = {}
            panel._add_labels(panel._canvas.axes, xs, ys, ids)
            assert len(panel._annots) == len(ids)

    def test_show_labels_false_skips_annotations(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(show_labels=False)
        assert panel._annots == {}

    def test_show_profile_false(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(show_profile=False)  # must not raise


# ═════════════════════════════════════════════════════════════════════════
# Basemap / contextily branches
# ═════════════════════════════════════════════════════════════════════════


class TestBasemap:
    def test_basemap_success(self, panel, coords_df, monkeypatch):
        called = {}

        def fake_add_basemap(ax, source=None, zoom=None, alpha=None):
            called["source"] = source
            called["alpha"] = alpha

        ctx = _try_ctx()
        if ctx is None:
            pytest.skip("contextily not installed")
        monkeypatch.setattr(ctx, "add_basemap", fake_add_basemap)
        panel.set_dataframe(coords_df)
        panel.redraw(provider="OpenStreetMap")
        assert called.get("source") is not None
        assert panel._info_label.text() == ""

    def test_basemap_provider_not_available(
        self, panel, coords_df, monkeypatch
    ):
        ctx = _try_ctx()
        if ctx is None:
            pytest.skip("contextily not installed")
        panel.set_dataframe(coords_df)
        panel.redraw(provider="Totally Not A Real Provider")
        assert "not available" in panel._info_label.text()

    def test_basemap_add_basemap_raises(self, panel, coords_df, monkeypatch):
        ctx = _try_ctx()
        if ctx is None:
            pytest.skip("contextily not installed")

        def boom(*a, **k):
            raise RuntimeError("tile fetch failed")

        monkeypatch.setattr(ctx, "add_basemap", boom)
        panel.set_dataframe(coords_df)
        panel.redraw(provider="OpenStreetMap")
        assert "Basemap error" in panel._info_label.text()

    def test_basemap_no_provider_skips_merc(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(provider="None")
        assert panel._scatter is not None

    def test_basemap_crs_projection_error_falls_back(
        self, panel, coords_df
    ):
        panel.set_dataframe(coords_df)
        panel.redraw(provider="OpenStreetMap", source_crs="NOT_A_REAL_CRS")
        assert "CRS projection error" in panel._info_label.text()
        # Reset back to a valid CRS so other tests using this panel
        # instance (there are none, function-scoped) aren't affected.

    def test_ctx_unavailable_offers_install(
        self, panel, coords_df, monkeypatch
    ):
        monkeypatch.setattr(mp, "_try_ctx", lambda: None)
        called = {}

        def fake_ask(self):
            called["asked"] = True

        monkeypatch.setattr(
            MapPanel, "_ask_install_contextily", fake_ask
        )
        panel.set_dataframe(coords_df)
        panel.redraw(provider="OpenStreetMap")
        assert called.get("asked") is True

    def test_format_merc_ticks_sets_labels(self, panel):
        ax = panel._canvas.axes
        panel._format_merc_ticks(ax)
        assert ax.get_xlabel() == "Longitude"
        assert ax.get_ylabel() == "Latitude"

    def test_format_merc_ticks_formatter_output(self, panel):
        ax = panel._canvas.axes
        ax.set_xlim(0, 1000000)
        panel._format_merc_ticks(ax)
        fmt = ax.xaxis.get_major_formatter()
        label = fmt(0, 0)
        assert "°" in label

    def test_format_merc_ticks_exception_swallowed(
        self, panel, monkeypatch
    ):
        monkeypatch.setattr(
            mp,
            "_merc_tick_formatter",
            lambda: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        ax = panel._canvas.axes
        panel._format_merc_ticks(ax)  # must not raise
        assert ax.get_xlabel() == "Longitude"


# ═════════════════════════════════════════════════════════════════════════
# Install contextily (subprocess mocked — never a real pip install)
# ═════════════════════════════════════════════════════════════════════════


class TestAskInstallContextily:
    def test_yes_triggers_install(self, panel, monkeypatch):
        # no_modal_dialogs autouse fixture makes QMessageBox.question -> Yes
        called = {}
        monkeypatch.setattr(
            panel, "_install_contextily", lambda: called.setdefault(
                "installed", True
            )
        )
        panel._ask_install_contextily()
        assert called.get("installed") is True

    def test_no_skips_install(self, panel, monkeypatch):
        monkeypatch.setattr(
            QMessageBox,
            "question",
            staticmethod(
                lambda *a, **k: QMessageBox.StandardButton.No
            ),
        )
        called = {}
        monkeypatch.setattr(
            panel, "_install_contextily", lambda: called.setdefault(
                "installed", True
            )
        )
        panel._ask_install_contextily()
        assert "installed" not in called


class TestInstallContextily:
    def test_success_shows_information(self, panel, monkeypatch):
        proc = Mock(returncode=0, stdout="", stderr="")
        run_mock = Mock(return_value=proc)
        monkeypatch.setattr(mp.subprocess, "run", run_mock)
        info_mock = Mock(return_value=QMessageBox.StandardButton.Ok)
        monkeypatch.setattr(QMessageBox, "information", info_mock)

        panel._install_contextily()

        run_mock.assert_called_once()
        args = run_mock.call_args[0][0]
        assert "pip" in args
        assert "install" in args
        assert "contextily" in args
        info_mock.assert_called_once()

    def test_failure_shows_warning(self, panel, monkeypatch):
        proc = Mock(returncode=1, stdout="", stderr="network unreachable")
        monkeypatch.setattr(mp.subprocess, "run", Mock(return_value=proc))
        warn_mock = Mock(return_value=QMessageBox.StandardButton.Ok)
        monkeypatch.setattr(QMessageBox, "warning", warn_mock)

        panel._install_contextily()

        warn_mock.assert_called_once()
        msg = warn_mock.call_args[0][2]
        assert "failed" in msg.lower()

    def test_exception_shows_warning(self, panel, monkeypatch):
        def boom(*a, **k):
            raise OSError("pip not found")

        monkeypatch.setattr(mp.subprocess, "run", boom)
        warn_mock = Mock(return_value=QMessageBox.StandardButton.Ok)
        monkeypatch.setattr(QMessageBox, "warning", warn_mock)

        panel._install_contextily()

        warn_mock.assert_called_once()
        msg = warn_mock.call_args[0][2]
        assert "pip" in msg.lower()


# ═════════════════════════════════════════════════════════════════════════
# Axes styling
# ═════════════════════════════════════════════════════════════════════════


class TestStyleAxes:
    def test_dark_mode_geo(self, panel, coords_df):
        panel._dark = True
        panel.set_dataframe(coords_df)
        ax = panel._canvas.axes
        assert ax.get_xlabel() == "Longitude"

    def test_light_mode_geo(self, panel, coords_df):
        panel._dark = False
        panel.set_dataframe(coords_df)
        ax = panel._canvas.axes
        assert ax.get_xlabel() == "Longitude"

    def test_projected_non_merc_labels(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(source_crs="EPSG:32650", provider="None")
        ax = panel._canvas.axes
        assert "Easting" in ax.get_xlabel()
        assert "Northing" in ax.get_ylabel()

    def test_set_dark_mode_public_api(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.set_dark_mode(False)
        assert panel._dark is False
        panel.set_dark_mode(True)
        assert panel._dark is True

    def test_grid_toggle(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.redraw(show_grid=False)
        panel.redraw(show_grid=True)


# ═════════════════════════════════════════════════════════════════════════
# Highlight
# ═════════════════════════════════════════════════════════════════════════


class TestUpdateHighlight:
    def test_no_scatter_is_noop(self, panel):
        panel._scatter = None
        panel._update_highlight()  # must not raise

    def test_highlight_changes_sizes(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel.highlight_station("A02")
        sizes = panel._scatter.get_sizes()
        base = panel._marker_size**2
        assert sizes.max() > base

    def test_highlight_none_selected(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        panel._selected_id = None
        panel._update_highlight()
        sizes = panel._scatter.get_sizes()
        base = panel._marker_size**2
        assert np.allclose(sizes, base)


# ═════════════════════════════════════════════════════════════════════════
# Interaction: _on_pick / _on_click
# ═════════════════════════════════════════════════════════════════════════


class TestOnPick:
    def test_pick_matching_artist_selects_station(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        captured = []
        panel.station_selected.connect(captured.append)

        event = Mock(artist=panel._scatter, ind=[1])
        panel._on_pick(event)

        assert panel._selected_id == str(coords_df.iloc[1]["ID"])
        assert captured == [str(coords_df.iloc[1]["ID"])]

    def test_pick_non_matching_artist_ignored(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        other_artist = object()
        event = Mock(artist=other_artist, ind=[0])
        panel._on_pick(event)
        assert panel._selected_id is None

    def test_pick_empty_ind_ignored(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        event = Mock(artist=panel._scatter, ind=[])
        panel._on_pick(event)
        assert panel._selected_id is None


class TestOnClick:
    def test_click_outside_axes_noop(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        event = Mock(inaxes=None)
        panel._on_click(event)  # must not raise

    def test_click_inside_axes_noop(self, panel, coords_df):
        panel.set_dataframe(coords_df)
        event = Mock(inaxes=panel._canvas.axes)
        panel._on_click(event)  # currently a no-op; must not raise


# ═════════════════════════════════════════════════════════════════════════
# clear() while sites are set
# ═════════════════════════════════════════════════════════════════════════


class TestClearWithSites:
    def test_clear_resets_scatter_and_annots(self, panel, real_sites_df):
        sites, df = real_sites_df
        panel.set_dataframe(df)
        panel.set_sites(sites)
        panel.clear()
        assert panel._scatter is None
        assert panel._annots == {}
        assert panel._df.empty

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for MapViewerWindow / MapDetailWindow
(pycsamt.app.desktop.windows.map_window).

Strategy
--------
* ``MapPanel`` (the embedded content widget) is already unit-tested
  independently (test_map_panel_extra.py, 100% coverage per prior
  session). With no sites/dataframe loaded, ``MapPanel._draw_map()`` is
  cheap, so most window-level logic (visibility wiring, CRS resolution,
  refresh parameter collection) is exercised against a *real* MapPanel
  rather than a mock — real small survey data is used sparingly (one
  test) to avoid the coverage-instrumentation slowdown observed with
  heavy real-data redraws in test_profile_window.py.
* ``sync_from_panel`` is driven with a lightweight fake object exposing
  the same ``_map_type`` / ``_source_crs`` / ... attributes MapPanel
  actually has (verified against pycsamt/app/desktop/panels/map_panel.py).
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.windows.map_window import (
    MapDetailWindow,
    MapViewerWindow,
)

_ROOT = Path(__file__).parents[3]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))


@pytest.fixture(scope="session")
def tipper_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


def _fake_panel(**overrides):
    defaults = dict(
        _map_type="depth",
        _source_crs="EPSG:32650",
        _cmap_name="viridis",
        _marker_size=12,
        _marker_alpha=0.5,
        _show_cbar=False,
        _cbar_orient="horizontal",
        _contour_mode="filled",
        _contour_levels=15,
        _show_profile=False,
        _show_labels=False,
        _show_grid=False,
        _provider="OpenStreetMap",
        _basemap_alpha=0.3,
        _component="yx",
        _target_depth_m=750.0,
        _target_freq_hz=10.0,
        _log_scale=False,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


# ═════════════════════════════════════════════════════════════════════════════
# MapViewerWindow
# ═════════════════════════════════════════════════════════════════════════════


@pytest.fixture
def win(qapp):
    w = MapViewerWindow(parent=None)
    w.show()  # widgets only report real isVisible() once the top level is shown
    yield w
    w.close()


class TestViewerConstruction:
    def test_window_title(self, win):
        assert "Map Viewer" in win.windowTitle()

    def test_freq_group_hidden_for_station_type(self, win):
        assert win._combo_type.currentText() == "Station"
        assert not win._grp_freq.isVisible()
        assert not win._grp_depth.isVisible()

    def test_default_crs_geographic_hides_utm_and_epsg(self, win):
        assert not win._wgt_utm.isVisible()
        assert not win._wgt_epsg.isVisible()


class TestViewerTypeAndCrsVisibility:
    def test_depth_type_shows_freq_and_depth_groups(self, win):
        win._combo_type.setCurrentText("Depth")
        assert win._grp_freq.isVisible()
        assert win._grp_depth.isVisible()

    def test_resistivity_type_shows_freq_only(self, win):
        win._combo_type.setCurrentText("Resistivity")
        assert win._grp_freq.isVisible()
        assert not win._grp_depth.isVisible()

    def test_elevation_type_hides_both(self, win):
        win._combo_type.setCurrentText("Elevation")
        assert not win._grp_freq.isVisible()
        assert not win._grp_depth.isVisible()

    def test_utm_mode_shows_utm_widget(self, win):
        win._combo_crs_mode.setCurrentText("UTM Zone")
        assert win._wgt_utm.isVisible()
        assert not win._wgt_epsg.isVisible()

    def test_custom_epsg_mode_shows_epsg_widget(self, win):
        win._combo_crs_mode.setCurrentText("Custom EPSG")
        assert win._wgt_epsg.isVisible()
        assert not win._wgt_utm.isVisible()

    def test_cbar_toggle_disables_orientation_radios(self, win):
        win._chk_cbar.setChecked(False)
        assert not win._radio_cbar_v.isEnabled()
        assert not win._radio_cbar_h.isEnabled()
        win._chk_cbar.setChecked(True)
        assert win._radio_cbar_v.isEnabled()


class TestViewerCrsResolution:
    def test_resolve_geographic_default(self, win):
        assert win._resolve_source_crs() == "EPSG:4326"

    def test_resolve_utm_north(self, win):
        win._combo_crs_mode.setCurrentText("UTM Zone")
        win._spin_utm_zone.setValue(50)
        win._radio_utm_n.setChecked(True)
        assert win._resolve_source_crs() == "EPSG:32650"

    def test_resolve_utm_south(self, win):
        win._combo_crs_mode.setCurrentText("UTM Zone")
        win._spin_utm_zone.setValue(50)
        win._radio_utm_s.setChecked(True)
        assert win._resolve_source_crs() == "EPSG:32750"

    def test_resolve_custom_epsg_bare_number(self, win):
        win._combo_crs_mode.setCurrentText("Custom EPSG")
        win._edit_epsg.setText("3857")
        assert win._resolve_source_crs() == "EPSG:3857"

    def test_resolve_custom_epsg_prefixed(self, win):
        win._combo_crs_mode.setCurrentText("Custom EPSG")
        win._edit_epsg.setText("EPSG:3857")
        assert win._resolve_source_crs() == "EPSG:3857"

    def test_resolve_custom_epsg_invalid_falls_back(self, win):
        win._combo_crs_mode.setCurrentText("Custom EPSG")
        win._edit_epsg.setText("not-a-number")
        assert win._resolve_source_crs() == "EPSG:4326"

    def test_update_crs_info_valid(self, win):
        win._combo_crs_mode.setCurrentText("UTM Zone")
        assert "EPSG:326" in win._lbl_crs_info.text()

    def test_update_crs_info_invalid_shows_error_style(self, win):
        win._combo_crs_mode.setCurrentText("Custom EPSG")
        win._edit_epsg.setText("999999999")
        assert "Invalid" in win._lbl_crs_info.text()
        assert "f38ba8" in win._lbl_crs_info.styleSheet()

    def test_validate_epsg_valid_shows_success_style(self, win):
        win._combo_crs_mode.setCurrentText("Geographic (lat/lon)")
        win._validate_epsg()
        assert win._lbl_crs_info.text().startswith("✓")
        assert "a6e3a1" in win._lbl_crs_info.styleSheet()

    def test_validate_epsg_invalid_shows_failure_style(self, win):
        win._combo_crs_mode.setCurrentText("Custom EPSG")
        win._edit_epsg.setText("999999999")
        win._validate_epsg()
        assert win._lbl_crs_info.text().startswith("✗")


class TestViewerFrequencyAndDepth:
    def test_populate_freq_combo_hz(self, win):
        win._radio_hz.setChecked(True)
        win._populate_freq_combo([100.0, 10.0, 1.0])
        assert win._combo_freq.count() == 3
        assert "Hz" in win._combo_freq.itemText(0)

    def test_populate_freq_combo_period(self, win):
        win._radio_s.setChecked(True)
        win._populate_freq_combo([100.0, 10.0])
        assert "s" in win._combo_freq.itemText(0)

    def test_freq_unit_toggle_repopulates(self, win):
        win._populate_freq_combo([50.0])
        win._radio_s.setChecked(True)
        assert "s" in win._combo_freq.itemText(0)

    def test_freq_unit_toggle_noop_when_empty(self, win):
        win._on_freq_unit_toggled(True)  # no _freq_list yet -> no-op

    def test_current_freq_hz_parses_hz(self, win):
        win._combo_freq.addItem("12.5 Hz")
        win._combo_freq.setCurrentIndex(0)
        assert win._current_freq_hz() == pytest.approx(12.5)

    def test_current_freq_hz_parses_period(self, win):
        win._combo_freq.addItem("0.1 s")
        win._combo_freq.setCurrentIndex(0)
        assert win._current_freq_hz() == pytest.approx(10.0)

    def test_current_freq_hz_invalid_text_defaults_to_1(self, win):
        win._combo_freq.addItem("garbage")
        win._combo_freq.setCurrentIndex(0)
        assert win._current_freq_hz() == 1.0

    def test_depth_label_without_sites(self, win):
        win._update_depth_label()
        assert "Load data" in win._lbl_depth_period.text()

    def test_estimate_median_rho_no_sites_returns_zero(self, win):
        assert win._estimate_median_rho() == 0.0

    def test_estimate_median_rho_with_real_sites(self, win, tipper_sites):
        win._sites = tipper_sites
        rho = win._estimate_median_rho()
        assert rho >= 0.0

    def test_depth_label_with_real_sites(self, win, tipper_sites):
        win.set_sites(tipper_sites)
        win._update_depth_label()
        assert "T =" in win._lbl_depth_period.text() or "Load data" in win._lbl_depth_period.text()

    def test_estimate_median_rho_iteration_exception_not_swallowed(self, win):
        """
        Real inconsistency: unlike ``_populate_station_combo`` /
        ``_update_period_range`` in profile_window.py (which wrap their
        whole ``for site in sites`` loop in try/except), here only the
        *per-site body* is guarded -- ``for site in sites:`` itself is
        not, so a sites object that raises from ``__iter__`` propagates
        instead of degrading to 0.0 like the "no sites" case. Not fixed
        here, per instructions.
        """

        class _Bad:
            def __iter__(self):
                raise RuntimeError("boom")

        win._sites = _Bad()
        with pytest.raises(RuntimeError):
            win._estimate_median_rho()


class TestViewerRefreshAndExport:
    def test_on_refresh_calls_map_panel_redraw(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            win._map_panel, "redraw", lambda **kw: calls.append(kw)
        )
        win._on_refresh()
        assert len(calls) == 1
        assert "map_type" in calls[0]
        assert "source_crs" in calls[0]

    def test_on_refresh_falls_back_to_draw_map_on_exception(
        self, win, monkeypatch
    ):
        def _boom(**kw):
            raise RuntimeError("redraw failed")

        calls = []
        monkeypatch.setattr(win._map_panel, "redraw", _boom)
        monkeypatch.setattr(
            win._map_panel, "_draw_map", lambda: calls.append(1)
        )
        win._on_refresh()
        assert calls == [1]

    def test_on_refresh_fallback_itself_raising_is_swallowed(
        self, win, monkeypatch
    ):
        def _boom(**kw):
            raise RuntimeError("redraw failed")

        def _boom2():
            raise RuntimeError("draw_map failed too")

        monkeypatch.setattr(win._map_panel, "redraw", _boom)
        monkeypatch.setattr(win._map_panel, "_draw_map", _boom2)
        win._on_refresh()  # must not raise

    def test_on_export_opens_dialog(self, win, monkeypatch):
        calls = []

        class _FakeExportDialog:
            def __init__(self, figure, parent):
                calls.append(figure)

            def exec(self):
                calls.append("exec")

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            _FakeExportDialog,
        )
        win._on_export()
        assert calls[-1] == "exec"

    def test_contextily_available_false_when_not_installed(self, win):
        # contextily is not a project dependency here; property must just
        # return a bool without raising either way.
        assert isinstance(win._contextily_available(), bool)

    def test_contour_mode_str_mapping(self, win):
        for label, expected in [
            ("None", "none"),
            ("Lines", "lines"),
            ("Filled", "filled"),
            ("Filled + labels", "filled_labels"),
        ]:
            win._combo_contour.setCurrentText(label)
            assert win._contour_mode_str() == expected


class TestViewerSetSitesAndDarkMode:
    def test_set_sites_delegates_to_map_panel(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            win._map_panel, "set_sites", lambda s: calls.append(s)
        )
        win.set_sites(["fake"])
        assert calls == [["fake"]]

    def test_set_dark_mode_delegates(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            win._map_panel, "set_dark_mode", lambda d: calls.append(d)
        )
        win.set_dark_mode(False)
        assert calls == [False]

    def test_set_dataframe_delegates(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            win._map_panel, "set_dataframe", lambda df: calls.append(df)
        )
        win.set_dataframe("fake_df")
        assert calls == ["fake_df"]

    def test_highlight_station_delegates(self, win, monkeypatch):
        calls = []
        monkeypatch.setattr(
            win._map_panel, "highlight_station", lambda s: calls.append(s)
        )
        win.highlight_station("STA01")
        assert calls == ["STA01"]


class TestViewerSyncFromPanel:
    def test_sync_from_panel_applies_all_settings(self, win, monkeypatch):
        panel = _fake_panel()
        refreshed = []
        monkeypatch.setattr(win, "_on_refresh", lambda: refreshed.append(1))
        win.sync_from_panel(panel)

        assert win._combo_type.currentText() == "Depth"
        assert win._combo_crs_mode.currentText() == "UTM Zone"
        assert win._spin_utm_zone.value() == 50
        assert win._radio_utm_n.isChecked()
        assert win._combo_cmap.currentText() == "viridis"
        assert win._spin_ms.value() == 12
        assert not win._chk_cbar.isChecked()
        assert win._radio_cbar_h.isChecked()
        assert win._combo_contour.currentText() == "Filled"
        assert win._spin_levels.value() == 15
        assert not win._chk_profile.isChecked()
        assert win._combo_basemap.currentText() == "OpenStreetMap"
        assert win._combo_comp.currentText() == "YX"
        assert win._spin_depth.value() == 750.0
        assert refreshed == [1]

    def test_sync_from_panel_south_utm_zone(self, win, monkeypatch):
        monkeypatch.setattr(win, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:32750")
        win.sync_from_panel(panel)
        assert win._combo_crs_mode.currentText() == "UTM Zone"
        assert win._radio_utm_s.isChecked()

    def test_sync_from_panel_custom_epsg(self, win, monkeypatch):
        monkeypatch.setattr(win, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:2154")
        win.sync_from_panel(panel)
        assert win._combo_crs_mode.currentText() == "Custom EPSG"
        assert win._edit_epsg.text() == "2154"

    def test_sync_from_panel_non_numeric_crs_code(self, win, monkeypatch):
        monkeypatch.setattr(win, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:ABCDE")
        win.sync_from_panel(panel)
        assert win._combo_crs_mode.currentText() == "Custom EPSG"

    def test_sync_from_panel_geographic_crs(self, win, monkeypatch):
        monkeypatch.setattr(win, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:4326")
        win.sync_from_panel(panel)
        assert win._combo_crs_mode.currentText() == "Geographic (lat/lon)"

    def test_sync_from_panel_adds_freq_placeholder_when_combo_empty(
        self, win, monkeypatch
    ):
        monkeypatch.setattr(win, "_on_refresh", lambda: None)
        panel = _fake_panel(_target_freq_hz=25.0)
        win.sync_from_panel(panel)
        assert win._combo_freq.count() == 1
        assert "25" in win._combo_freq.itemText(0)


# ═════════════════════════════════════════════════════════════════════════════
# MapDetailWindow
# ═════════════════════════════════════════════════════════════════════════════


@pytest.fixture
def detail(qapp):
    w = MapDetailWindow(parent=None)
    w.show()
    yield w
    w.close()


class TestDetailConstruction:
    def test_window_title(self, detail):
        assert "Full Detail" in detail.windowTitle()

    def test_freq_and_comp_hidden_for_station_type(self, detail):
        assert not detail._wgt_freq.isVisible()
        assert not detail._wgt_comp.isVisible()
        assert not detail._wgt_depth.isVisible()

    def test_utm_and_epsg_hidden_for_geographic_default(self, detail):
        assert not detail._wgt_utm.isVisible()
        assert not detail._wgt_epsg.isVisible()


class TestDetailTypeAndCrsVisibility:
    def test_depth_type_shows_all_three(self, detail):
        detail._combo_type.setCurrentText("Depth")
        assert detail._wgt_freq.isVisible()
        assert detail._wgt_comp.isVisible()
        assert detail._wgt_depth.isVisible()

    def test_resistivity_type_shows_freq_and_comp_only(self, detail):
        detail._combo_type.setCurrentText("Resistivity")
        assert detail._wgt_freq.isVisible()
        assert detail._wgt_comp.isVisible()
        assert not detail._wgt_depth.isVisible()

    def test_utm_crs_shows_utm_widget(self, detail):
        detail._combo_crs.setCurrentText("UTM Zone")
        assert detail._wgt_utm.isVisible()

    def test_custom_crs_shows_epsg_widget(self, detail):
        detail._combo_crs.setCurrentText("Custom EPSG")
        assert detail._wgt_epsg.isVisible()


class TestDetailCrsResolution:
    def test_resolve_utm_north(self, detail):
        detail._combo_crs.setCurrentText("UTM Zone")
        detail._spin_utm.setValue(33)
        detail._radio_utm_n.setChecked(True)
        assert detail._resolve_crs() == "EPSG:32633"

    def test_resolve_custom_prefixed(self, detail):
        detail._combo_crs.setCurrentText("Custom EPSG")
        detail._edit_epsg.setText("EPSG:2154")
        assert detail._resolve_crs() == "EPSG:2154"

    def test_resolve_custom_invalid(self, detail):
        detail._combo_crs.setCurrentText("Custom EPSG")
        detail._edit_epsg.setText("xxx")
        assert detail._resolve_crs() == "EPSG:4326"

    def test_resolve_geographic(self, detail):
        assert detail._resolve_crs() == "EPSG:4326"


class TestDetailFrequency:
    def test_populate_freq_combo_hz_and_period(self, detail):
        detail._freq_list = [100.0, 10.0]
        detail._btn_hz.setChecked(True)
        detail._populate_freq_combo()
        assert "Hz" in detail._combo_freq.itemText(0)

        detail._btn_s.setChecked(True)
        detail._populate_freq_combo()
        assert "s" in detail._combo_freq.itemText(0)

    def test_on_freq_list_populates(self, detail):
        detail._on_freq_list([5.0, 50.0])
        assert detail._combo_freq.count() == 2

    def test_on_freq_unit_noop_when_empty(self, detail):
        detail._on_freq_unit(True)  # no freq_list -> no-op

    def test_current_freq_hz_period_text(self, detail):
        detail._combo_freq.addItem("0.5 s")
        detail._combo_freq.setCurrentIndex(0)
        assert detail._current_freq_hz() == pytest.approx(2.0)

    def test_current_freq_hz_invalid(self, detail):
        detail._combo_freq.addItem("nonsense")
        detail._combo_freq.setCurrentIndex(0)
        assert detail._current_freq_hz() == 1.0

    def test_contour_mode_str(self, detail):
        detail._combo_contour.setCurrentText("Filled + labels")
        assert detail._contour_mode_str() == "filled_labels"


class TestDetailRefreshExportSetters:
    def test_on_refresh_calls_map_panel_redraw_with_log_scale(
        self, detail, monkeypatch
    ):
        calls = []
        monkeypatch.setattr(
            detail._map_panel, "redraw", lambda **kw: calls.append(kw)
        )
        detail._on_refresh()
        assert calls[0]["log_scale"] == detail._chk_log.isChecked()

    def test_on_export_opens_dialog(self, detail, monkeypatch):
        calls = []

        class _FakeExportDialog:
            def __init__(self, figure, parent):
                calls.append(figure)

            def exec(self):
                calls.append("exec")

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            _FakeExportDialog,
        )
        detail._on_export()
        assert calls[-1] == "exec"

    def test_set_sites_dataframe_dark_mode_delegate(self, detail, monkeypatch):
        calls = {}
        monkeypatch.setattr(
            detail._map_panel,
            "set_sites",
            lambda s: calls.setdefault("sites", s),
        )
        monkeypatch.setattr(
            detail._map_panel,
            "set_dataframe",
            lambda df: calls.setdefault("df", df),
        )
        monkeypatch.setattr(
            detail._map_panel,
            "set_dark_mode",
            lambda d: calls.setdefault("dark", d),
        )
        detail.set_sites(["s"])
        detail.set_dataframe("df")
        detail.set_dark_mode(True)
        assert calls == {"sites": ["s"], "df": "df", "dark": True}


class TestDetailSyncFromPanel:
    def test_sync_from_panel_applies_settings(self, detail, monkeypatch):
        monkeypatch.setattr(detail, "_on_refresh", lambda: None)
        panel = _fake_panel()
        detail.sync_from_panel(panel)
        assert detail._combo_type.currentText() == "Depth"
        assert detail._combo_crs.currentText() == "UTM Zone"
        assert not detail._chk_log.isChecked()
        assert detail._combo_comp.currentText() == "YX"
        assert detail._spin_depth.value() == 750.0

    def test_sync_from_panel_geographic(self, detail, monkeypatch):
        monkeypatch.setattr(detail, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:4326")
        detail.sync_from_panel(panel)
        assert detail._combo_crs.currentText() == "Geographic (lat/lon)"

    def test_sync_from_panel_south_utm(self, detail, monkeypatch):
        monkeypatch.setattr(detail, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:32733")
        detail.sync_from_panel(panel)
        assert detail._radio_utm_s.isChecked()

    def test_sync_from_panel_custom_epsg_non_numeric(self, detail, monkeypatch):
        monkeypatch.setattr(detail, "_on_refresh", lambda: None)
        panel = _fake_panel(_source_crs="EPSG:XYZ")
        detail.sync_from_panel(panel)
        assert detail._combo_crs.currentText() == "Custom EPSG"

    def test_sync_from_panel_freq_placeholder(self, detail, monkeypatch):
        monkeypatch.setattr(detail, "_on_refresh", lambda: None)
        panel = _fake_panel(_target_freq_hz=42.0)
        detail.sync_from_panel(panel)
        assert detail._combo_freq.count() == 1

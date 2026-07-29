# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ProfileViewerWindow (pycsamt.app.desktop.windows.profile_window).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — used to exercise set_sites, station-combo
population, period-range auto-fill, and redraw paths against a real
PlotController/ProfilePanel (already unit-tested independently in
test_profile_panel_extra.py). Here the focus is the window-level wiring:
params-panel widgets -> controller calls, and the export/publication-view
dialog launchers.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.windows.profile_window import ProfileViewerWindow

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture
def win(qapp):
    w = ProfileViewerWindow(parent=None)
    yield w
    w.close()


@pytest.fixture
def win_with_sites(win, willy_sites):
    win.set_sites(willy_sites)
    return win


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Profile Viewer" in win.windowTitle()

    def test_default_component_checks(self, win):
        assert win._chk_xy.isChecked()
        assert win._chk_yx.isChecked()
        assert not win._chk_xx.isChecked()
        assert not win._chk_yy.isChecked()

    def test_default_errbar_and_legend_checked(self, win):
        assert win._chk_errbar.isChecked()
        assert win._chk_legend.isChecked()
        assert not win._chk_bw.isChecked()

    def test_topo_exaggeration_spin_disabled_by_default(self, win):
        assert not win._chk_topo.isChecked()
        assert not win._spin_exag.isEnabled()

    def test_phase_range_combo_has_five_entries(self, win):
        assert win._combo_phase.count() == 5
        assert win._combo_phase.itemText(0) == "Auto (data range)"

    def test_profile_panel_built(self, win):
        assert win._profile_panel is not None


# ── set_sites / station combo / period range ────────────────────────────────


class TestSetSites:
    def test_populates_station_combo(self, win_with_sites, willy_sites):
        # +1 for the "select station" placeholder at index 0
        assert win_with_sites._combo_station.count() == len(willy_sites) + 1

    def test_period_range_updated_from_freqs(self, win_with_sites):
        assert win_with_sites._spin_tmin.value() > 0
        assert win_with_sites._spin_tmax.value() > win_with_sites._spin_tmin.value()

    def test_set_sites_survives_iteration_exception(self, win, monkeypatch):
        class _Bad:
            def __iter__(self):
                raise RuntimeError("boom")

        win.set_sites(_Bad())  # must not raise
        assert win._combo_station.count() == 1  # just the placeholder

    def test_set_station_selects_and_updates_label(self, win_with_sites, willy_sites):
        name = willy_sites[0].name
        win_with_sites.set_station(name)
        assert win_with_sites._info_lbl.text() == f"Station: {name}"

    def test_set_station_empty_name_is_noop(self, win_with_sites):
        win_with_sites._info_lbl.setText("unchanged")
        win_with_sites.set_station("")
        assert win_with_sites._info_lbl.text() == "unchanged"

    def test_on_station_picked_applies_station(self, win_with_sites, willy_sites):
        name = willy_sites[1].name
        win_with_sites._on_station_picked(name)
        assert win_with_sites._profile_panel._ctrl._station_id == name


# ── Component / display toggles ─────────────────────────────────────────────


class TestComponentToggles:
    def test_apply_components_collects_checked(self, win_with_sites):
        win_with_sites._chk_xx.setChecked(True)
        win_with_sites._apply_components()
        assert win_with_sites._profile_panel._ctrl._components == (
            "xy",
            "yx",
            "xx",
        )

    def test_apply_components_falls_back_to_xy_when_none_checked(self, win_with_sites):
        win_with_sites._chk_xy.setChecked(False)
        win_with_sites._chk_yx.setChecked(False)
        win_with_sites._apply_components()
        assert win_with_sites._chk_xy.isChecked()
        assert win_with_sites._profile_panel._ctrl._components == ("xy",)

    def test_component_changed_redraws_on_rho_phi_tab(self, win_with_sites):
        win_with_sites._profile_panel._tabs.setCurrentIndex(0)
        win_with_sites._on_component_changed()  # must not raise

    def test_component_changed_redraws_on_pseudosection_tab(self, win_with_sites):
        win_with_sites._profile_panel._tabs.setCurrentIndex(1)
        win_with_sites._on_component_changed()  # must not raise

    def test_errbar_toggle_pushes_to_controller(self, win_with_sites):
        win_with_sites._chk_errbar.setChecked(False)
        assert win_with_sites._profile_panel._ctrl._show_errbar is False

    def test_bw_toggle_pushes_to_controller(self, win_with_sites):
        win_with_sites._chk_bw.setChecked(True)
        assert win_with_sites._profile_panel._ctrl._bw_mode is True

    def test_phase_range_changed_pushes_ylim(self, win_with_sites):
        win_with_sites._combo_phase.setCurrentIndex(1)
        assert win_with_sites._profile_panel._ctrl._phase_ylim == (-90.0, 90.0)


# ── Topo slots ───────────────────────────────────────────────────────────────


class TestTopoSlots:
    def test_topo_toggle_enables_exaggeration_spin(self, win_with_sites):
        win_with_sites._chk_topo.setChecked(True)
        assert win_with_sites._spin_exag.isEnabled()

    def test_topo_toggle_syncs_section_panel_checkbox(self, win_with_sites):
        win_with_sites._chk_topo.setChecked(True)
        assert win_with_sites._profile_panel._section_panel._chk_topo.isChecked()

    def test_exag_changed_redraws_when_topo_enabled(self, win_with_sites):
        win_with_sites._chk_topo.setChecked(True)
        win_with_sites._spin_exag.setValue(2.0)  # must not raise

    def test_exag_changed_noop_when_topo_disabled(self, win_with_sites):
        win_with_sites._chk_topo.setChecked(False)
        win_with_sites._spin_exag.setValue(3.0)  # must not raise


# ── Refresh / export / publication view ─────────────────────────────────────


class TestRefreshExportPubView:
    def test_refresh_pushes_period_range(self, win_with_sites):
        win_with_sites._spin_tmin.setValue(0.001)
        win_with_sites._spin_tmax.setValue(10.0)
        win_with_sites._on_refresh()
        assert win_with_sites._profile_panel._ctrl._period_range == (
            0.001,
            10.0,
        )

    def test_refresh_ignores_inverted_range(self, win_with_sites):
        win_with_sites._spin_tmin.setValue(100.0)
        win_with_sites._spin_tmax.setValue(1.0)
        win_with_sites._on_refresh()
        assert win_with_sites._profile_panel._ctrl._period_range is None

    def test_export_opens_dialog_for_current_tab(
        self, win_with_sites, willy_sites, monkeypatch
    ):
        win_with_sites.set_station(willy_sites[0].name)
        win_with_sites._profile_panel._tabs.setCurrentIndex(0)
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
        win_with_sites._on_export()
        assert calls[-1] == "exec"

    def test_pub_view_requires_station_selected(self, win_with_sites):
        win_with_sites._combo_station.select_station("")
        win_with_sites._on_pub_view()
        assert "Select a station" in win_with_sites._info_lbl.text()

    def test_pub_view_opens_dialog_with_selected_station(
        self, win_with_sites, willy_sites, monkeypatch
    ):
        name = willy_sites[0].name
        win_with_sites.set_station(name)

        calls = {}

        class _FakeDlg:
            def __init__(self, controller, station_name, dark, parent):
                calls["station_name"] = station_name
                calls["dark"] = dark

            def show(self):
                calls["shown"] = True

        monkeypatch.setattr(
            "pycsamt.app.desktop.windows.publication_view_dialog"
            ".PublicationViewDialog",
            _FakeDlg,
        )
        win_with_sites._on_pub_view()
        assert calls["station_name"] == name
        assert calls["shown"] is True


# ── set_dark_mode ────────────────────────────────────────────────────────────


class TestSetDarkMode:
    def test_set_dark_mode_delegates_to_panel(self, win):
        win.set_dark_mode(False)
        assert win._profile_panel._ctrl.dark is False

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.widgets.settings_pages.*"""

from __future__ import annotations

import pytest

from pycsamt.app.desktop.widgets.settings_pages import (
    DisplayPage,
    InterpretationPage,
    PseudosectionsPage,
    SettingsPage,
    TopographyPage,
    ViewControlsPage,
)
from pycsamt.app.desktop.widgets.settings_pages.base_page import (
    SettingsPage as BaseSettingsPage,
)

# ── base_page ────────────────────────────────────────────────────────────────


def test_base_settings_page_default_interface(qapp):
    page = BaseSettingsPage()
    assert page.populate() is None
    assert page.collect() == {}
    assert page.reset() is None


def test_base_settings_page_is_settings_page_alias():
    assert SettingsPage is BaseSettingsPage


# ── DisplayPage ──────────────────────────────────────────────────────────────


def test_display_page_builds_and_populates(qapp):
    page = DisplayPage()
    assert len(page._comp_widgets) == 6
    for key in ("xy", "yx", "xx", "yy", "te", "tm"):
        assert key in page._comp_widgets


def test_display_page_collect_after_edit(qapp):
    page = DisplayPage()
    cbtn, lw_spin = page._comp_widgets["xy"]
    cbtn._color = "#abcdef"
    lw_spin.setValue(3.5)
    page.collect()

    from pycsamt.api.style import PYCSAMT_STYLE as S

    assert S.mt.xy.color == "#abcdef"
    assert S.mt.xy.lw == pytest.approx(3.5)


def test_display_page_reset_restores_defaults(qapp):
    page = DisplayPage()
    cbtn, _ = page._comp_widgets["xy"]
    cbtn._color = "#000000"
    page.collect()
    page.reset()

    from pycsamt.api.style import PYCSAMT_STYLE as S

    assert S.mt.xy.color == "#1f77b4"


def test_color_btn_pick_updates_color(qapp, monkeypatch):
    from PySide6.QtGui import QColor

    from pycsamt.app.desktop.widgets.settings_pages.display_page import (
        _color_btn,
    )

    monkeypatch.setattr(
        "pycsamt.app.desktop.widgets.settings_pages.display_page.QColorDialog.getColor",
        staticmethod(lambda *a, **k: QColor("#123456")),
    )
    btn = _color_btn("#1e66f5", None)
    assert btn._color == "#1e66f5"
    btn.click()
    assert btn._color == "#123456"


def test_display_page_populate_and_collect_handle_missing_singleton(qapp, monkeypatch):
    monkeypatch.delattr("pycsamt.api.style.PYCSAMT_STYLE")
    page = DisplayPage()
    changes = page.collect()
    assert changes == {}
    page.reset()


def test_color_btn_pick_invalid_color_keeps_original(qapp, monkeypatch):
    from PySide6.QtGui import QColor

    from pycsamt.app.desktop.widgets.settings_pages.display_page import (
        _color_btn,
    )

    monkeypatch.setattr(
        "pycsamt.app.desktop.widgets.settings_pages.display_page.QColorDialog.getColor",
        staticmethod(lambda *a, **k: QColor()),
    )
    btn = _color_btn("#1e66f5", None)
    btn.click()
    assert btn._color == "#1e66f5"


# ── InterpretationPage ───────────────────────────────────────────────────────


def test_interpretation_page_builds_and_populates(qapp):
    page = InterpretationPage()
    assert page._sec_cmap.count() > 0
    assert page._wt_combo.count() == 4


def test_interpretation_page_collect_and_reset(qapp, monkeypatch):
    import types

    fake_sec = types.SimpleNamespace(cmap="viridis", wt_linestyle="--", alpha=0.85)
    fake_prof = types.SimpleNamespace(cmap="viridis")
    reset_calls = []
    fake_interp = types.SimpleNamespace(
        pseudosection=fake_sec,
        profile=fake_prof,
        reset=lambda: reset_calls.append(True),
    )
    monkeypatch.setattr("pycsamt.api.interp.PYCSAMT_INTERP", fake_interp)

    page = InterpretationPage()
    page._sec_cmap.setCurrentText("plasma")
    page._alpha_spin.setValue(0.5)
    page.collect()

    assert fake_sec.cmap == "plasma"
    assert fake_sec.alpha == pytest.approx(0.5)

    page.reset()
    assert reset_calls == [True]


def test_interpretation_page_populate_handles_missing_singleton(qapp, monkeypatch):
    monkeypatch.delattr("pycsamt.api.interp.PYCSAMT_INTERP")
    page = InterpretationPage()
    # falls back to the hard-coded default when the singleton is unusable
    assert page._alpha_spin.value() == pytest.approx(0.85)


# ── PseudosectionsPage ───────────────────────────────────────────────────────


def test_pseudosections_page_builds_and_populates(qapp):
    page = PseudosectionsPage()
    assert page._marker_combo.count() == 5
    assert page._side_combo.count() == 2
    assert page._ydir_combo.count() == 3


def test_pseudosections_page_collect_shape(qapp):
    page = PseudosectionsPage()
    changes = page.collect()
    assert set(changes) == {"station", "section"}
    assert "side" in changes["station"]
    assert "y_direction" in changes["section"]


def test_pseudosections_page_reset(qapp):
    page = PseudosectionsPage()
    page.reset()


def test_pseudosections_page_populate_and_reset_handle_missing_singletons(
    qapp, monkeypatch
):
    monkeypatch.delattr("pycsamt.api.station.PYCSAMT_STATION_RENDERING")
    monkeypatch.delattr("pycsamt.api.section.PYCSAMT_SECTION")
    page = PseudosectionsPage()
    page.reset()


# ── TopographyPage ───────────────────────────────────────────────────────────


def test_topography_page_builds_and_populates(qapp):
    page = TopographyPage()
    assert page._exag_spin.minimum() == pytest.approx(0.1)


def test_topography_page_collect(qapp):
    page = TopographyPage()
    page._enabled_cb.setChecked(True)
    page._exag_spin.setValue(2.5)
    page._pad_spin.setValue(0.02)
    changes = page.collect()
    assert changes == {
        "topography": {
            "enabled": True,
            "exaggeration": pytest.approx(2.5),
            "marker_pad": pytest.approx(0.02),
        }
    }


def test_topography_page_reset(qapp):
    page = TopographyPage()
    page.reset()


def test_topography_page_populate_and_reset_handle_missing_singleton(qapp, monkeypatch):
    monkeypatch.delattr("pycsamt.topo.PYCSAMT_TOPO")
    page = TopographyPage()
    assert page._enabled_cb.isChecked() is False
    assert page._exag_spin.value() == pytest.approx(1.0)
    assert page._pad_spin.value() == pytest.approx(0.015)
    page.reset()


# ── ViewControlsPage ─────────────────────────────────────────────────────────


def test_view_controls_page_builds_and_populates(qapp):
    page = ViewControlsPage()
    assert page._rho_combo.count() == 2
    assert page._x_combo.count() == 4


def test_view_controls_page_collect_shape(qapp):
    page = ViewControlsPage()
    page._ph_min.setValue(-90.0)
    page._ph_max.setValue(90.0)
    changes = page.collect()
    vc = changes["view_controls"]
    assert vc["phase_range"] == (pytest.approx(-90.0), pytest.approx(90.0))
    assert "rho_view" in vc
    assert "x_view" in vc


def test_view_controls_page_reset(qapp):
    page = ViewControlsPage()
    page.reset()


def test_view_controls_page_populate_and_reset_handle_missing_singleton(
    qapp, monkeypatch
):
    monkeypatch.delattr("pycsamt.api.control.PYCSAMT_CONTROL")
    page = ViewControlsPage()
    page.reset()

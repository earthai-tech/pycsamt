# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Extra tests for ProfilePanel, on top of test_profile_panel.py (which
covers construction and no-data draw calls). This file drives the
real-data paths: set_sites with real EDI Sites, per-tab
set_selected_station branches, the phase-tensor redraw cache fast-path,
invalidate_phase_tensor, set_inversion_result, and the exception-swallow
branches in set_sites.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.profile_panel import ProfilePanel

_ROOT = Path(__file__).parents[3]
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))


@pytest.fixture(scope="session")
def tipper_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


@pytest.fixture
def panel(qapp):
    p = ProfilePanel()
    yield p
    p.close()


def _first_station(sites):
    return sites.as_list()[0].station


# ── set_sites with real data ─────────────────────────────────────────────


def test_set_sites_populates_freq_range(panel, tipper_sites):
    panel.set_sites(tipper_sites)
    assert panel._freq_sel is not None  # range_changed emission side-effect
    assert panel._ctrl._sites is tipper_sites


def test_set_sites_resets_pt_cache_key(panel, tipper_sites):
    panel._pt_last_key = ("stale",)
    panel.set_sites(tipper_sites)
    # set_sites clears it up front, then _redraw_all repopulates it via
    # _redraw_phase_tensor — either way it must not still be the stale value.
    assert panel._pt_last_key != ("stale",)


def test_set_sites_iteration_exception_is_swallowed(panel):
    class _BadIter:
        def __iter__(self):
            raise RuntimeError("boom")

    panel.set_sites(_BadIter())  # must not raise


def test_set_sites_redraw_exception_is_swallowed(panel, tipper_sites, monkeypatch):
    monkeypatch.setattr(
        panel, "_redraw_all", mock.Mock(side_effect=RuntimeError("boom"))
    )
    panel.set_sites(tipper_sites)  # must not raise


# ── set_selected_station per-tab branches ────────────────────────────────


@pytest.mark.parametrize("tab_index", [0, 1, 2, 3, 4, 5])
def test_set_selected_station_per_tab(panel, tipper_sites, tab_index):
    panel.set_sites(tipper_sites)
    panel._tabs.setCurrentIndex(tab_index)
    station = _first_station(tipper_sites)
    panel.set_selected_station(station)  # must not raise
    assert panel._ctrl._station_id == station


# ── Phase-tensor redraw cache fast path ──────────────────────────────────


def test_redraw_phase_tensor_cache_hit_skips_recompute(panel, tipper_sites):
    panel.set_sites(tipper_sites)
    panel._tabs.setCurrentIndex(4)
    panel._redraw_phase_tensor()
    with mock.patch.object(panel._ctrl, "draw_phase_tensor") as m:
        panel._redraw_phase_tensor()  # same key, axes present -> fast path
        m.assert_not_called()


def test_redraw_phase_tensor_cache_miss_after_invalidate(panel, tipper_sites):
    panel.set_sites(tipper_sites)
    panel._redraw_phase_tensor()
    panel.invalidate_phase_tensor()
    with mock.patch.object(
        panel._ctrl, "draw_phase_tensor", wraps=panel._ctrl.draw_phase_tensor
    ) as m:
        panel._redraw_phase_tensor()
        m.assert_called_once()


def test_invalidate_phase_tensor_resets_key_and_delegates(panel):
    panel._pt_last_key = ("x",)
    with mock.patch.object(panel._ctrl, "invalidate_phase_tensor") as m:
        panel.invalidate_phase_tensor()
        m.assert_called_once()
    assert panel._pt_last_key is None


# ── set_dark_mode delegates to section panel ─────────────────────────────


def test_set_dark_mode_delegates_to_section_panel(panel):
    with mock.patch.object(panel._section_panel, "set_dark_mode") as m:
        panel.set_dark_mode(False)
        m.assert_called_once_with(False)
    assert panel._pt_last_key is None


# ── set_inversion_result ─────────────────────────────────────────────────


def test_set_inversion_result_delegates_and_switches_tab(panel):
    fake_result = object()
    with mock.patch.object(panel._section_panel, "set_result") as m:
        panel.set_inversion_result(fake_result)
        m.assert_called_once_with(fake_result)
    assert panel._tabs.currentWidget() is panel._section_panel


# ── _on_freq_range_changed edge cases ─────────────────────────────────────


def test_on_freq_range_changed_zero_lo_hz(panel):
    # lo_hz<=0 -> T_max is None going into set_period_range, which replaces
    # a None T_max with its "no upper limit" sentinel 1e9 (`T_max or 1e9`).
    # T_min still derives from hi_hz.
    panel._on_freq_range_changed(0.0, 100.0)
    assert panel._ctrl._period_range[1] == pytest.approx(1e9)
    assert panel._ctrl._period_range[0] == pytest.approx(1.0 / 100.0)


def test_on_freq_range_changed_zero_hi_hz(panel):
    # hi_hz<=0 -> T_min is None going into set_period_range, which replaces
    # a None T_min with its "no lower limit" sentinel 0.0 (`T_min or 0.0`).
    # T_max still derives from lo_hz.
    panel._on_freq_range_changed(1.0, 0.0)
    assert panel._ctrl._period_range[0] == pytest.approx(0.0)
    assert panel._ctrl._period_range[1] == pytest.approx(1.0 / 1.0)

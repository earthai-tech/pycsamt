# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for Edit ▸ Frequencies ▸ Point Editor: the point edits on real
kap03 EDIs and Gabbs Valley EMTF-XML (USGS, doi:10.5066/P9GZ9Z56), the
dialog, and the main-window apply / undo round trip."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers import point_edits as pe

ROOT = Path(__file__).resolve().parents[4]
KAP03 = ROOT / "data" / "MT" / "kap03lmt_edis"
GV_XML = ROOT / "data" / "gv_data" / "xml"


def _load(path):
    if not path.is_dir():
        pytest.skip(f"{path.name} missing")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(path))


@pytest.fixture(scope="module")
def kap03():
    return _load(KAP03)


@pytest.fixture(scope="module", params=["kap03", "gv_xml"])
def site(request):
    sites = _load(KAP03 if request.param == "kap03" else GV_XML)
    return list(sites)[1]


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.CRITICAL)  # NaN-row recompute noise
    yield
    logging.disable(logging.NOTSET)


# ── point edits (Qt-free) ───────────────────────────────────────────────────


class TestPointEdits:
    def test_mask_one_component(self, site):
        c0 = pe.curves(site)
        m = pe.mask(site, [3, 4], ("xy",))
        c = pe.curves(m)
        assert not c["xy"]["valid"][[3, 4]].any()
        assert (c["yx"]["valid"] == c0["yx"]["valid"]).all()
        assert "masked XY" in pe.info_log(m)
        assert c0["xy"]["valid"][[3, 4]].all()  # input untouched

    def test_interpolate_lies_between_neighbours(self, site):
        c0 = pe.curves(site)
        it = pe.interpolate(pe.mask(site, [5], ("xy",)), [5], ("xy",))
        rho = pe.curves(it)["xy"]["rho"]
        lo, hi = sorted((c0["xy"]["rho"][4], c0["xy"]["rho"][6]))
        assert lo * 0.999 <= rho[5] <= hi * 1.001
        assert "interpolated XY" in pe.info_log(it)

    def test_restore_brings_back_the_original(self, site):
        m = pe.mask(site, [2, 3], ("xy", "yx"))
        r = pe.restore(m, site, [2, 3])
        assert np.allclose(np.asarray(r.z), np.asarray(site.z),
                           equal_nan=True)

    def test_delete_rows(self, site):
        n = pe.curves(site)["freq"].size
        d = pe.delete(site, [0, 1])
        assert pe.curves(d)["freq"].size == n - 2
        assert "deleted" in pe.info_log(d)
        with pytest.raises(ValueError):
            pe.delete(site, range(n))

    def test_static_shift_keeps_phase(self, site):
        c0 = pe.curves(site)
        sh = pe.static_shift(site, "yx", 2.0)
        c = pe.curves(sh)
        ok = np.isfinite(c0["yx"]["rho"])
        assert np.allclose(c["yx"]["rho"][ok] / c0["yx"]["rho"][ok], 2.0)
        assert np.allclose(c["yx"]["phi"][ok], c0["yx"]["phi"][ok])
        assert np.allclose(c["xy"]["rho"][ok], c0["xy"]["rho"][ok],
                           equal_nan=True)
        with pytest.raises(ValueError):
            pe.static_shift(site, "xy", 0)

    def test_nearest_rows(self):
        assert pe.nearest_rows([100.0, 10.0, 1.0], [10.05, 5.0]) == [1]


# ── dialog ──────────────────────────────────────────────────────────────────


@pytest.fixture
def dlg(qapp, kap03):
    from pycsamt.app.desktop.dialogs.point_editor import PointEditorDialog

    d = PointEditorDialog(kap03)
    d.show()
    yield d
    d.close()


class TestDialog:
    def test_opens_on_the_given_station(self, qapp, kap03):
        from pycsamt.app.desktop.dialogs.point_editor import (
            PointEditorDialog,
        )

        name = list(kap03)[4].name
        d = PointEditorDialog(kap03, station=name)
        assert d.station == name
        assert not d._btn_apply.isEnabled()
        d.close()

    def test_select_mask_undo(self, dlg):
        dlg.select([2, 3], ["xy"])
        assert "2 point(s) selected" in dlg._sel_lbl.text()
        dlg.edit("mask")
        cur = dlg._cur[dlg.station]
        assert not pe.curves(cur)["xy"]["valid"][[2, 3]].any()
        assert dlg.edited_stations() == [dlg.station]
        assert "(1 station)" in dlg._btn_apply.text()
        dlg.undo()
        assert dlg.edited_stations() == []

    def test_all_stations(self, dlg):
        f0 = pe.curves(dlg._cur[dlg.station])["freq"][0]
        having = [n for n in dlg._order
                  if pe.nearest_rows(pe.curves(dlg._cur[n])["freq"], [f0])]
        dlg.select([0], ["xy", "yx"])
        dlg._all_stations.setChecked(True)
        dlg.edit("mask")
        # every station recording that frequency, and only those
        assert dlg.edited_stations() == having
        missing = len(dlg._order) - len(having)
        if missing:
            assert f"{missing} station(s) do not have" in dlg._log_lbl.text()

    def test_edit_without_selection_explains(self, dlg):
        dlg.edit("delete")
        assert "Select points first" in dlg._log_lbl.text()
        assert not dlg.edited_stations()

    def test_interpolated_points_are_marked(self, dlg):
        dlg.select([4], ["xy"])
        dlg.edit("interp")
        assert (4, "xy") in dlg._interp[dlg.station]
        dlg._draw()  # hollow diamonds drawn without error

    def test_shift_and_navigation(self, dlg):
        first = dlg.station
        dlg.shift("xy", 0.5)
        assert "Static shift XY" in dlg._log_lbl.text()
        dlg.step_station(+1)
        assert dlg.station != first
        dlg.step_station(-1)
        assert dlg.station == first
        assert dlg.edited_stations() == [first]


# ── main window ─────────────────────────────────────────────────────────────


@pytest.fixture
def window(qapp, monkeypatch, kap03):
    from pycsamt.app.desktop.main_window import MainWindow
    from pycsamt.app.desktop.models.session import SessionState

    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    win = MainWindow()
    assert win._adopt_full_dataset(kap03, record=False, merge=False,
                                   lines={s.name: "Survey" for s in kap03})
    win._history.reset(win._survey_state())
    yield win
    win._history.mark_saved()
    win.close()


def test_apply_point_edits_is_one_undo_step(window, kap03):
    from pycsamt.app.desktop.dialogs.point_editor import PointEditorDialog

    d = PointEditorDialog(window._all_sites)
    d.select([1], ["xy"])
    d.edit("mask")
    name = d.station
    assert window.apply_point_edits(d.edited_sites(), d.edited_stations())
    assert window._act_undo.text() == f"&Undo Point Editor: {name}"
    site = next(s for s in window._all_sites if s.name == name)
    assert not pe.curves(site)["xy"]["valid"][1]
    window.undo_edit()
    site = next(s for s in window._all_sites if s.name == name)
    assert pe.curves(site)["xy"]["valid"][1]
    df = window._all_dataframe
    assert (df["Line"] == "Survey").all()
    d.close()

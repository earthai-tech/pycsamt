# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for Edit ▸ Frequencies / Tensor whole-survey operations: the two
library fixes they rely on (no values invented outside a station's band
when regridding; decimate keeping z / z_err / freq aligned), every
operation on the real kap03 EDIs, the dialog, and the main-window undo."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from pycsamt.app.desktop.controllers import point_edits as pe
from pycsamt.app.desktop.controllers import survey_ops as so

KAP03 = Path(__file__).resolve().parents[4] / "data" / "MT" / "kap03lmt_edis"


@pytest.fixture(scope="module")
def kap03():
    if not KAP03.is_dir():
        pytest.skip("kap03 EDIs missing")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(KAP03))


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


def _nf(site) -> int:
    return int(np.isfinite(pe.curves(site)["freq"]).sum())


# ── library fixes ───────────────────────────────────────────────────────────


def test_regrid_leaves_out_of_band_points_empty():
    from pycsamt.emtools.frequency import _regrid_t, _regrid_z

    fr = np.array([100.0, 10.0, 1.0])
    z = np.ones((3, 2, 2), dtype=complex) * (1 + 1j)
    fnew = np.array([1000.0, 100.0, 3.0, 1.0, 0.01])
    for method in ("nearest", "linear"):
        out = _regrid_z(z, fr, fnew, method=method)
        assert np.isnan(out[[0, 4]]).all()  # outside 1-100 Hz
        assert np.isfinite(out[1:4]).all()
        t = _regrid_t(np.ones((3, 2), dtype=complex), fr, fnew,
                      method=method)
        assert np.isnan(t[[0, 4]]).all() and np.isfinite(t[1:4]).all()


def test_decimate_keeps_arrays_aligned(kap03):
    from pycsamt.emtools.frequency import decimate_step

    out = decimate_step(kap03, step=2)
    for s0, s1 in zip(kap03, out):
        n0 = _nf(s0)
        assert _nf(s1) == (n0 + 1) // 2  # was 0: the freq write failed
        assert np.asarray(s1.z).shape[0] == np.asarray(s1.z_err).shape[0]
        assert np.allclose(pe.curves(s1)["freq"], pe.curves(s0)["freq"][::2])


# ── operations ──────────────────────────────────────────────────────────────


class TestOperations:
    def test_trim(self, kap03):
        out = so.run("trim", kap03, {"fmin": 0.001, "fmax": 0.05})
        for s in out:
            f = pe.curves(s)["freq"]
            assert ((f >= 0.001 * 0.999) & (f <= 0.05 * 1.001)).all()
        with pytest.raises(ValueError):
            so.run("trim", kap03, {"fmin": 0.0, "fmax": 0.0})
        with pytest.raises(ValueError):
            so.run("trim", kap03, {"fmin": 1.0, "fmax": 0.1})

    def test_regrid_is_log_spaced(self, kap03):
        out = so.run("regrid", kap03, {"per_decade": 5, "fmin": 0.001,
                                       "fmax": 1.0, "method": "nearest"})
        f = np.sort(pe.curves(next(iter(out)))["freq"])
        steps = np.diff(np.log10(f))
        assert np.allclose(steps, steps[0])

    def test_align_union_invents_nothing(self, kap03):
        out = so.run("align", kap03, {"mode": "union", "ref": "",
                                      "method": "nearest"})
        grids = {tuple(np.round(pe.curves(s)["freq"], 12)) for s in out}
        assert len(grids) == 1  # one common grid
        for s0, s1 in zip(kap03, out):
            f0 = pe.curves(s0)["freq"]
            c1 = pe.curves(s1)
            outside = (c1["freq"] < f0.min() * 0.999) | (
                c1["freq"] > f0.max() * 1.001)
            assert not c1["xy"]["valid"][outside].any()

    def test_align_reference_needs_a_station(self, kap03):
        with pytest.raises(ValueError, match="reference"):
            so.run("align", kap03, {"mode": "ref", "ref": "",
                                    "method": "nearest"})

    def test_fill_gaps_only_between_measured(self, kap03):
        from pycsamt.site.base import Sites

        site = list(kap03)[1]
        n = _nf(site)
        m = pe.mask(site, [0, 5, 6], ("xy",))  # an edge and an interior
        out = list(so.fill_gaps(Sites([m]), ("xy",)))[0]
        valid = pe.curves(out)["xy"]["valid"]
        assert valid[[5, 6]].all()
        assert not valid[0]  # band edge: never extrapolated
        assert valid.sum() == n - 1

    def test_rotate(self, kap03):
        out = so.run("rotate", kap03, {"angle": 30.0})
        z0 = np.asarray(next(iter(kap03)).z)
        z1 = np.asarray(next(iter(out)).z)
        assert not np.allclose(z0, z1, equal_nan=True)
        with pytest.raises(ValueError):
            so.run("rotate", kap03, {"angle": 0.0})

    def test_strike_preview_reports_angles(self, kap03):
        out = so.run("strike", kap03, {"method": "swift"})
        t = so.preview_table(kap03, out, with_strike=True)
        assert t["strike"].notna().sum() > len(t) // 2

    def test_input_untouched(self, kap03):
        z0 = np.asarray(next(iter(kap03)).z).copy()
        so.run("decimate", kap03, {"step": 3})
        so.run("rotate", kap03, {"angle": 45.0})
        assert np.allclose(np.asarray(next(iter(kap03)).z), z0,
                           equal_nan=True)


# ── dialog ──────────────────────────────────────────────────────────────────


@pytest.fixture
def dlg(qapp, kap03):
    pytest.importorskip("PySide6")
    from pycsamt.app.desktop.dialogs.survey_ops import SurveyOpsDialog

    d = SurveyOpsDialog(kap03, op="decimate")
    yield d
    d.close()


class TestDialog:
    def test_preview_and_apply(self, dlg):
        assert dlg.op_key == "decimate"
        assert dlg._table.rowCount() == len(list(dlg._sites))
        assert "frequencies 518 →" in dlg._status.text()
        assert dlg._btn_apply.isEnabled()
        dlg._on_apply()
        assert dlg.result is not None
        assert _nf(next(iter(dlg.result))) < _nf(next(iter(dlg._sites)))

    def test_nothing_to_change_disables_apply(self, dlg):
        dlg.select_op("dedupe")
        assert "nothing would change" in dlg._status.text()
        assert not dlg._btn_apply.isEnabled()

    def test_invalid_settings_are_explained(self, dlg):
        dlg.select_op("trim")
        assert "⚠" in dlg._status.text()
        assert not dlg._btn_apply.isEnabled()
        dlg.set_values({"fmin": 0.001, "fmax": 0.05})
        assert dlg._btn_apply.isEnabled()

    def test_reference_station_choices(self, dlg):
        dlg.select_op("align")
        combo = dlg._forms["align"].widgets["ref"]
        names = [combo.itemData(i) for i in range(combo.count())]
        assert names[0] == "" and set(names[1:]) == {
            s.name for s in dlg._sites}


# ── main window ─────────────────────────────────────────────────────────────


def test_apply_survey_op_is_one_undo_step(qapp, monkeypatch, kap03):
    pytest.importorskip("PySide6")
    from pycsamt.app.desktop.main_window import MainWindow
    from pycsamt.app.desktop.models.session import SessionState

    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    win = MainWindow()
    try:
        assert win._adopt_full_dataset(kap03, record=False, merge=False,
                                       lines={s.name: "Survey" for s in kap03})
        win._history.reset(win._survey_state())
        n0 = int(win._all_dataframe["N_freq"].sum())
        out = so.run("decimate", win._all_sites, {"step": 2})
        assert win.apply_survey_op(out, "Decimate")
        assert int(win._all_dataframe["N_freq"].sum()) < n0
        assert win._act_undo.text() == "&Undo Decimate"
        win.undo_edit()
        assert int(win._all_dataframe["N_freq"].sum()) == n0
        assert (win._all_dataframe["Line"] == "Survey").all()
        texts = [a.text() for a in win._edit_tensor_menu.actions()]
        assert "Rotate to Strike…" in texts
    finally:
        win._history.mark_saved()
        win.close()

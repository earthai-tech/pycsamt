# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the TDEM studio (controllers/tdem_studio + windows/tdem_window).

Real data: two profiles (TEM100, TEM1020) of the bundled JIANGSU TEMAVG
survey, copied to a temporary folder so a load takes a fraction of the
whole 55-profile survey.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers import tdem_studio as st
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.controllers.inversion_engines import defaults

_ROOT = Path(__file__).resolve().parents[4]
_JIANGSU = _ROOT / "data" / "TEMAVG" / "JIANGSU"
_STEMS = ("TEM100", "TEM1020")
pytestmark = pytest.mark.skipif(
    not (_JIANGSU / "TEM100.AVG").exists(), reason="JIANGSU data missing")
# The station coordinate table is not redistributed (gitignored); without it
# the coordinate-based views legitimately refuse to draw.
_HAS_COORDS = bool(list(_JIANGSU.glob("Coordinate*")))
_COORD_VIEWS = {"map", "elevation", "overview"}


def _views():
    return [v for v in st.VIEWS
            if _HAS_COORDS or v.key not in _COORD_VIEWS]


@pytest.fixture(scope="module")
def folder(tmp_path_factory):
    out = tmp_path_factory.mktemp("tem")
    for stem in _STEMS:
        for ext in (".AVG", ".Z", ".LOG"):
            src = _JIANGSU / f"{stem}{ext}"
            if src.exists():
                shutil.copy(src, out / src.name)
    for xls in _JIANGSU.glob("Coordinate*"):  # the station map needs it
        shutil.copy(xls, out / xls.name)
    return out


@pytest.fixture(scope="module")
def data(folder):
    return st.load("temavg", str(folder), defaults(st.source("temavg")
                                                   .fields))


@pytest.fixture
def win(qapp):
    from pycsamt.app.desktop.windows.tdem_window import TDEMWindow

    w = TDEMWindow(parent=None)
    w.show()
    yield w
    w.close()
    plt.close("all")


# ── studio (Qt-free) ──────────────────────────────────────────────────────
class TestStudio:
    def test_load_groups_profiles(self, data):
        assert set(data.profiles) == set(_STEMS)
        assert data.n == sum(len(v) for v in data.profiles.values())
        assert data.survey is not None

    def test_every_view_and_choice_renders(self, data):
        sel = data.profiles["TEM100"][:3]
        for v in _views():
            base = defaults(list(v.fields))
            variants = [base] + [dict(base, **{f.key: c})
                                 for f in v.fields if f.kind == "choice"
                                 for c, _l in f.choices[1:]]
            for vals in variants:
                fig = st.render(v.key, data, selected=sel,
                                profile="TEM100", values=vals)
                assert figure_blank_reason(fig) is None, (v.key, vals)
                plt.close("all")

    def test_whole_survey_sections(self, data):
        for key in ("avg_section", "z_section", "gates"):
            fig = st.render(key, data, profile="")
            assert figure_blank_reason(fig) is None
            plt.close("all")

    def test_curve_views_need_a_selection(self, data):
        with pytest.raises(ValueError, match="select"):
            st.render("decay", data, selected=())

    def test_dashboard_needs_a_profile(self, data):
        with pytest.raises(ValueError, match="profile"):
            st.render("dashboard", data, selected=[0], profile="")

    def test_survey_views_need_temavg(self, data):
        single = st.TDEMData(source="xyz", soundings=data.soundings[:2],
                             profiles={"survey": [0, 1]})
        with pytest.raises(ValueError, match="TEMAVG"):
            st.render("map", single)
        fig = st.render("decay", single, selected=[0, 1])
        assert figure_blank_reason(fig) is None
        plt.close("all")

    @pytest.mark.parametrize("wf", ["none", "square", "ramp", "halfsine"])
    @pytest.mark.parametrize("method", ["late_time", "fourier"])
    def test_convert(self, data, wf, method):
        snds = [data.soundings[i] for i in data.profiles["TEM100"][:3]]
        sites = st.convert(snds, dict(defaults(list(st.CONVERT_FIELDS)),
                                      waveform=wf, method=method))
        assert len(sites) == 3
        assert all(len(s.freq) > 3 for s in sites)
        # a waveform is attached to copies, never to the loaded soundings
        assert all(s.waveform is None for s in snds) or wf == "none"

    def test_convert_nothing_raises(self):
        with pytest.raises(ValueError):
            st.convert([])

    def test_selection_every_keeps_ends(self):
        assert st.selection_every(range(10), 4) == [0, 4, 8, 9]
        assert st.selection_every([3, 4], 5) == [3, 4]
        assert st.selection_every(range(5), 1) == [0, 1, 2, 3, 4]

    def test_unknown_source(self, folder):
        with pytest.raises(ValueError, match="unknown source"):
            st.load("nope", str(folder))

    def test_xyz_source(self, tmp_path):
        import numpy as np

        t = np.logspace(-5, -2, 20)
        p = tmp_path / "s.txt"
        np.savetxt(p, np.c_[t, 1e-6 * t ** -2.5])
        d = st.load("xyz", str(p), defaults(list(st.source("xyz").fields)))
        assert d.n == 1 and d.survey is None


# ── window ────────────────────────────────────────────────────────────────
class TestWindow:
    def test_empty_state(self, win):
        assert not win.has_data
        assert not win._btn_convert.isEnabled()
        assert not win._btn_send.isEnabled()
        assert "No TDEM" in win._data_status.text()

    def test_source_switches_acquisition_form(self, win):
        win.set_source("walktem")
        assert win._acq_stack.currentWidget() is win._acq_forms["walktem"]
        assert win._btn_load.text() == "Load file…"
        win.set_source("temavg")
        assert win._btn_load.text() == "Load folder…"

    def test_load_populates_profiles_and_draws(self, win, folder):
        n = win.load(str(folder), source="temavg")
        assert n > 0 and win.has_data
        items = [win._profile_combo.itemData(i)
                 for i in range(win._profile_combo.count())]
        assert items[0] == "" and set(items[1:]) == set(_STEMS)
        assert win.profile() in _STEMS  # first profile, not the survey
        assert 0 < len(win.selected()) <= 6
        assert win._snd_table.rowCount() == n
        assert win._canvas_view.canvas.figure.axes  # decay drawn

    def test_profile_change_lists_its_soundings(self, win, folder):
        win.load(str(folder), source="temavg")
        win.select_profile("TEM1020")
        assert win._snd_list.count() == len(win._data.profiles["TEM1020"])
        win.select_profile("")
        assert win._snd_list.count() == win._data.n

    def test_check_every_all_none(self, win, folder):
        win.load(str(folder), source="temavg")
        idx = win._data.profiles[win.profile()]
        win.check_every(3)
        assert win.selected() == st.selection_every(idx, 3)
        win._check_all(False)
        assert win.selected() == []
        assert not win._btn_convert.isEnabled()
        win._check_all(True)
        assert len(win.selected()) == len(idx)

    def test_every_view_draws(self, win, folder):
        win.load(str(folder), source="temavg")
        for v in _views():
            win.select_view(v.key)
            win._on_draw()
            assert win._status.text().startswith("✓"), (v.key,
                                                         win._status.text())

    def test_view_error_shows_reason_card(self, win, folder):
        win.load(str(folder), source="temavg")
        win._check_all(False)
        win.select_view("decay")
        win._on_draw()
        assert win._status.text().startswith("✕")

    def test_convert_save_send(self, win, folder, tmp_path):
        win.load(str(folder), source="temavg")
        win.check_every(10)
        sites = win.convert()
        assert len(sites) == len(win.selected())
        assert win._conv_table.rowCount() == len(sites)
        assert win._btn_save.isEnabled() and win._btn_send.isEnabled()
        paths = win.save_edis(str(tmp_path / "edi"))
        assert len(paths) == len(sites)
        got = []
        win.send_to_survey.connect(lambda s, m: got.append((len(s), m)))
        win.send("replace")
        assert got == [(len(sites), "replace")]

    def test_clear(self, win, folder):
        win.load(str(folder), source="temavg")
        win._on_clear()
        assert not win.has_data and win._profile_combo.count() == 0

    def test_load_dialog_cancel_is_noop(self, win, monkeypatch):
        from PySide6.QtWidgets import QFileDialog

        monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                            staticmethod(lambda *a, **k: ""))
        win._on_load()
        assert not win.has_data

    def test_background_load(self, win, folder, monkeypatch, qapp):
        from PySide6.QtWidgets import QFileDialog

        monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                            staticmethod(lambda *a, **k: str(folder)))
        win._on_load()
        win._worker.wait(60000)
        qapp.processEvents()
        assert win.has_data

    def test_background_load_failure(self, win, tmp_path, monkeypatch,
                                     qapp):
        from PySide6.QtWidgets import QFileDialog

        monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                            staticmethod(lambda *a, **k: str(tmp_path)))
        win._on_load()
        win._worker.wait(60000)
        qapp.processEvents()
        assert "Load failed" in win._data_status.text()


class TestMainWindowSend:
    def test_append_and_replace_are_undoable(self, qapp, data):
        from pycsamt.app.desktop.main_window import MainWindow

        mw = MainWindow()
        try:
            snds = [data.soundings[i] for i in data.profiles["TEM100"][:3]]
            sites = st.convert(snds)
            mw._on_tdem_sites(sites, "replace")
            assert mw._controller.n_stations == 3
            more = st.convert([data.soundings[i]
                               for i in data.profiles["TEM1020"][:2]])
            mw._on_tdem_sites(more, "append")
            assert mw._controller.n_stations == 5
            mw.undo_edit()
            assert mw._controller.n_stations == 3
        finally:
            mw.close()

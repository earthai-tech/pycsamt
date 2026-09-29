# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for Edit ▸ Stations: the Qt-free helpers (station_edits), the
Station Editor dialog and the main-window apply / undo round trip, on the
real kap03 EDIs."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers import station_edits as se

KAP03 = Path(__file__).resolve().parents[4] / "data" / "MT" / "kap03lmt_edis"


@pytest.fixture(scope="module")
def kap03():
    if not KAP03.is_dir():
        pytest.skip("kap03 EDIs missing")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(KAP03))


# ── helpers (Qt-free) ───────────────────────────────────────────────────────


class TestHelpers:
    def test_rename_preview_rules(self):
        names = ["kap103", "kap104", "L1-3"]
        assert se.rename_preview(names, prefix="K-") == [
            "K-kap103", "K-kap104", "K-L1-3"]
        assert se.rename_preview(names, find="kap", replace="S") == [
            "S103", "S104", "L1-3"]
        assert se.rename_preview(["L1-3"], pad=3) == ["L1-003"]
        assert se.rename_preview(["ab12"], find=r"[a-z]+", replace="X",
                                 regex=True, case="lower") == ["x12"]

    def test_check_names(self):
        assert se.check_names(["a", "b"]) == []
        assert "duplicate names: a" in se.check_names(["a", "A"])
        assert "empty name" in se.check_names(["a", " "])

    def test_table_changes_and_carry_lines(self):
        before = pd.DataFrame({"station": ["a", "b"], "line": ["L1", "L1"],
                               "lat": [1.0, 2.0]})
        after = before.copy()
        after.loc[0, "station"] = "a2"
        after.loc[1, "lat"] = 2.5
        ch = se.table_changes(before, after)
        assert ch == {"a": {"name": "a2"}, "b": {"lat": 2.5}}
        assert se.carry_lines({"a": "L1", "b": "L1"}, ch) == {"a2": "L1",
                                                               "b": "L1"}
        assert se.carry_lines({"a": "L1"}, {"a": {"line": "L9"}}) == {
            "a": "L9"}

    def test_parse_pasted(self):
        assert se.parse_pasted("1\t2\n3\t4\n") == [["1", "2"], ["3", "4"]]
        assert se.parse_pasted("1,2") == [["1", "2"]]
        assert se.parse_pasted("") == []

    def test_projection_round_trip(self):
        pytest.importorskip("pyproj")
        t = pd.DataFrame({"lat": [-32.13], "lon": [20.46]})
        epsg = se.utm_epsg(-32.13, 20.46)
        assert epsg == 32734
        back = se.unproject_table(se.project_table(t, epsg), epsg)
        assert back["lat"][0] == pytest.approx(-32.13, abs=1e-7)
        assert back["lon"][0] == pytest.approx(20.46, abs=1e-7)

    def test_detect_lines(self):
        assert se.detect_lines(["18-001", "18-002", "19-001"]) == {
            "18-001": "L18", "18-002": "L18", "19-001": "L19"}

    def test_apply_changes_on_real_sites(self, kap03):
        name = next(iter(kap03)).name
        new, n, lines = se.apply_changes(kap03, {
            name: {"lat": -32.0, "elev": 111.0, "name": name + "b",
                   "acqby": "Field crew", "line": "L7"}})
        s = next(iter(new))
        assert (s.name, s.coords[0], s.coords[2]) == (name + "b", -32.0,
                                                      111.0)
        assert se.station_table([s])["acqby"][0] == "Field crew"
        assert n == 1 and lines == {name + "b": "L7"}
        assert next(iter(kap03)).name == name  # input untouched

    def test_apply_changes_rejects_bad_numbers(self, kap03):
        name = next(iter(kap03)).name
        with pytest.raises(ValueError, match="lat must be a number"):
            se.apply_changes(kap03, {name: {"lat": "north"}})
        with pytest.raises(Exception):  # outside [-90, 90]: validated
            se.apply_changes(kap03, {name: {"lat": 120.0}})


# ── dialog ──────────────────────────────────────────────────────────────────


@pytest.fixture
def dlg(qapp, kap03):
    from pycsamt.app.desktop.dialogs.station_editor import (
        StationEditorDialog,
    )

    d = StationEditorDialog(kap03, {s.name: "Survey" for s in kap03})
    yield d
    d.close()


class TestDialog:
    def test_no_changes_disables_apply(self, dlg):
        assert not dlg._btn_apply.isEnabled()
        assert "No changes" in dlg._summary.text()

    def test_coordinate_edit_is_highlighted_and_counted(self, dlg):
        it = dlg._coords.item(0, 2)
        it.setText("-31.5")
        assert it.background().color().alpha() > 0
        assert dlg.compute_changes() == {dlg._base["station"][0]:
                                         {"lat": -31.5}}
        assert "1 change(s) in 1 station(s)" in dlg._summary.text()
        assert dlg._btn_apply.isEnabled()

    def test_paste_block(self, dlg, qapp):
        from PySide6.QtWidgets import QApplication

        QApplication.clipboard().setText("-31.1\t20.1\n-31.2\t20.2")
        dlg._coords.setCurrentCell(0, 2)
        dlg._coords.paste_block()
        ch = dlg.compute_changes()
        names = list(dlg._base["station"])
        assert ch[names[0]] == {"lat": -31.1, "lon": 20.1}
        assert ch[names[1]] == {"lat": -31.2, "lon": 20.2}

    def test_rename_rules_and_duplicate_block(self, dlg):
        dlg._prefix.setText("SA-")
        dlg._preview_names()
        ch = dlg.compute_changes()
        assert all(d["name"].startswith("SA-") for d in ch.values())
        dlg._names.item(1, 1).setText(dlg._names.item(0, 1).text())
        assert not dlg._btn_apply.isEnabled()
        assert "duplicate" in dlg._names_problems.text()

    def test_projected_coordinates(self, dlg):
        pytest.importorskip("pyproj")
        dlg._crs.setCurrentIndex(1)  # suggested UTM
        assert "easting" in dlg._coords.keys
        c = dlg._coords.keys.index("easting")
        e0 = float(dlg._coords.item(0, c).text())
        dlg._coords.item(0, c).setText(f"{e0 + 1000:.3f}")
        d = dlg.compute_changes()[dlg._base["station"][0]]
        assert "lon" in d and d["lon"] > dlg._base["lon"][0]

    def test_import_csv(self, dlg, tmp_path):
        names = list(dlg._base["station"])
        p = tmp_path / "c.csv"
        p.write_text("Station,Latitude,Longitude,Elevation\n"
                     f"{names[0]},-31.9,20.3,555\n"
                     f"{names[2].upper()},-31.8,20.4,556\nnope,0,0,0\n")
        assert dlg.import_coordinates(p) == 2
        ch = dlg.compute_changes()
        assert ch[names[0]]["elev"] == 555.0
        assert ch[names[2]]["lat"] == -31.8

    def test_lines_auto_and_rename(self, dlg):
        dlg._auto_lines()
        dlg._line_from.setCurrentText(dlg._line_from.itemText(0))
        dlg._line_to.setText("Profile A")
        dlg._rename_line()
        lines = {d.get("line") for d in dlg.compute_changes().values()}
        assert "Profile A" in lines

    def test_header_fill_all(self, dlg):
        dlg._fill_field.setCurrentIndex(dlg._fill_field.findData("project"))
        dlg._fill_value.setText("Kaapvaal")
        dlg._fill_header()
        ch = dlg.compute_changes()
        assert len(ch) == len(dlg._base)
        assert all(d == {"project": "Kaapvaal"} for d in ch.values())
        dlg.changes = ch
        assert dlg.change_label() == "Edit stations: header"


# ── main window round trip ──────────────────────────────────────────────────


@pytest.fixture
def window(qapp, monkeypatch, kap03):
    from pycsamt.app.desktop.main_window import MainWindow
    from pycsamt.app.desktop.models.session import SessionState

    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    win = MainWindow()
    lines = {s.name: "Survey" for s in kap03}
    assert win._adopt_full_dataset(kap03, record=False, lines=lines,
                                   merge=False)
    win._history.reset(win._survey_state())
    yield win
    win._history.mark_saved()
    win.close()


class TestMainWindow:
    def test_rename_keeps_line_and_undoes(self, window):
        sites, lines = window._survey_state()
        name = next(iter(sites)).name
        assert window.apply_station_changes(
            {name: {"name": "NEW-1", "elev": 999.0}}, "Edit stations: names")
        df = window._all_dataframe
        row = df[df["ID"] == "NEW-1"].iloc[0]
        assert row["Line"] == "Survey" and row["Elevation"] == 999.0
        assert window._act_undo.text() == "&Undo Edit stations: names"
        window.undo_edit()
        df = window._all_dataframe
        assert name in set(df["ID"]) and "NEW-1" not in set(df["ID"])
        assert df[df["ID"] == name].iloc[0]["Line"] == "Survey"

    def test_line_edit_is_undoable(self, window):
        sites, _ = window._survey_state()
        name = next(iter(sites)).name
        window.apply_station_changes({name: {"line": "L9"}},
                                     "Edit stations: lines")
        df = window._all_dataframe
        assert df[df["ID"] == name].iloc[0]["Line"] == "L9"
        window.undo_edit()
        df = window._all_dataframe
        assert df[df["ID"] == name].iloc[0]["Line"] == "Survey"

    def test_invalid_change_is_reported_not_applied(self, window):
        sites, _ = window._survey_state()
        name = next(iter(sites)).name
        assert not window.apply_station_changes({name: {"lat": 200.0}})
        assert "ERROR" in window._log_panel.text()
        assert not window._history.can_undo

    def test_menu_entries(self, window):
        texts = [a.text() for a in window._edit_station_menu.actions()]
        assert "Edit &Coordinates…" in texts and "&Rename Stations…" in texts
        assert "Assign &Lines…" in texts

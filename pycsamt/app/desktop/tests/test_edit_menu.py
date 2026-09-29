# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the Edit menu: survey edit history (undo / redo / history /
revert), Find Station, Save Edited Survey, the menu layout (Help at the
right) and shortcut uniqueness.  Real kap03 EDIs drive the survey tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.edit_history import EditHistory

KAP03 = Path(__file__).resolve().parents[4] / "data" / "MT" / "kap03lmt_edis"


# ── EditHistory (Qt-free) ───────────────────────────────────────────────────


class TestEditHistory:
    def test_undo_redo_and_labels(self):
        h = EditHistory()
        h.reset([0])
        h.record("A", [0], [1])
        h.record("B", [1], [2])
        assert h.undo_label == "B" and h.dirty
        assert h.undo() == [1] and h.redo_label == "B"
        assert h.undo() == [0] and not h.can_undo
        assert h.redo() == [1]
        assert [e[2] for e in h.entries()] == [True, False]

    def test_new_step_drops_the_redo_branch(self):
        h = EditHistory()
        h.reset([0])
        h.record("A", [0], [1])
        h.record("B", [1], [2])
        h.undo()
        h.record("C", [1], [3])
        assert [e[1] for e in h.entries()] == ["A", "C"]
        assert not h.can_redo

    def test_jump_original_and_saved(self):
        h = EditHistory()
        h.reset([0])
        for i in range(1, 4):
            h.record(f"s{i}", [i - 1], [i])
        assert h.jump(1) == [1] and h.position == 1
        assert h.jump(0) == [0]
        assert h.original() == [0]
        h.jump(3)
        h.mark_saved()
        assert not h.dirty
        h.undo()
        assert h.dirty

    def test_snapshots_are_independent(self):
        h = EditHistory()
        state = {"z": [1, 2]}
        h.reset(state)
        new = {"z": [1, 2, 3]}
        h.record("edit", state, new)
        state["z"].append(99)  # a tool editing in place afterwards
        assert h.undo() == {"z": [1, 2]}

    def test_bounded(self):
        h = EditHistory(max_steps=3)
        h.reset(0)
        for i in range(1, 6):
            h.record(str(i), i - 1, i)
        assert len(h) == 3 and h.position == 3
        assert [e[1] for e in h.entries()] == ["3", "4", "5"]


# ── main window ─────────────────────────────────────────────────────────────


@pytest.fixture
def window(qapp, monkeypatch):
    from pycsamt.app.desktop.models.session import SessionState

    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    from pycsamt.app.desktop.main_window import MainWindow

    win = MainWindow()
    yield win
    win._history.mark_saved()
    win.close()


@pytest.fixture(scope="module")
def kap03():
    if not KAP03.is_dir():
        pytest.skip("kap03 EDIs missing")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(KAP03))


@pytest.fixture
def loaded(window, kap03):
    window._history.reset(kap03)
    assert window._adopt_full_dataset(kap03, record=False)
    return window


def _menus(window) -> list[str]:
    return [a.text().replace("&", "") for a in window.menuBar().actions()]


def _menu_of(action):
    try:
        return action.menu()
    except RuntimeError:  # a rebuilt, already deleted menu
        return None


def _all_actions(menu):
    out = []
    for a in menu.actions():
        try:
            sub = a.menu()
            if sub is not None:
                out += _all_actions(sub)
            elif not a.isSeparator():
                out.append(a)
        except RuntimeError:  # a rebuilt submenu (Recent Files)
            continue
    return out


class TestMenuLayout:
    def test_edit_after_file_and_help_at_the_right(self, window):
        from PySide6.QtCore import Qt

        assert _menus(window)[:4] == ["File", "Edit", "View", "Tools"]
        assert "Help" not in _menus(window)
        corner = window.menuBar().cornerWidget(Qt.Corner.TopRightCorner)
        assert [a.text() for a in corner.actions()] == ["&Help"]

    def test_edit_menu_contents(self, window):
        texts = [a.text().replace("&", "").split("\t")[0]
                 for a in _all_actions(window._edit_menu)]
        for want in ("Undo", "Redo", "History…", "Revert to As Loaded",
                     "Find Station…", "Frequency Editor…", "Recompute EDIs…",
                     "Elevation Enrichment…", "Preferences…"):
            assert want in texts

    def test_editors_left_tools(self, window):
        # keep the action list alive: a menu reached through a temporary
        # list of actions is deleted with it (PySide ownership quirk)
        bar_actions = window.menuBar().actions()
        tools = next(_menu_of(a) for a in bar_actions if a.text() == "&Tools")
        texts = [a.text() for a in _all_actions(tools)]
        assert not any("Frequency Editor" in t or "Recompute" in t
                       or "Elevation" in t for t in texts)

    def test_no_two_actions_share_a_shortcut(self, window):
        from PySide6.QtCore import Qt

        corner = window.menuBar().cornerWidget(Qt.Corner.TopRightCorner)
        bar_actions = [*window.menuBar().actions(), *corner.actions()]
        menus = [m for a in bar_actions if (m := _menu_of(a)) is not None]
        seen: dict[str, str] = {}
        for m in menus:
            for a in _all_actions(m):
                for k in a.shortcuts():
                    key = k.toString()
                    if not key:
                        continue
                    assert key not in seen, (key, seen[key], a.text())
                    seen[key] = a.text()

    def test_nothing_to_undo_on_start(self, window):
        assert not window._act_undo.isEnabled()
        assert not window._act_save_survey.isEnabled()
        assert window.windowTitle() == "pycsamt"


class TestSurveyEdits:
    def _z(self, window):
        return np.asarray(next(iter(window._all_sites)).z)

    def test_edit_is_undoable_and_redoable(self, loaded):
        from pycsamt.site.edit import rotate_all

        z0 = self._z(loaded).copy()
        loaded._apply_modified_sites(rotate_all(loaded._all_sites, 30.0),
                                     source="Rotate 30°")
        z1 = self._z(loaded).copy()
        assert not np.allclose(z0, z1, equal_nan=True)
        assert loaded._act_undo.text() == "&Undo Rotate 30°"
        assert "● edited" in loaded.windowTitle()
        loaded.undo_edit()
        assert np.allclose(self._z(loaded), z0, equal_nan=True)
        assert loaded._act_redo.isEnabled()
        loaded.redo_edit()
        assert np.allclose(self._z(loaded), z1, equal_nan=True)

    def test_history_jump_and_revert(self, loaded):
        from pycsamt.site.edit import rotate_all

        z0 = self._z(loaded).copy()
        for ang in (10.0, 20.0):
            loaded._apply_modified_sites(
                rotate_all(loaded._all_sites, ang), source=f"Rotate {ang}")
        loaded.jump_to_history(0)
        assert np.allclose(self._z(loaded), z0, equal_nan=True)
        loaded.jump_to_history(2)
        loaded.revert_to_loaded()
        assert np.allclose(self._z(loaded), z0, equal_nan=True)
        assert loaded._history.undo_label == "Revert to as loaded"

    def test_failed_edit_is_reported(self, loaded):
        loaded._apply_modified_sites([], source="Broken tool")
        text = loaded._log_panel.text()
        assert "Broken tool" in text and "nothing was applied" in text
        assert not loaded._history.can_undo

    def test_find_station(self, loaded):
        name = loaded._station_ids()[3]
        assert loaded.find_station(name.upper() if name.isalpha()
                                   else name)
        assert loaded._controller.selected_station == name
        assert not loaded.find_station("no-such-station")

    def test_save_edited_survey(self, loaded, tmp_path):
        from pycsamt.site.edit import rotate_all

        loaded._apply_modified_sites(rotate_all(loaded._all_sites, 15.0),
                                     source="Rotate")
        paths = loaded.save_edited_survey(tmp_path / "edi")
        assert len(paths) == len(loaded._station_ids())
        assert not loaded._history.dirty
        assert "● edited" not in loaded.windowTitle()
        xml = loaded.save_edited_survey(tmp_path / "xml", fmt="xml")
        assert all(Path(p).suffix == ".xml" for p in xml)


def test_open_menu_says_load_data(window):
    assert window._act_open.text() == "&Open / Load Data…"


def test_toolbar_overflow_button_follows_the_theme(window):
    from PySide6.QtGui import QColor
    from PySide6.QtWidgets import QToolButton

    btn = window._main_toolbar.findChild(QToolButton,
                                         "qt_toolbar_ext_button")
    assert btn is not None

    def dot_colour():
        img = btn.icon().pixmap(20, 20).toImage()
        return QColor(img.pixel(img.width() // 2, img.height() // 2))

    window._apply_theme("dark")
    assert dot_colour().lightness() > 150  # light dots on the dark theme
    window._apply_theme("light")
    assert dot_colour().lightness() < 120
    assert btn.toolTip() == "More tools"

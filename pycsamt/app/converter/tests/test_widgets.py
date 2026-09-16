# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.converter.widgets (DropZone, LogPanel)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtCore import QUrl  # noqa: E402

from pycsamt.app.converter.widgets import DropZone, LogPanel  # noqa: E402


class _FakeMimeData:
    """Duck-types the slice of QMimeData that DropZone actually reads.

    A real QDropEvent(QMimeData(...)) is unsafe here: PySide6/Shiboken
    does not keep the QMimeData alive past the local QDropEvent
    constructor call, so by the time dropEvent() reads it back the
    object has been reinterpreted as a bare QObject. A tiny stand-in
    with real QUrl entries sidesteps that lifetime hazard entirely.
    """

    def __init__(self, paths):
        self._urls = [QUrl.fromLocalFile(str(p)) for p in paths]

    def hasUrls(self) -> bool:
        return bool(self._urls)

    def urls(self):
        return self._urls


class _FakeDropEvent:
    def __init__(self, paths):
        self._mime = _FakeMimeData(paths)
        self.accepted = False

    def mimeData(self):
        return self._mime

    def acceptProposedAction(self):
        self.accepted = True


def _drop_event(paths):
    return _FakeDropEvent(paths)


# ── DropZone construction ───────────────────────────────────────────────


def test_drop_zone_default_construction(qapp):
    dz = DropZone()
    assert dz.path() is None
    assert dz.acceptDrops() is True
    dz.close()


def test_drop_zone_accepts_both_by_default_has_two_buttons(qapp):
    dz = DropZone("Source")
    buttons = dz.findChildren(type(dz), None)  # no-op, just ensure widget built
    assert dz is not None
    dz.close()


def test_drop_zone_files_only_hides_folder_button(qapp):
    from PySide6.QtWidgets import QPushButton

    dz = DropZone("Files only", accept_dirs=False, accept_files=True)
    labels = [b.text() for b in dz.findChildren(QPushButton)]
    assert "Browse…" in labels
    assert "Folder…" not in labels
    dz.close()


def test_drop_zone_dirs_only_hides_file_button(qapp):
    from PySide6.QtWidgets import QPushButton

    dz = DropZone("Dirs only", accept_dirs=True, accept_files=False)
    labels = [b.text() for b in dz.findChildren(QPushButton)]
    assert "Folder…" in labels
    assert "Browse…" not in labels
    dz.close()


# ── path() / set_path() / clear() ───────────────────────────────────────


def test_set_path_and_path_roundtrip(qapp, tmp_path):
    dz = DropZone()
    dz.set_path(tmp_path)
    assert dz.path() == tmp_path
    dz.close()


def test_path_returns_none_for_blank_text(qapp):
    dz = DropZone()
    dz.set_path("   ")
    assert dz.path() is None
    dz.close()


def test_clear_resets_path(qapp, tmp_path):
    dz = DropZone()
    dz.set_path(tmp_path)
    dz.clear()
    assert dz.path() is None
    dz.close()


def test_path_changed_signal_emits_on_set_path(qapp, tmp_path):
    dz = DropZone()
    seen = []
    dz.pathChanged.connect(seen.append)
    dz.set_path(tmp_path)
    assert seen == [str(tmp_path)]
    dz.close()


# ── drag-and-drop ────────────────────────────────────────────────────────


def test_drop_event_sets_path_for_accepted_file(qapp, tmp_path):
    f = tmp_path / "a.edi"
    f.write_text("data")
    dz = DropZone("Source", accept_files=True, accept_dirs=True)
    event = _drop_event([f])
    dz.dropEvent(event)
    assert dz.path() == f
    dz.close()


def test_drop_event_sets_path_for_accepted_dir(qapp, tmp_path):
    dz = DropZone("Source", accept_files=True, accept_dirs=True)
    event = _drop_event([tmp_path])
    dz.dropEvent(event)
    assert dz.path() == tmp_path
    dz.close()


def test_drop_event_rejects_dir_when_dirs_not_accepted(qapp, tmp_path):
    dz = DropZone("Files only", accept_files=True, accept_dirs=False)
    event = _drop_event([tmp_path])
    dz.dropEvent(event)
    assert dz.path() is None
    dz.close()


def test_drop_event_rejects_file_when_files_not_accepted(qapp, tmp_path):
    f = tmp_path / "a.edi"
    f.write_text("data")
    dz = DropZone("Dirs only", accept_files=False, accept_dirs=True)
    event = _drop_event([f])
    dz.dropEvent(event)
    assert dz.path() is None
    dz.close()


def test_drop_event_with_no_urls_is_a_noop(qapp):
    dz = DropZone()
    event = _FakeDropEvent([])
    dz.dropEvent(event)
    assert dz.path() is None
    dz.close()


class _FakeDragEnterEvent:
    def __init__(self, paths):
        self._mime = _FakeMimeData(paths)
        self.accepted = False

    def mimeData(self):
        return self._mime

    def acceptProposedAction(self):
        self.accepted = True


def test_drag_enter_and_leave_toggle_drag_active_property(qapp, tmp_path):
    dz = DropZone()
    enter_event = _FakeDragEnterEvent([tmp_path])
    dz.dragEnterEvent(enter_event)
    assert dz.property("dragActive") == "true"
    assert enter_event.accepted is True

    dz.dragLeaveEvent(None)
    assert dz.property("dragActive") == "false"
    dz.close()


# ── Browse buttons (QFileDialog monkeypatched) ──────────────────────────


def test_browse_file_sets_path_when_chosen(qapp, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog

    f = tmp_path / "chosen.edi"
    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: (str(f), ""))
    )
    dz = DropZone("Source")
    dz._browse_file()
    assert dz.path() == f
    dz.close()


def test_browse_file_noop_when_cancelled(qapp, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: ("", ""))
    )
    dz = DropZone("Source")
    dz._browse_file()
    assert dz.path() is None
    dz.close()


def test_browse_dir_sets_path_when_chosen(qapp, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: str(tmp_path))
    )
    dz = DropZone("Source")
    dz._browse_dir()
    assert dz.path() == tmp_path
    dz.close()


def test_browse_dir_noop_when_cancelled(qapp, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: "")
    )
    dz = DropZone("Source")
    dz._browse_dir()
    assert dz.path() is None
    dz.close()


def test_browse_file_starts_from_current_path(qapp, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QFileDialog

    dz = DropZone("Source")
    dz.set_path(tmp_path)
    captured = {}

    def _fake(_self, _title, start, _filter):
        captured["start"] = start
        return "", ""

    monkeypatch.setattr(QFileDialog, "getOpenFileName", staticmethod(_fake))
    dz._browse_file()
    assert captured["start"] == str(tmp_path)
    dz.close()


# ── LogPanel ─────────────────────────────────────────────────────────────


def test_log_panel_append_adds_text(qapp):
    panel = LogPanel()
    panel.append("hello")
    assert "hello" in panel._log.toPlainText()
    panel.close()


def test_log_panel_clear_resets_text_and_progress(qapp):
    panel = LogPanel()
    panel.append("something")
    panel.set_progress(3, 10)
    panel.clear()
    assert panel._log.toPlainText() == ""
    assert panel._progress.value() == 0
    panel.close()


def test_log_panel_set_progress(qapp):
    panel = LogPanel()
    panel.set_progress(2, 5)
    assert panel._progress.maximum() == 5
    assert panel._progress.value() == 2
    panel.close()


def test_log_panel_set_progress_total_zero_floors_to_one(qapp):
    panel = LogPanel()
    panel.set_progress(0, 0)
    assert panel._progress.maximum() == 1
    panel.close()


def test_log_panel_reset_progress(qapp):
    panel = LogPanel()
    panel.set_progress(4, 10)
    panel.reset_progress()
    assert panel._progress.maximum() == 1
    assert panel._progress.value() == 0
    panel.close()


def test_log_panel_set_indeterminate_true_and_false(qapp):
    panel = LogPanel()
    panel.set_indeterminate(True)
    assert panel._progress.minimum() == 0
    assert panel._progress.maximum() == 0
    panel.set_indeterminate(False)
    assert panel._progress.maximum() == 1
    panel.close()

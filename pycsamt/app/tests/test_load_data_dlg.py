# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for LoadDataDialog and its _DropZone helper widget."""

from __future__ import annotations

from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QDialog

from pycsamt.app.desktop.dialogs.load_data_dlg import (
    LoadDataDialog,
    _DropZone,
)


@pytest.fixture
def dlg(qapp):
    d = LoadDataDialog()
    yield d
    d.close()


def _fake_drop_event(paths, local=True):
    urls = []
    for p in paths:
        u = mock.Mock()
        u.isLocalFile.return_value = local
        u.toLocalFile.return_value = p
        urls.append(u)
    mime = mock.Mock()
    mime.hasUrls.return_value = bool(urls)
    mime.urls.return_value = urls
    event = mock.Mock()
    event.mimeData.return_value = mime
    return event


# ── _DropZone ─────────────────────────────────────────────────────────────


def test_drop_zone_initial_text(qapp):
    dz = _DropZone()
    assert dz.text() == _DropZone._TEXT_IDLE


def test_drag_enter_with_urls_sets_hover(qapp):
    dz = _DropZone()
    event = _fake_drop_event(["/a/b.edi"])
    dz.dragEnterEvent(event)
    event.acceptProposedAction.assert_called_once()
    assert dz.text() == _DropZone._TEXT_HOVER
    assert dz.property("drag_over") == "true"


def test_drag_enter_without_urls_noop(qapp):
    dz = _DropZone()
    event = _fake_drop_event([])
    dz.dragEnterEvent(event)
    event.acceptProposedAction.assert_not_called()
    assert dz.text() == _DropZone._TEXT_IDLE


def test_drag_leave_resets_hover(qapp):
    dz = _DropZone()
    dz._set_drag_over(True)
    dz.dragLeaveEvent(mock.Mock())
    assert dz.text() == _DropZone._TEXT_IDLE
    assert dz.property("drag_over") == "false"


def test_drop_event_emits_local_paths(qapp):
    dz = _DropZone()
    received = []
    dz.raw_paths_dropped.connect(received.append)
    event = _fake_drop_event(["/a/b.edi", "/c/d.edi"])
    dz.dropEvent(event)
    assert received == [["/a/b.edi", "/c/d.edi"]]
    event.acceptProposedAction.assert_called_once()


def test_drop_event_no_local_files_no_emit(qapp):
    dz = _DropZone()
    received = []
    dz.raw_paths_dropped.connect(received.append)
    event = _fake_drop_event(["remote://x"], local=False)
    dz.dropEvent(event)
    assert received == []


# ── LoadDataDialog construction ──────────────────────────────────────────


def test_creates(dlg):
    assert dlg.windowTitle() == "Open Survey Data"


def test_no_recomputed_button_by_default(dlg):
    from PySide6.QtWidgets import QPushButton

    labels = [b.text() for b in dlg.findChildren(QPushButton)]
    assert not any("Recomputed" in t for t in labels)


def test_recomputed_button_shown_when_dir_exists(qapp, tmp_path):
    from PySide6.QtWidgets import QPushButton

    d = LoadDataDialog(recomputed_dir=tmp_path)
    labels = [b.text() for b in d.findChildren(QPushButton)]
    assert any("Recomputed" in t for t in labels)
    d.close()


def test_recomputed_button_hidden_when_dir_missing(qapp, tmp_path):
    from PySide6.QtWidgets import QPushButton

    d = LoadDataDialog(recomputed_dir=tmp_path / "nope")
    labels = [b.text() for b in d.findChildren(QPushButton)]
    assert not any("Recomputed" in t for t in labels)
    d.close()


def test_ok_button_initially_disabled(dlg):
    assert not dlg._ok_btn.isEnabled()


def test_remove_clear_buttons_initially_disabled(dlg):
    assert not dlg._btn_remove.isEnabled()
    assert not dlg._btn_clear.isEnabled()


# ── _on_dropped ───────────────────────────────────────────────────────────


def test_on_dropped_file_matching_format(dlg, tmp_path):
    f = tmp_path / "A.edi"
    f.write_text("x")
    dlg._on_dropped([str(f)])
    assert dlg._file_list.count() == 1
    assert dlg._file_list.item(0).text() == str(f)


def test_on_dropped_file_not_matching_format_ignored(dlg, tmp_path):
    f = tmp_path / "A.avg"
    f.write_text("x")
    dlg._on_dropped([str(f)])
    assert dlg._file_list.count() == 0


def test_on_dropped_folder_recurses(dlg, tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "S1.edi").write_text("x")
    (tmp_path / "S2.edi").write_text("x")
    (tmp_path / "ignore.txt").write_text("x")
    dlg._on_dropped([str(tmp_path)])
    assert dlg._file_list.count() == 2


def test_on_dropped_nothing_matched_shows_warning_text(dlg, tmp_path):
    f = tmp_path / "A.txt"
    f.write_text("x")
    dlg._on_dropped([str(f)])
    assert "No" in dlg._drop_zone.text()


def test_on_dropped_empty_list_noop(dlg):
    dlg._on_dropped([])
    assert dlg._file_list.count() == 0


# ── _load_recomputed ──────────────────────────────────────────────────────


def test_load_recomputed_populates_list(qapp, tmp_path):
    (tmp_path / "S1.edi").write_text("x")
    (tmp_path / "S2.edi").write_text("x")
    d = LoadDataDialog(recomputed_dir=tmp_path)
    d._load_recomputed()
    assert d._file_list.count() == 2
    d.close()


def test_load_recomputed_empty_dir_shows_warning(qapp, tmp_path):
    d = LoadDataDialog(recomputed_dir=tmp_path)
    d._load_recomputed()
    assert "No EDI files" in d._drop_zone.text()
    d.close()


def test_load_recomputed_no_dir_noop(dlg):
    dlg._load_recomputed()  # _recomputed_dir is None -> early return, no raise


# ── Browse slots ──────────────────────────────────────────────────────────


def test_browse_files_sets_paths(dlg, tmp_path, monkeypatch):
    f1 = tmp_path / "A.edi"
    f2 = tmp_path / "B.edi"
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileNames",
        staticmethod(lambda *a, **k: ([str(f1), str(f2)], "")),
    )
    dlg._browse_files()
    assert dlg._file_list.count() == 2
    assert dlg._last_dir == str(tmp_path)


def test_browse_files_cancelled_noop(dlg, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileNames",
        staticmethod(lambda *a, **k: ([], "")),
    )
    dlg._browse_files()
    assert dlg._file_list.count() == 0


def test_browse_folder_sets_paths(dlg, tmp_path, monkeypatch):
    (tmp_path / "A.edi").write_text("x")
    (tmp_path / "B.edi").write_text("x")
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    dlg._browse_folder()
    assert dlg._file_list.count() == 2
    assert dlg._last_dir == str(tmp_path)


def test_browse_folder_cancelled_noop(dlg, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: ""),
    )
    dlg._browse_folder()
    assert dlg._file_list.count() == 0


# ── File-list mutations ────────────────────────────────────────────────────


def test_add_paths_skips_duplicates(dlg):
    dlg._add_paths(["/a.edi", "/b.edi"])
    dlg._add_paths(["/b.edi", "/c.edi"])
    texts = [dlg._file_list.item(i).text() for i in range(dlg._file_list.count())]
    assert texts == ["/a.edi", "/b.edi", "/c.edi"]


def test_remove_selected(dlg):
    dlg._set_paths(["/a.edi", "/b.edi", "/c.edi"])
    dlg._file_list.item(1).setSelected(True)
    dlg._remove_selected()
    texts = [dlg._file_list.item(i).text() for i in range(dlg._file_list.count())]
    assert texts == ["/a.edi", "/c.edi"]


def test_clear_all_resets_dropzone_text(dlg):
    dlg._set_paths(["/a.edi"])
    dlg._drop_zone.setText("something else")
    dlg._clear_all()
    assert dlg._file_list.count() == 0
    assert dlg._drop_zone.text() == _DropZone._TEXT_IDLE


def test_refresh_ui_updates_count_label_and_buttons(dlg):
    dlg._set_paths(["/a.edi", "/b.edi"])
    assert dlg._count_lbl.text() == "Selected files: 2"
    assert dlg._btn_clear.isEnabled()
    assert dlg._ok_btn.isEnabled()


def test_refresh_ui_zero_files_disables_buttons(dlg):
    dlg._set_paths([])
    assert dlg._count_lbl.text() == "Selected files: 0"
    assert not dlg._btn_clear.isEnabled()
    assert not dlg._ok_btn.isEnabled()


def test_on_sel_changed_enables_remove(dlg):
    dlg._set_paths(["/a.edi"])
    dlg._file_list.item(0).setSelected(True)
    dlg._on_sel_changed()
    assert dlg._btn_remove.isEnabled()


# ── Accept ────────────────────────────────────────────────────────────────


def test_on_accepted_sets_selected_paths(dlg):
    dlg._set_paths(["/a.edi", "/b.edi"])
    dlg._on_accepted()
    assert dlg.selected_paths == ["/a.edi", "/b.edi"]
    assert dlg.result() == QDialog.DialogCode.Accepted


def test_ok_button_click_triggers_accept_flow(dlg):
    dlg._set_paths(["/a.edi"])
    dlg._ok_btn.click()
    assert dlg.selected_paths == ["/a.edi"]
    assert dlg.result() == QDialog.DialogCode.Accepted

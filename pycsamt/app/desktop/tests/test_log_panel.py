# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for LogPanel (Phase 1) — runs headless via Qt offscreen platform."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")


@pytest.fixture
def log_panel(qapp):
    from pycsamt.app.desktop.panels.log_panel import LogPanel

    panel = LogPanel()
    yield panel
    panel.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_log_panel_creates_without_error(qapp):
    from pycsamt.app.desktop.panels.log_panel import LogPanel

    panel = LogPanel()
    assert panel is not None
    panel.close()


# ── append_line ───────────────────────────────────────────────────────────


def test_append_line_adds_text(log_panel):
    log_panel.append_line("hello world")
    text = log_panel._text.toPlainText()
    assert "hello world" in text


def test_append_line_adds_timestamp(log_panel):
    log_panel.append_line("timestamped")
    text = log_panel._text.toPlainText()
    # Timestamp format is [HH:MM:SS]
    assert "[" in text and "]" in text


def test_multiple_lines_are_all_present(log_panel):
    for i in range(5):
        log_panel.append_line(f"line {i}")
    text = log_panel._text.toPlainText()
    for i in range(5):
        assert f"line {i}" in text


# ── clear ─────────────────────────────────────────────────────────────────


def test_clear_removes_all_text(log_panel):
    log_panel.append_line("will be cleared")
    log_panel.clear()
    assert log_panel._text.toPlainText() == ""


# ── block cap ─────────────────────────────────────────────────────────────


def test_max_block_count_is_set(log_panel):
    assert log_panel._text.maximumBlockCount() == 2000


# ── read-only ─────────────────────────────────────────────────────────────


def test_text_widget_is_read_only(log_panel):
    assert log_panel._text.isReadOnly()


# ── levels, filter, theme (v2.6) ───────────────────────────────────────────


@pytest.mark.parametrize("text, level", [
    ("ERROR: disk full", "error"), ("Load failed: bad header", "error"),
    ("✕  no model", "error"), ("Warning: 3 stations skipped", "warn"),
    ("Could not read x.edi", "warn"), ("Session saved.", "ok"),
    ("✓  done", "ok"), ("pycsamt ready — load EDI files", "info"),
])
def test_level_is_read_from_the_text(text, level):
    from pycsamt.app.desktop.panels.log_panel import classify

    assert classify(text) == level


def test_explicit_level_and_signal(log_panel):
    got = []
    log_panel.message_logged.connect(lambda lv, t: got.append((lv, t)))
    log_panel.append_line("plain words", level="warn")
    assert got == [("warn", "plain words")]
    assert log_panel.records()[-1][1] == "warn"


def test_level_filter_and_find(log_panel):
    for t in ("loaded 3 sites", "Warning: gap", "ERROR: boom", "other"):
        log_panel.append_line(t)
    log_panel._level.setCurrentIndex(log_panel._level.findData("error"))
    assert log_panel.text().count("\n") == 0 and "boom" in log_panel.text()
    log_panel._level.setCurrentIndex(log_panel._level.findData("warn"))
    assert "gap" in log_panel.text() and "other" not in log_panel.text()
    log_panel._level.setCurrentIndex(0)
    log_panel._find.setText("sites")
    assert log_panel.text().strip().endswith("loaded 3 sites")
    assert "2 lines" not in log_panel._counts.text()
    assert "1 error" in log_panel._counts.text()


def test_dark_and_light_colours(log_panel):
    log_panel.append_line("ERROR: boom")

    def err_colour():
        block = log_panel._text.document().lastBlock()
        it = block.begin()
        colours = []
        while not it.atEnd():
            colours.append(it.fragment().charFormat().foreground().color()
                           .name())
            it += 1
        return colours[-1]

    assert err_colour() == "#b42318"
    log_panel.set_dark(True)
    assert err_colour() == "#ff7b72"
    assert "#161b22" in log_panel._text.styleSheet()
    log_panel.set_dark(False)
    assert "#ffffff" in log_panel._text.styleSheet()


def test_save_to(log_panel, tmp_path):
    log_panel.append_line("Session saved.")
    log_panel.append_line("ERROR: x")
    text = log_panel.save_to(tmp_path / "log.txt").read_text()
    assert "OK    Session saved." in text and "ERROR ERROR: x" in text

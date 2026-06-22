# -*- coding: utf-8 -*-
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

# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for SessionState (Phase 1)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pycsamt.app.desktop.models.session import SessionState


# ── Defaults ──────────────────────────────────────────────────────────────

def test_default_theme():
    s = SessionState()
    assert s.theme == "dark"


def test_default_recent_files_is_empty_list():
    s = SessionState()
    assert s.recent_files == []


def test_default_dock_fields_are_none():
    s = SessionState()
    assert s.dock_geometry is None
    assert s.dock_state is None


# ── Round-trip save / load ────────────────────────────────────────────────

def test_save_and_load_round_trip(tmp_path):
    p = tmp_path / "session.json"
    s = SessionState(
        theme="light",
        recent_files=["/data/survey1.edi", "/data/survey2.edi"],
        last_data_dir="/data",
        selected_station="S05",
        freq_min_hz=0.1,
        freq_max_hz=1000.0,
        overlay="Phase",
        dock_geometry="abc123==",
        dock_state="xyz456==",
    )
    s.save(p)
    s2 = SessionState.load(p)

    assert s2.theme == "light"
    assert s2.recent_files == ["/data/survey1.edi", "/data/survey2.edi"]
    assert s2.last_data_dir == "/data"
    assert s2.selected_station == "S05"
    assert s2.freq_min_hz == pytest.approx(0.1)
    assert s2.freq_max_hz == pytest.approx(1000.0)
    assert s2.overlay == "Phase"
    assert s2.dock_geometry == "abc123=="
    assert s2.dock_state == "xyz456=="


def test_load_returns_default_when_file_missing(tmp_path):
    s = SessionState.load(tmp_path / "no_such_file.json")
    assert s.theme == "dark"
    assert s.dock_geometry is None


def test_load_ignores_unknown_keys(tmp_path):
    p = tmp_path / "session.json"
    data = {"theme": "light", "unknown_key": "should_be_ignored"}
    p.write_text(json.dumps(data), encoding="utf-8")
    s = SessionState.load(p)
    assert s.theme == "light"
    assert not hasattr(s, "unknown_key")


def test_load_returns_default_on_corrupt_json(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{not valid json", encoding="utf-8")
    s = SessionState.load(p)
    assert s.theme == "dark"


def test_save_creates_parent_directories(tmp_path):
    p = tmp_path / "a" / "b" / "session.json"
    SessionState().save(p)
    assert p.exists()


# ── Mutability ────────────────────────────────────────────────────────────

def test_session_fields_are_mutable():
    s = SessionState()
    s.theme = "light"
    s.recent_files.append("/data/new.edi")
    assert s.theme == "light"
    assert "/data/new.edi" in s.recent_files


def test_two_instances_do_not_share_recent_files():
    s1 = SessionState()
    s2 = SessionState()
    s1.recent_files.append("/x.edi")
    assert "/x.edi" not in s2.recent_files

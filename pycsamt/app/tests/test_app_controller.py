# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AppController (Phase 2) — no Qt required."""

from __future__ import annotations

from pycsamt.app.desktop.controllers.app_controller import (
    AppController,
)
from pycsamt.app.desktop.models.session import SessionState

# ── Construction ──────────────────────────────────────────────────────────

def test_creates_with_default_session(monkeypatch):
    # Prevent reading the real on-disk session (which may have persisted state).
    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    ctrl = AppController()
    assert ctrl.session is not None
    assert ctrl.sites is None
    assert ctrl.selected_station is None


def test_creates_with_given_session():
    sess = SessionState(theme="light", selected_station="S07")
    ctrl = AppController(session=sess)
    assert ctrl.session.theme == "light"
    assert ctrl.selected_station == "S07"


# ── set_sites ────────────────────────────────────────────────────────────

def test_set_sites_stores_object():
    ctrl = AppController()
    sentinel = object()
    ctrl.set_sites(sentinel)
    assert ctrl.sites is sentinel


def test_set_sites_fires_callbacks():
    ctrl = AppController()
    received = []
    ctrl.on_data_loaded(received.append)
    sentinel = object()
    ctrl.set_sites(sentinel)
    assert received == [sentinel]


def test_set_sites_fires_multiple_callbacks():
    ctrl = AppController()
    log1, log2 = [], []
    ctrl.on_data_loaded(log1.append)
    ctrl.on_data_loaded(log2.append)
    ctrl.set_sites("x")
    assert log1 == ["x"]
    assert log2 == ["x"]


def test_set_sites_ignores_callback_exception():
    ctrl = AppController()
    def bad_cb(sites):
        raise RuntimeError("oops")
    ctrl.on_data_loaded(bad_cb)
    ctrl.set_sites(object())   # must not raise


# ── select_station ────────────────────────────────────────────────────────

def test_select_station_stores_id():
    ctrl = AppController()
    ctrl.select_station("S03")
    assert ctrl.selected_station == "S03"


def test_select_station_updates_session():
    ctrl = AppController()
    ctrl.select_station("S10")
    assert ctrl.session.selected_station == "S10"


def test_select_station_fires_callbacks():
    ctrl = AppController()
    received = []
    ctrl.on_station_selected(received.append)
    ctrl.select_station("A99")
    assert received == ["A99"]


# ── post_status ───────────────────────────────────────────────────────────

def test_post_status_fires_callbacks():
    ctrl = AppController()
    msgs = []
    ctrl.on_status_message(msgs.append)
    ctrl.post_status("hello")
    assert msgs == ["hello"]


# ── add_recent_file ───────────────────────────────────────────────────────

def test_add_recent_file_prepends():
    ctrl = AppController()
    ctrl.add_recent_file("/a.edi")
    ctrl.add_recent_file("/b.edi")
    assert ctrl.session.recent_files[0] == "/b.edi"


def test_add_recent_file_deduplicates():
    ctrl = AppController()
    ctrl.add_recent_file("/a.edi")
    ctrl.add_recent_file("/b.edi")
    ctrl.add_recent_file("/a.edi")
    assert ctrl.session.recent_files.count("/a.edi") == 1
    assert ctrl.session.recent_files[0] == "/a.edi"


def test_add_recent_file_caps_at_20():
    ctrl = AppController()
    for i in range(25):
        ctrl.add_recent_file(f"/{i}.edi")
    assert len(ctrl.session.recent_files) == 20


# ── n_stations ────────────────────────────────────────────────────────────

def test_n_stations_is_zero_when_no_data():
    assert AppController().n_stations == 0


def test_n_stations_reflects_sites_len():
    ctrl = AppController()
    ctrl.set_sites([1, 2, 3])   # a list counts as mock Sites
    assert ctrl.n_stations == 3

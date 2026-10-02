# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.recompute (EDI recompute workflow).

Real data
---------
Uses the bundled ``data/3edis`` EDI folder (skipped if absent) to
exercise ``_worker``/``EDIRecomputer`` end to end rather than mocking
the recompute pipeline itself. Callback functions registered via
``@app.callback(...)`` (not bare ``dash.callback``) are visible in
``web_app.callback_map`` immediately after ``create_app()`` — no
``_setup_server()`` dance needed here, unlike the global-callback
modules documented elsewhere in this test tree.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path

import pytest
from dash import no_update

from pycsamt.app.web.callbacks import recompute as recompute_mod
from pycsamt.app.web.layout import IDs

_PROJECT_ROOT = Path(__file__).resolve().parents[4]
_EDI_DIR = _PROJECT_ROOT / "data" / "3edis"
_HAS_EDIS = _EDI_DIR.exists() and any(_EDI_DIR.glob("*.edi"))


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_by_input(web_app, output_substr, input_id_substr):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(
            input_id_substr in str(i.get("id"))
            for i in v.get("inputs", [])
        ):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} "
        f"input~={input_id_substr!r}"
    )


def _set_triggered(prop_id):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(
        AttributeDict(triggered_inputs=[{"prop_id": prop_id}])
    )


def _clear_triggered():
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    cc.context_value.set(AttributeDict(triggered_inputs=[]))


@pytest.fixture(autouse=True)
def _reset_jobs():
    with recompute_mod._LOCK:
        recompute_mod._JOBS.clear()
    yield
    with recompute_mod._LOCK:
        recompute_mod._JOBS.clear()


@pytest.fixture
def sync_thread(monkeypatch):
    """Run the worker synchronously instead of on a background thread."""

    class _SyncThread:
        def __init__(self, target=None, args=(), daemon=None):
            self._target = target
            self._args = args

        def start(self):
            self._target(*self._args)

    monkeypatch.setattr(recompute_mod.threading, "Thread", _SyncThread)


# ─────────────────────────────────────────────────────────────────────────────
# pure helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_collect_edis_finds_real_files():
    if not _HAS_EDIS:
        pytest.skip("3edis dataset not found")
    paths = recompute_mod._collect_edis(str(_EDI_DIR))
    assert len(paths) > 0
    assert all(p.lower().endswith(".edi") for p in paths)
    assert paths == sorted(paths)


def test_collect_edis_empty_folder(tmp_path):
    paths = recompute_mod._collect_edis(str(tmp_path))
    assert paths == []


def test_safe_name_fallback_on_exception():
    class _Bad:
        def __getattr__(self, item):
            raise RuntimeError("boom")

    name = recompute_mod._safe_name(_Bad(), 3)
    assert name == "site_0004"


def test_safe_name_uses_station_name(monkeypatch):
    monkeypatch.setattr(
        "pycsamt.site.utils.station_name", lambda ed: "MYSTATION"
    )

    class _Edi:
        pass

    name = recompute_mod._safe_name(_Edi(), 0)
    assert name == "MYSTATION"


def test_log_chips_truncation_and_order():
    entries = [f"✓ s{i}" for i in range(20)]
    chips = recompute_mod._log_chips(entries)
    assert len(chips) == recompute_mod._LOG_SHOWN
    # newest entry first
    assert chips[0].children == "✓ s19"
    assert "recompute-log-ok" in chips[0].className


def test_log_chips_marks_failures():
    entries = ["✓ ok_one", "✗ bad_one: some error"]
    chips = recompute_mod._log_chips(entries)
    classes = {c.children: c.className for c in chips}
    assert "recompute-log-ok" in classes["✓ ok_one"]
    assert "recompute-log-err" in classes["✗ bad_one: some error"]


# ─────────────────────────────────────────────────────────────────────────────
# _sites_to_store
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def real_sites():
    if not _HAS_EDIS:
        pytest.skip("3edis dataset not found")
    from pycsamt.agents import MTLoaderAgent

    r = MTLoaderAgent().execute({"path": str(_EDI_DIR)})
    if r.status != "success":
        pytest.skip("Could not load 3edis sites.")
    return r["sites"]


def test_sites_to_store_basic_fields(real_sites):
    store = recompute_mod._sites_to_store(real_sites)
    assert store["recomputed"] is True
    assert store["data_dir"] == "[recomputed]"
    assert store["n_stations"] == len(store["station_records"])
    assert "n_lines" in store
    assert "line_counts" in store


def test_sites_to_store_merges_old_line_mapping(real_sites):
    old_ids = [
        r["ID"] for r in recompute_mod._sites_to_store(real_sites)[
            "station_records"
        ]
    ]
    old_store = {
        "station_records": [
            {"ID": sid, "Line": "L99"} for sid in old_ids
        ]
    }
    store = recompute_mod._sites_to_store(real_sites, old_store=old_store)
    assert store["n_stations"] == len(old_ids)


# ─────────────────────────────────────────────────────────────────────────────
# _worker
# ─────────────────────────────────────────────────────────────────────────────


def test_worker_path_source_completes(real_sites):
    session_id = "sess-worker-path"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "starting",
            "phase": "loading",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 0,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }

    recompute_mod._worker(session_id, str(_EDI_DIR), {})

    job = recompute_mod._JOBS[session_id]
    assert job["status"] == "done"
    assert job["sites"] is not None
    assert job["done"] == job["total"]
    assert job["total"] > 0


def test_worker_current_source_completes(real_sites):
    session_id = "sess-worker-current"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "starting",
            "phase": "recomputing",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 0,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }

    recompute_mod._worker(
        session_id,
        real_sites,
        {"comps": ["Z"], "resphase": ["on"]},
    )

    job = recompute_mod._JOBS[session_id]
    assert job["status"] == "done"
    assert job["sites"] is not None


def test_worker_path_with_no_edis_errors(tmp_path):
    session_id = "sess-worker-empty"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "starting",
            "phase": "loading",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 0,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }

    recompute_mod._worker(session_id, str(tmp_path), {})

    job = recompute_mod._JOBS[session_id]
    assert job["status"] == "error"
    assert "No EDI files found" in job["err_msg"]


def test_worker_outer_exception_sets_error_status(monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("kaboom")

    monkeypatch.setattr(recompute_mod, "Sites", None, raising=False)
    monkeypatch.setattr(
        "pycsamt.site.base.Sites", _boom
    )

    session_id = "sess-worker-boom"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "starting",
            "phase": "loading",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 0,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }

    if not _HAS_EDIS:
        pytest.skip("3edis dataset not found")

    recompute_mod._worker(session_id, str(_EDI_DIR), {})

    job = recompute_mod._JOBS[session_id]
    assert job["status"] == "error"
    assert "kaboom" in job["err_msg"]


# ─────────────────────────────────────────────────────────────────────────────
# toggle_canvas / cancel_confirm
# ─────────────────────────────────────────────────────────────────────────────


def test_toggle_canvas_flips_state(web_app):
    fn = _cb(web_app, f'{{"index":"{IDs.RECOMPUTE_CANVAS}","type":"?"}}.is_open') \
        if False else None
    fn = _cb_by_input(
        web_app, f"{IDs.RECOMPUTE_CANVAS}.is_open", IDs.BTN_RECOMPUTE_OPEN
    )
    assert fn(1, False) is True
    assert fn(1, True) is False
    assert fn(0, False) is no_update


def test_cancel_confirm_hides_on_click(web_app):
    fn = _cb_by_input(
        web_app,
        f"{IDs.RECOMPUTE_CONFIRM_SEC}.style",
        IDs.BTN_RECOMPUTE_CANCEL,
    )
    assert fn(1) == {"display": "none"}
    assert fn(0) is no_update


# ─────────────────────────────────────────────────────────────────────────────
# start_recompute
# ─────────────────────────────────────────────────────────────────────────────


def _start_recompute_fn(web_app):
    return _cb_by_input(
        web_app, f"{IDs.RECOMPUTE_PROG_SEC}.style", IDs.BTN_RECOMPUTE
    )


def test_start_recompute_no_trigger_skips(web_app):
    fn = _start_recompute_fn(web_app)
    _clear_triggered()
    result = fn(
        None, None, "current", "", [], None, [], None, None, None,
        "sess-1", None,
    )
    assert result == (no_update,) * 12


def test_start_recompute_no_session_id(web_app):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    result = fn(
        1, None, "current", "", [], None, [], None, None, None,
        None, None,
    )
    assert result[1] == "No session"
    assert result[2] == "danger"
    assert result[8] is True


def test_start_recompute_confirm_gate(web_app):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    store_data = {"recomputed": True}
    result = fn(
        1, None, "current", "", [], None, [], None, None, None,
        "sess-2", store_data,
    )
    assert result[0] == {"display": "none"}
    assert result[10] == {"display": "block"}


def test_start_recompute_yes_bypasses_confirm_gate(web_app, monkeypatch):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE_YES}.n_clicks")
    monkeypatch.setattr(recompute_mod, "cache_get", lambda sid: None)
    result = fn(
        None, 1, "current", "", [], None, [], None, None, None,
        "sess-3", {"recomputed": True},
    )
    assert result[1] == "No data"


def test_start_recompute_current_source_no_cache_no_store(web_app, monkeypatch):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    monkeypatch.setattr(recompute_mod, "cache_get", lambda sid: None)
    result = fn(
        1, None, "current", "", [], None, [], None, None, None,
        "sess-4", None,
    )
    assert "Load Data first" in result[9]


def test_start_recompute_current_source_cache_expired(web_app, monkeypatch):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    monkeypatch.setattr(recompute_mod, "cache_get", lambda sid: None)
    result = fn(
        1, None, "current", "", [], None, [], None, None, None,
        "sess-5", {"n_stations": 3},
    )
    assert "cache expired" in result[9]


def test_start_recompute_path_source_no_path(web_app):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    result = fn(
        1, None, "path", "   ", [], None, [], None, None, None,
        "sess-6", None,
    )
    assert result[1] == "No path"


def test_start_recompute_path_source_not_found(web_app):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    result = fn(
        1, None, "path", "/no/such/dir/xyz", [], None, [], None, None,
        None, "sess-7", None,
    )
    assert result[1] == "Not found"


def test_start_recompute_current_source_starts_job(
    web_app, monkeypatch, sync_thread, real_sites
):
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")
    monkeypatch.setattr(recompute_mod, "cache_get", lambda sid: real_sites)

    session_id = "sess-8"
    result = fn(
        1, None, "current", "", ["on"], "30", ["Z"], "1", "100", None,
        session_id, None,
    )
    assert result[0] == {"display": "block"}
    assert result[8] is False  # interval enabled
    job = recompute_mod._JOBS[session_id]
    assert job["status"] == "done"


def test_start_recompute_path_source_starts_job(web_app, sync_thread):
    if not _HAS_EDIS:
        pytest.skip("3edis dataset not found")
    fn = _start_recompute_fn(web_app)
    _set_triggered(f"{IDs.BTN_RECOMPUTE}.n_clicks")

    session_id = "sess-9"
    result = fn(
        1, None, "path", str(_EDI_DIR), [], None, [], None, None, None,
        session_id, None,
    )
    assert result[0] == {"display": "block"}
    job = recompute_mod._JOBS[session_id]
    assert job["status"] == "done"


# ─────────────────────────────────────────────────────────────────────────────
# poll_recompute
# ─────────────────────────────────────────────────────────────────────────────


def _poll_fn(web_app):
    return _cb_by_input(
        web_app,
        f"{IDs.RECOMPUTE_TICKER}.children",
        IDs.RECOMPUTE_INTERVAL,
    )


def test_poll_no_session_id(web_app):
    fn = _poll_fn(web_app)
    assert fn(1, None, None) == (no_update,) * 11


def test_poll_missing_job(web_app):
    fn = _poll_fn(web_app)
    assert fn(1, "sess-missing", None) == (no_update,) * 11


def test_poll_loading_phase_scanning(web_app):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-loading"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "starting",
            "phase": "loading",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 0,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }
    ticker, pct, animated, striped, log, badge, color, disabled, *_ = fn(
        1, session_id, None
    )
    assert "Scanning folder" in ticker
    assert pct == 0
    assert disabled is False


def test_poll_loading_phase_found_some(web_app):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-found"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "running",
            "phase": "loading",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 5,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }
    ticker = fn(1, session_id, None)[0]
    assert "5 EDI files found" in ticker


def test_poll_recomputing_with_progress_and_failures(web_app):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-running"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "running",
            "phase": "recomputing",
            "total": 10,
            "done": 3,
            "failed": 1,
            "found": 10,
            "current": "HBH07",
            "log": ["✓ a", "✗ b: err"],
            "sites": None,
            "err_msg": "",
            "t_start": time.time() - 5,
        }
    ticker, pct = fn(1, session_id, None)[:2]
    assert "HBH07" in ticker
    assert "3 / 10" in ticker
    assert "1 ✗" in ticker
    assert pct == 30


def test_poll_recomputing_zero_done_no_eta(web_app):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-zero"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "starting",
            "phase": "recomputing",
            "total": 5,
            "done": 0,
            "failed": 0,
            "found": 5,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "",
            "t_start": time.time(),
        }
    ticker, pct, _, _, _, badge_txt = fn(1, session_id, None)[:6]
    assert "ETA" not in ticker
    assert badge_txt == "Starting"


def test_poll_done_updates_store_and_cache(web_app, monkeypatch, real_sites):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-done"
    cached = {}
    monkeypatch.setattr(
        recompute_mod, "cache_set", lambda sid, sites: cached.__setitem__(sid, sites)
    )

    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "done",
            "phase": "recomputing",
            "total": 3,
            "done": 3,
            "failed": 0,
            "found": 3,
            "current": "",
            "log": ["✓ a", "✓ b", "✓ c"],
            "sites": real_sites,
            "err_msg": "",
            "t_start": time.time() - 2,
        }

    result = fn(1, session_id, None)
    ticker, pct, animated, striped, log_children, badge_txt, color, disabled, \
        new_store, feedback, icon = result

    assert pct == 100
    assert disabled is True
    assert color == "success"
    assert new_store["recomputed"] is True
    assert session_id in cached
    assert "recomputed" in feedback


def test_poll_done_with_failures_message(web_app, monkeypatch, real_sites):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-done-fail"
    monkeypatch.setattr(recompute_mod, "cache_set", lambda sid, sites: None)

    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "done",
            "phase": "recomputing",
            "total": 4,
            "done": 4,
            "failed": 2,
            "found": 4,
            "current": "",
            "log": [],
            "sites": real_sites,
            "err_msg": "",
            "t_start": time.time() - 1,
        }

    result = fn(1, session_id, None)
    badge_txt = result[5]
    feedback = result[9]
    assert "failed" in badge_txt
    assert "failed" in feedback


def test_poll_done_sites_to_store_exception_falls_back(
    web_app, monkeypatch, real_sites
):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-done-exc"
    monkeypatch.setattr(recompute_mod, "cache_set", lambda sid, sites: None)

    def _boom(sites, old_store=None):
        raise RuntimeError("store boom")

    monkeypatch.setattr(recompute_mod, "_sites_to_store", _boom)

    old_store = {"n_stations": 42}
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "done",
            "phase": "recomputing",
            "total": 1,
            "done": 1,
            "failed": 0,
            "found": 1,
            "current": "",
            "log": [],
            "sites": real_sites,
            "err_msg": "",
            "t_start": time.time(),
        }

    result = fn(1, session_id, old_store)
    new_store = result[8]
    assert new_store == old_store


def test_poll_error_status(web_app):
    fn = _poll_fn(web_app)
    session_id = "sess-poll-error"
    with recompute_mod._LOCK:
        recompute_mod._JOBS[session_id] = {
            "status": "error",
            "phase": "recomputing",
            "total": 0,
            "done": 0,
            "failed": 0,
            "found": 0,
            "current": "",
            "log": [],
            "sites": None,
            "err_msg": "disk exploded",
            "t_start": time.time(),
        }
    result = fn(1, session_id, None)
    assert "disk exploded" in result[0]
    assert result[5] == "✗ Error"
    assert result[6] == "danger"
    assert result[7] is True

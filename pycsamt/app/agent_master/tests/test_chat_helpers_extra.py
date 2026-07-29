# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for smaller/pure helper functions in
pycsamt.app.agent_master.callbacks.chat that were not already covered by
test_chat_ui.py / test_data_overview.py / test_metric_line_resolution.py."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


# ── figure helpers ────────────────────────────────────────────────────────────


def _tiny_fig():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1, 1))
    ax.plot([0, 1], [0, 1])
    return fig


class TestFigToB64:
    def test_returns_base64_png(self):
        import base64

        from pycsamt.app.agent_master.callbacks.chat import _fig_to_b64

        b64 = _fig_to_b64(_tiny_fig())
        raw = base64.b64decode(b64)
        assert raw[:8] == b"\x89PNG\r\n\x1a\n"


class TestFigThumbItem:
    def test_truncates_long_title(self):
        from pycsamt.app.agent_master.callbacks.chat import _fig_thumb_item

        long_title = "x" * 40
        el = _fig_thumb_item("k1", long_title, "b64data")
        text = str(el)
        assert "..." in text

    def test_short_title_not_truncated(self):
        from pycsamt.app.agent_master.callbacks.chat import _fig_thumb_item

        el = _fig_thumb_item("k1", "short", "b64data")
        assert "..." not in str(el)


class TestFigAccordion:
    def test_shows_figure_count(self):
        from pycsamt.app.agent_master.callbacks.chat import _fig_accordion

        figs = {
            "f1": {"title": "Fig 1", "b64": "aaa"},
            "f2": {"title": "Fig 2", "b64": "bbb"},
        }
        text = str(_fig_accordion(figs))
        assert "Figures (2)" in text


# ── message bubble builders ───────────────────────────────────────────────────


class TestMid:
    def test_generates_prefixed_stable_id(self):
        from pycsamt.app.agent_master.callbacks.chat import _mid

        m = _mid()
        assert m.startswith("am-msg-")
        assert len(m) == len("am-msg-") + 8


class TestUserBubbleWithMid:
    def test_includes_pin_button_and_id_when_mid_given(self):
        from pycsamt.app.agent_master.callbacks.chat import _user_bubble

        div = _user_bubble("hi", mid="am-msg-abc123")
        assert div.id == "am-msg-abc123"
        assert "am-pin-btn" in str(div)

    def test_no_id_or_pin_button_without_mid(self):
        from pycsamt.app.agent_master.callbacks.chat import _user_bubble

        div = _user_bubble("hi")
        assert not hasattr(div, "id") or div.id is None
        assert "am-pin-btn" not in str(div)


class TestExecStepRow:
    def test_done_status(self):
        from pycsamt.app.agent_master.callbacks.chat import _exec_step_row

        row = _exec_step_row("Load data", "done")
        assert "bi-check-lg" in str(row)

    def test_error_status(self):
        from pycsamt.app.agent_master.callbacks.chat import _exec_step_row

        row = _exec_step_row("Load data", "error")
        assert "bi-exclamation-lg" in str(row)

    def test_running_status(self):
        from pycsamt.app.agent_master.callbacks.chat import _exec_step_row

        row = _exec_step_row("Load data", "running")
        assert "am-tl-spin" in str(row)

    def test_pending_status(self):
        from pycsamt.app.agent_master.callbacks.chat import _exec_step_row

        row = _exec_step_row("Load data", "pending")
        assert "am-tl-step pending" in str(row)


class TestKindHeader:
    def test_unknown_kind_returns_none(self):
        from pycsamt.app.agent_master.callbacks.chat import _kind_header

        assert _kind_header("not-a-real-kind") is None

    def test_known_kind_builds_chip(self, monkeypatch):
        # _HEADER_KINDS is currently an empty frozenset in production, so
        # _kind_header() always returns None for every real kind — this
        # exercises the chip-building branch directly to keep it covered.
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        kind = next(iter(chat_mod._KIND_HEADER))
        monkeypatch.setattr(chat_mod, "_HEADER_KINDS", frozenset({kind}))
        header = chat_mod._kind_header(kind)
        assert header is not None
        assert f"am-kind-{kind}" in header.className


class TestAgentBubbleWithMidAndSteps:
    def test_bubble_with_steps_figs_code_and_mid(self):
        from pycsamt.app.agent_master.callbacks.chat import _agent_bubble

        div = _agent_bubble(
            "Done.",
            steps=[{"label": "Load", "status": "done"}],
            figs={"f1": {"title": "Fig", "b64": "aaa"}},
            code="print(1)",
            mid="am-msg-xyz",
        )
        text = str(div)
        assert div.id == "am-msg-xyz"
        assert "am-pin-btn" in text
        assert "1 figure" in text
        assert "am-trace" in text

    def test_bubble_without_mid_has_empty_toolbar(self):
        from pycsamt.app.agent_master.callbacks.chat import _agent_bubble

        div = _agent_bubble("Done.")
        assert not hasattr(div, "id") or div.id is None


# ── waiting bubbles ────────────────────────────────────────────────────────────


class TestLineWaitingBubble:
    def test_has_expected_id_and_prompt(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _line_waiting_bubble,
        )

        div = _line_waiting_bubble()
        assert div.id == "am-line-waiting-bubble"
        assert "Multiple survey" in str(div)


# ── line reference extraction ─────────────────────────────────────────────────


class TestExtractLineRefOrdinals:
    def test_ordinal_word_resolves_to_number(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _extract_line_ref,
        )

        assert _extract_line_ref("run qc on the first line", {}) == "1"
        assert _extract_line_ref("plot the second profile", {}) == "2"

    def test_no_match_returns_none(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _extract_line_ref,
        )

        assert _extract_line_ref("run qc", {}) is None

    def test_line_n_pattern_strips_punctuation(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _extract_line_ref,
        )

        assert _extract_line_ref("run qc on line L22.", {}) == "l22"


class TestMatchGroup:
    def test_exact_key_match(self):
        from pycsamt.app.agent_master.callbacks.chat import _match_group

        assert _match_group("L1", {"L1": [], "L2": []}) == "L1"

    def test_case_insensitive_match(self):
        from pycsamt.app.agent_master.callbacks.chat import _match_group

        assert _match_group("l1", {"L1": []}) == "L1"

    def test_numeric_ref_unique_embedded_number(self):
        from pycsamt.app.agent_master.callbacks.chat import _match_group

        assert _match_group("22", {"L22PLT": [], "L18PLT": []}) == "L22PLT"

    def test_numeric_ref_ambiguous_returns_none(self):
        from pycsamt.app.agent_master.callbacks.chat import _match_group

        assert _match_group("1", {"L1A": [], "L1B": []}) is None

    def test_named_substring_match(self):
        from pycsamt.app.agent_master.callbacks.chat import _match_group

        assert _match_group("plt", {"L22PLT": []}) == "L22PLT"


class TestIsPinnOrHybrid:
    def test_matches_pinn_keyword(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _is_pinn_or_hybrid,
        )

        assert _is_pinn_or_hybrid("use a PINN model")

    def test_matches_hybrid_keyword(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _is_pinn_or_hybrid,
        )

        assert _is_pinn_or_hybrid("warm start the run")

    def test_no_match(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _is_pinn_or_hybrid,
        )

        assert not _is_pinn_or_hybrid("run standard qc")


# ── application launchers ─────────────────────────────────────────────────────


class TestDetectAppRequest:
    def test_desktop_request(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _detect_app_request,
        )

        assert _detect_app_request("open the desktop app") == ("desktop", "")

    def test_mapview_request(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _detect_app_request,
        )

        assert _detect_app_request("open the map") == ("mapview", "")

    def test_mapview_viz_redirect(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _VIZ_REDIRECT_REASON,
            _detect_app_request,
        )

        assert _detect_app_request("show a 3d map") == (
            "mapview",
            _VIZ_REDIRECT_REASON,
        )

    def test_web_app_request(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _detect_app_request,
        )

        assert _detect_app_request("open web app") == ("web", "")

    def test_complex_viz_redirect(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _VIZ_REDIRECT_REASON,
            _detect_app_request,
        )

        assert _detect_app_request("browse edis") == (
            "web",
            _VIZ_REDIRECT_REASON,
        )

    def test_no_match_returns_none(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _detect_app_request,
        )

        assert _detect_app_request("run qc") is None


class TestPortInUse:
    def test_free_port_not_in_use(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _free_port,
            _port_in_use,
        )

        port = _free_port(0)
        assert _port_in_use(port) is False

    def test_bound_port_is_in_use(self):
        import socket

        from pycsamt.app.agent_master.callbacks.chat import _port_in_use

        srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        srv.bind(("127.0.0.1", 0))
        srv.listen(1)
        port = srv.getsockname()[1]
        try:
            assert _port_in_use(port) is True
        finally:
            srv.close()


class TestFreePort:
    def test_returns_preferred_when_available(self):
        import socket

        from pycsamt.app.agent_master.callbacks.chat import _free_port

        # bind briefly to find a genuinely free ephemeral port to prefer
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("127.0.0.1", 0))
        preferred = s.getsockname()[1]
        s.close()
        assert _free_port(preferred) == preferred

    def test_falls_back_on_bind_error(self, monkeypatch):
        import socket

        from pycsamt.app.agent_master.callbacks.chat import _free_port

        orig_socket_cls = socket.socket
        calls = {"n": 0}

        class _FlakySocket(orig_socket_cls):
            def bind(self, addr):
                calls["n"] += 1
                if calls["n"] == 1:
                    raise OSError("address in use")
                return super().bind(addr)

        monkeypatch.setattr(socket, "socket", _FlakySocket)
        got = _free_port(12345)
        assert isinstance(got, int) and got > 0


class TestEnsureMapview:
    def test_reuses_running_process(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        class _FakeProc:
            def poll(self):
                return None

        chat_mod._EXT_APPS["mapview"] = {
            "proc": _FakeProc(),
            "url": "http://127.0.0.1:9999",
        }
        try:
            assert chat_mod._ensure_mapview() == "http://127.0.0.1:9999"
        finally:
            chat_mod._EXT_APPS["mapview"] = {"proc": None, "url": None}

    def test_points_at_already_listening_port(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        chat_mod._EXT_APPS["mapview"] = {"proc": None, "url": None}
        monkeypatch.setattr(chat_mod, "_port_in_use", lambda port: True)
        try:
            url = chat_mod._ensure_mapview()
            assert url == f"http://127.0.0.1:{chat_mod._MAPVIEW_PORT}"
        finally:
            chat_mod._EXT_APPS["mapview"] = {"proc": None, "url": None}

    def test_spawns_new_process(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        chat_mod._EXT_APPS["mapview"] = {"proc": None, "url": None}
        monkeypatch.setattr(chat_mod, "_port_in_use", lambda port: False)
        monkeypatch.setattr(chat_mod, "_free_port", lambda preferred: 9123)

        class _FakePopen:
            def __init__(self, *a, **k):
                pass

        import subprocess as sp

        monkeypatch.setattr(sp, "Popen", _FakePopen)
        try:
            url = chat_mod._ensure_mapview()
            assert url == "http://127.0.0.1:9123"
        finally:
            chat_mod._EXT_APPS["mapview"] = {"proc": None, "url": None}


class TestEnsureDesktop:
    def test_already_running(self):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        class _FakeProc:
            def poll(self):
                return None

        chat_mod._EXT_APPS["desktop"] = {"proc": _FakeProc()}
        try:
            ok, note = chat_mod._ensure_desktop()
            assert ok is True
            assert "already running" in note
        finally:
            chat_mod._EXT_APPS["desktop"] = {"proc": None}

    def test_popen_failure_reports_error(self, monkeypatch):
        import subprocess as sp

        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        chat_mod._EXT_APPS["desktop"] = {"proc": None}

        def boom(*a, **k):
            raise OSError("no exe")

        monkeypatch.setattr(sp, "Popen", boom)
        try:
            ok, note = chat_mod._ensure_desktop()
            assert ok is False
            assert "could not be started" in note
        finally:
            chat_mod._EXT_APPS["desktop"] = {"proc": None}

    def test_process_exits_immediately(self, monkeypatch):
        import subprocess as sp

        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        chat_mod._EXT_APPS["desktop"] = {"proc": None}

        class _FakeProc:
            def poll(self):
                return 1  # already exited

        monkeypatch.setattr(sp, "Popen", lambda *a, **k: _FakeProc())
        monkeypatch.setattr(chat_mod.time, "sleep", lambda *_a: None)
        try:
            ok, note = chat_mod._ensure_desktop()
            assert ok is False
            assert "exited immediately" in note
        finally:
            chat_mod._EXT_APPS["desktop"] = {"proc": None}

    def test_process_stays_alive(self, monkeypatch):
        import subprocess as sp

        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        chat_mod._EXT_APPS["desktop"] = {"proc": None}

        class _FakeProc:
            def poll(self):
                return None  # still running

        monkeypatch.setattr(sp, "Popen", lambda *a, **k: _FakeProc())
        monkeypatch.setattr(chat_mod.time, "sleep", lambda *_a: None)
        try:
            ok, note = chat_mod._ensure_desktop()
            assert ok is True
            assert "should appear shortly" in note
        finally:
            chat_mod._EXT_APPS["desktop"] = {"proc": None}


class TestEnsureWebApp:
    def test_reuses_running_server(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        monkeypatch.setattr(
            chat_mod,
            "_webapp_state",
            {"running": True, "url": "http://127.0.0.1:1234"},
        )
        assert chat_mod._ensure_web_app() == "http://127.0.0.1:1234"

    def test_starts_new_server_thread(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        monkeypatch.setattr(chat_mod, "_webapp_state", {"running": False, "url": None})
        monkeypatch.setattr(chat_mod, "_free_port", lambda preferred: 8199)

        started = {}

        class _FakeThread:
            def __init__(self, target, daemon):
                started["target"] = target

            def start(self):
                started["started"] = True

        monkeypatch.setattr(chat_mod.threading, "Thread", _FakeThread)
        url = chat_mod._ensure_web_app()
        assert url == "http://127.0.0.1:8199"
        assert started.get("started") is True
        assert chat_mod._webapp_state["running"] is True


class TestLaunchBubble:
    def test_success_with_url(self):
        from pycsamt.app.agent_master.callbacks.chat import _launch_bubble

        div = _launch_bubble("web", url="http://127.0.0.1:8051")
        text = str(div)
        assert "Launching pyCSAMT Web App" in text
        assert "http://127.0.0.1:8051" in text
        assert "Server is starting" in text

    def test_failure_shows_could_not_launch(self):
        from pycsamt.app.agent_master.callbacks.chat import _launch_bubble

        div = _launch_bubble("desktop", ok=False, note="boom")
        text = str(div)
        assert "Could not launch" in text
        assert "boom" in text

    def test_custom_reason_overrides_default_desc(self):
        from pycsamt.app.agent_master.callbacks.chat import _launch_bubble

        div = _launch_bubble("mapview", reason="Custom reason here")
        assert "Custom reason here" in str(div)


# ── capability text helpers ───────────────────────────────────────────────────


class TestApiKeyHint:
    def test_mentions_settings_and_providers(self):
        from pycsamt.app.agent_master.callbacks.chat import _api_key_hint

        text = _api_key_hint()
        assert "Settings" in text
        assert "Anthropic" in text


class TestCorrectionCapabilityBlock:
    def test_returns_bullet_list_grouped_by_category(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _correction_capability_block,
        )

        text = _correction_capability_block()
        assert isinstance(text, str)


# ── prep-workflow helpers ─────────────────────────────────────────────────────


class TestPrepAllLinesRequested:
    def test_matches_all_lines_phrase(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _prep_all_lines_requested,
        )

        assert _prep_all_lines_requested("prep all survey lines")
        assert _prep_all_lines_requested("run it on every profile")

    def test_no_match(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _prep_all_lines_requested,
        )

        assert not _prep_all_lines_requested("prep line 1")


class TestSafeDirname:
    def test_sanitizes_special_chars(self):
        from pycsamt.app.agent_master.callbacks.chat import _safe_dirname

        assert _safe_dirname("Line/Name*1") == "Line_Name_1"

    def test_spaces_become_underscores(self):
        from pycsamt.app.agent_master.callbacks.chat import _safe_dirname

        assert _safe_dirname("My Line") == "My_Line"

    def test_empty_falls_back_to_line(self):
        from pycsamt.app.agent_master.callbacks.chat import _safe_dirname

        assert _safe_dirname("   ") == "line"
        assert _safe_dirname("") == "line"


class TestFmtFileSize:
    def test_kb_size(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_file_size

        f = tmp_path / "small.bin"
        f.write_bytes(b"x" * 2048)
        assert _fmt_file_size(f).endswith("KB")

    def test_mb_size(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_file_size

        f = tmp_path / "big.bin"
        f.write_bytes(b"x" * (2 * 1024 * 1024))
        assert _fmt_file_size(f).endswith("MB")

    def test_missing_file_returns_empty(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_file_size

        assert _fmt_file_size(tmp_path / "ghost.bin") == ""


class TestPrepRagNote:
    def test_degrades_silently_when_rag_unavailable(self, monkeypatch):
        import sys

        from pycsamt.app.agent_master.callbacks.chat import _prep_rag_note

        monkeypatch.setitem(
            sys.modules,
            "pycsamt.assistant.rag.context_builder",
            None,
        )
        assert _prep_rag_note("some query") == ""


# ── frequency formatting ──────────────────────────────────────────────────────


class TestFmtHz:
    def test_below_1khz(self):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_hz

        assert _fmt_hz(0.125) == "0.125 Hz"

    def test_above_1khz(self):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_hz

        assert _fmt_hz(8190) == "8.19 kHz"


class TestFmtFreq:
    def test_empty_returns_na(self):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_freq

        assert _fmt_freq(None) == "n/a"
        assert _fmt_freq([]) == "n/a"

    def test_range_formats_both_ends(self):
        from pycsamt.app.agent_master.callbacks.chat import _fmt_freq

        assert _fmt_freq([0.1, 1000.0]) == "0.1 Hz – 1 kHz"


# ── session / recent-runs helpers ─────────────────────────────────────────────


class TestSessionHasData:
    def test_false_when_session_unavailable(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        monkeypatch.setattr(chat_mod, "_session", lambda: None)
        assert chat_mod._session_has_data() is False

    def test_true_when_edi_path_set(self, monkeypatch):
        import types

        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        fake = types.SimpleNamespace(edi_path="/some/path")
        monkeypatch.setattr(chat_mod, "_session", lambda: fake)
        assert chat_mod._session_has_data() is True


class TestResetSessionNoop:
    def test_noop_when_session_unavailable(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        monkeypatch.setattr(chat_mod, "_session", lambda: None)
        chat_mod._reset_session()  # should not raise


class TestNamesRegistryLine:
    def test_false_when_registry_unavailable(self, monkeypatch):
        import sys

        from pycsamt.app.agent_master.callbacks.chat import (
            _names_registry_line,
        )

        monkeypatch.setitem(
            sys.modules, "pycsamt.assistant.tools.project_registry", None
        )
        assert _names_registry_line("run on line L22") is False

    def test_false_when_registry_returns_none(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        class _FakeRegistry:
            @staticmethod
            def from_default():
                return None

        import pycsamt.assistant.tools.project_registry as pr_mod

        monkeypatch.setattr(pr_mod, "ProjectRegistry", _FakeRegistry)
        assert chat_mod._names_registry_line("run on line L22") is False

    def test_true_when_line_resolves(self, monkeypatch):
        from pycsamt.app.agent_master.callbacks import chat as chat_mod

        class _FakeRegistry:
            def find_line_in_text(self, text):
                return "L22"

            def resolve_line(self, line):
                return {"exists": True}

            @staticmethod
            def from_default():
                return _FakeRegistry()

        import pycsamt.assistant.tools.project_registry as pr_mod

        monkeypatch.setattr(pr_mod, "ProjectRegistry", _FakeRegistry)
        assert chat_mod._names_registry_line("run on line L22") is True


class TestRecentRunsUnavailable:
    def test_returns_empty_list_on_error(self, monkeypatch):
        import sys

        from pycsamt.app.agent_master.callbacks.chat import _recent_runs

        monkeypatch.setitem(
            sys.modules,
            "pycsamt.assistant.memory.workflow_history",
            None,
        )
        assert _recent_runs() == []


class TestRunItem:
    def test_ok_status_uses_check_icon(self):
        from pycsamt.app.agent_master.callbacks.chat import _run_item

        div = _run_item(
            {
                "workflow": "qc",
                "status": "success",
                "timestamp": "2026-01-01T10:30:00",
                "summary": "ok",
            }
        )
        assert "bi-check-circle-fill" in str(div)

    def test_failed_status_uses_warning_icon(self):
        from pycsamt.app.agent_master.callbacks.chat import _run_item

        div = _run_item({"workflow": "qc", "status": "failed"})
        assert "bi-exclamation-triangle-fill" in str(div)


def test_record_run_swallows_import_errors(monkeypatch):
    import sys

    from pycsamt.app.agent_master.callbacks.chat import _record_run

    monkeypatch.setitem(sys.modules, "pycsamt.assistant.memory.workflow_history", None)
    _record_run(
        workflow="qc",
        path="/x",
        output_dir="/out",
        status="success",
        summary="ok",
        n_figures=0,
    )  # should not raise

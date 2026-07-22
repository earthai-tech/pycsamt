# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AgentWorker (Phase 4).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — used for the real processing-mode success
path (a fast emtools plot function), so ``_run_processing``'s ``_call``
helper and its ``result_plot`` follow-up branch are exercised for real.

LLM mode is exercised against a fake agent class monkeypatched onto
``pycsamt.agents`` (constructing/calling a real LLM agent would need a
network round-trip and an API key), matching the ``AGENT_REGISTRY``
metadata for a real entry (``MT Loader`` / ``MTLoaderAgent``) so the
worker's own kwarg-filtering logic runs unmodified.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import logging

from pycsamt.app.desktop.workers.agent_worker import (
    AgentWorker,
    ChatWorker,
    _has_figure,
    _is_sites_like,
    _SignalLogHandler,
    _StreamCapture,
)

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


class _FakeLLMAgent:
    last_instance = None

    def __init__(self, recursive=None, on_dup=None, api_key=None):
        self.recursive = recursive
        self.on_dup = on_dup
        self.api_key = api_key
        _FakeLLMAgent.last_instance = self

    def execute(self, input_data):
        return SimpleNamespace(
            status="success",
            elapsed_seconds=0.1,
            summary="Loaded OK",
            warnings=["one warning"],
            llm_interpretation="interpretation text",
        )

    def query_llm(self, prompt):
        return f"reply to: {prompt}"

# ── _SignalLogHandler ─────────────────────────────────────────────────────


def test_signal_log_handler_emits_line():
    received = []
    handler = _SignalLogHandler(received.append)
    handler.setFormatter(logging.Formatter("%(message)s"))
    record = logging.LogRecord("test", logging.INFO, "", 0, "hello", (), None)
    handler.emit(record)
    assert any("hello" in line for line in received)


def test_signal_log_handler_ignores_exceptions():
    def bad_fn(line):
        raise RuntimeError("sink error")

    handler = _SignalLogHandler(bad_fn)
    record = logging.LogRecord("test", logging.INFO, "", 0, "msg", (), None)
    handler.emit(record)  # must not propagate


# ── AgentWorker construction ──────────────────────────────────────────────


def test_agent_worker_creates(qapp):
    w = AgentWorker(
        agent_name="QC Quicklook",
        sites=None,
        params={},
    )
    assert w is not None


def test_agent_worker_unknown_name_emits_error(qapp):
    errors = []
    w = AgentWorker(agent_name="UNKNOWN_AGENT", sites=None, params={})
    w.error.connect(errors.append)
    w.run()  # run synchronously in test
    assert len(errors) == 1
    assert "UNKNOWN_AGENT" in errors[0]


def test_agent_worker_cancel_flag(qapp):
    w = AgentWorker(agent_name="QC Quicklook", sites=None, params={})
    assert not w._cancelled
    w.cancel()
    assert w._cancelled


# ── Processing mode — no-data graceful failure ────────────────────────────


def test_agent_worker_processing_no_sites_emits_error(qapp):
    errors = []
    w = AgentWorker(
        agent_name="QC Quicklook",
        sites=None,
        params={},
    )
    w.error.connect(errors.append)
    w.run()
    # Should emit an error because sites=None
    assert len(errors) == 1


# ── Signals exist ─────────────────────────────────────────────────────────


def test_agent_worker_has_required_signals(qapp):
    w = AgentWorker("QC Quicklook", sites=None, params={})
    assert hasattr(w, "log_line")
    assert hasattr(w, "progress")
    assert hasattr(w, "result")
    assert hasattr(w, "error")


# ── _StreamCapture ─────────────────────────────────────────────────────────


def test_stream_capture_forwards_nonempty_lines():
    lines = []
    cap = _StreamCapture(lines.append)
    n = cap.write("first\n\nsecond  \n")
    assert lines == ["first", "second"]
    assert n == len("first\n\nsecond  \n")


def test_stream_capture_flush_is_noop():
    cap = _StreamCapture(lambda line: None)
    cap.flush()  # must not raise


# ── _has_figure / _is_sites_like ────────────────────────────────────────────


class TestHasFigureAndSitesLike:
    def test_has_figure_none(self):
        assert _has_figure(None) is False

    def test_has_figure_savefig_object(self):
        obj = SimpleNamespace(savefig=lambda *a, **k: None)
        assert _has_figure(obj) is True

    def test_has_figure_get_figure_object(self):
        obj = SimpleNamespace(get_figure=lambda: None)
        assert _has_figure(obj) is True

    def test_has_figure_nested_data_dict(self):
        fig_like = SimpleNamespace(savefig=lambda *a, **k: None)
        obj = SimpleNamespace(data={"k": fig_like})
        assert _has_figure(obj) is True

    def test_has_figure_plain_object_false(self):
        assert _has_figure(SimpleNamespace(data={"k": 1})) is False

    def test_is_sites_like_none(self):
        assert _is_sites_like(None) is False

    def test_is_sites_like_figure_is_false(self):
        obj = SimpleNamespace(savefig=lambda *a, **k: None)
        assert _is_sites_like(obj) is False

    def test_is_sites_like_structural_fallback(self):
        obj = SimpleNamespace(edic={}, by_index=lambda i: None, as_list=lambda: [])
        assert _is_sites_like(obj) is True

    def test_is_sites_like_plain_object_false(self):
        assert _is_sites_like(SimpleNamespace()) is False

    def test_is_sites_like_real_sites(self, willy_sites):
        assert _is_sites_like(willy_sites) is True


# ── Processing mode — real success path ─────────────────────────────────────


class TestProcessingModeReal:
    def test_qc_quicklook_success(self, willy_sites):
        results = []
        w = AgentWorker(
            agent_name="QC Quicklook", sites=willy_sites, params={}
        )
        w.result.connect(results.append)
        w.run()
        assert len(results) == 1

    def test_dimensionality_with_separate_result_plot(self, willy_sites):
        """
        'Dimensionality' has fn_name='classify_dimensionality' and a
        *different* result_plot='plot_dimensionality_psection' -- exercises
        the extra plotting branch in _run_processing (the classification
        table returned by the main fn isn't itself renderable, so the
        worker calls the separate plot function and swaps in its Axes).
        """
        results = []
        w = AgentWorker(
            agent_name="Dimensionality",
            sites=willy_sites,
            params={"skew_th": 3.0, "ellipt_th": 0.2},
        )
        w.result.connect(results.append)
        w.run()
        assert len(results) == 1

    def test_unknown_function_reports_error(self, willy_sites, monkeypatch):
        import pycsamt.app.desktop.workers.agent_worker as aw_mod

        monkeypatch.setattr(
            aw_mod,
            "get_entry",
            lambda name: {
                "type": "processing",
                "fn_name": "definitely_not_a_real_fn",
                "params": {},
            },
        )
        errors = []
        w = AgentWorker(agent_name="Bogus", sites=willy_sites, params={})
        w.error.connect(errors.append)
        w.run()
        assert len(errors) == 1

    def test_result_plot_failure_logged_and_swallowed(
        self, willy_sites, monkeypatch
    ):
        import pycsamt.emtools as et

        def _boom(*a, **k):
            raise RuntimeError("plot boom")

        monkeypatch.setattr(et, "plot_dimensionality_psection", _boom)
        results = []
        lines = []
        w = AgentWorker(
            agent_name="Dimensionality",
            sites=willy_sites,
            params={"skew_th": 3.0, "ellipt_th": 0.2},
        )
        w.result.connect(results.append)
        w.log_line.connect(lines.append)
        w.run()
        assert len(results) == 1
        assert any("Plot step failed" in line for line in lines)

    def test_cancelled_before_completion_suppresses_result(
        self, willy_sites
    ):
        w = AgentWorker(agent_name="QC Quicklook", sites=willy_sites, params={})
        w.cancel()
        results = []
        w.result.connect(results.append)
        w.run()
        assert results == []


# ── LLM mode ─────────────────────────────────────────────────────────────────


class TestLLMMode:
    def test_llm_mode_success(self, monkeypatch):
        import pycsamt.agents as ag

        monkeypatch.setattr(ag, "MTLoaderAgent", _FakeLLMAgent, raising=False)
        results = []
        ready = []
        w = AgentWorker(
            agent_name="MT Loader",
            sites=object(),
            params={"recursive": True, "on_dup": "replace"},
            api_key="sk-test",
        )
        w.result.connect(results.append)
        w.agent_ready.connect(ready.append)
        w.run()
        assert len(results) == 1
        assert results[0].status == "success"
        assert len(ready) == 1
        assert ready[0].api_key == "sk-test"
        assert ready[0].recursive is True

    def test_llm_mode_forwards_extra_params_as_input_data(self, monkeypatch):
        import pycsamt.agents as ag

        captured = {}

        class _CapturingAgent(_FakeLLMAgent):
            def execute(self, input_data):
                captured.update(input_data)
                return super().execute(input_data)

        monkeypatch.setattr(
            ag, "MTLoaderAgent", _CapturingAgent, raising=False
        )
        w = AgentWorker(
            agent_name="MT Loader",
            sites="the-sites",
            params={
                "recursive": True,
                "on_dup": "replace",
                "user_prompt": "extra question",
            },
        )
        w.run()
        assert captured["sites"] == "the-sites"
        assert captured["user_prompt"] == "extra question"

    def test_llm_mode_warnings_logged(self, monkeypatch):
        import pycsamt.agents as ag

        monkeypatch.setattr(ag, "MTLoaderAgent", _FakeLLMAgent, raising=False)
        lines = []
        w = AgentWorker(agent_name="MT Loader", sites=object(), params={})
        w.log_line.connect(lines.append)
        w.run()
        assert any("one warning" in line for line in lines)

    def test_llm_mode_exception_reports_error(self, monkeypatch):
        import pycsamt.agents as ag

        class _BoomAgent(_FakeLLMAgent):
            def execute(self, input_data):
                raise RuntimeError("agent execution failed")

        monkeypatch.setattr(ag, "MTLoaderAgent", _BoomAgent, raising=False)
        errors = []
        w = AgentWorker(agent_name="MT Loader", sites=object(), params={})
        w.error.connect(errors.append)
        w.run()
        assert errors == ["agent execution failed"]


# ── ChatWorker ────────────────────────────────────────────────────────────────


class TestChatWorker:
    def test_success_without_context(self, qapp):
        agent = SimpleNamespace(
            api_key="sk-x", query_llm=lambda prompt: f"answer: {prompt}"
        )
        replies = []
        w = ChatWorker(agent, "what is this?")
        w.reply_done.connect(replies.append)
        w.run()
        assert replies == ["answer: what is this?"]

    def test_success_with_context_prefixes_prompt(self, qapp):
        captured = {}

        def _query(prompt):
            captured["prompt"] = prompt
            return "ok"

        agent = SimpleNamespace(api_key="sk-x", query_llm=_query)
        w = ChatWorker(agent, "follow up?", context="earlier summary")
        w.run()
        assert "earlier summary" in captured["prompt"]
        assert "follow up?" in captured["prompt"]

    def test_no_api_key_emits_error(self, qapp):
        agent = SimpleNamespace(api_key=None, query_llm=lambda p: "unused")
        errors = []
        w = ChatWorker(agent, "question")
        w.error.connect(errors.append)
        w.run()
        assert "No API key" in errors[0]

    def test_empty_reply_uses_placeholder(self, qapp):
        agent = SimpleNamespace(api_key="sk-x", query_llm=lambda p: "")
        replies = []
        w = ChatWorker(agent, "question")
        w.reply_done.connect(replies.append)
        w.run()
        assert replies == ["(no response from LLM)"]

    def test_query_exception_reports_error(self, qapp):
        def _boom(prompt):
            raise RuntimeError("LLM down")

        agent = SimpleNamespace(api_key="sk-x", query_llm=_boom)
        errors = []
        w = ChatWorker(agent, "question")
        w.error.connect(errors.append)
        w.run()
        assert errors == ["LLM down"]

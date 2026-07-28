# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.agents (AI agent execution panel).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — cached via ``pycsamt.app.web.cache.cache_set``
so the "processing" agent branch of ``_exec_agent`` runs a real emtools
function ("QC Quicklook" / "Dimensionality") end to end.

LLM agents are exercised against a fake class monkeypatched onto
``pycsamt.agents`` (a real LLM call needs a network round trip + API
key), matching real ``AGENT_REGISTRY`` metadata so the function's own
kwarg-filtering logic runs unmodified.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from pycsamt.app.web.callbacks.agents import _exec_agent
from pycsamt.app.web.layout import IDs

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


@pytest.fixture
def cached_session(willy_sites):
    from pycsamt.app.web.cache import cache_set

    session_id = "test-agents-session"
    cache_set(session_id, willy_sites)
    return session_id


class _FakeLLMAgent:
    last_instance = None

    def __init__(self, recursive=None, on_dup=None, api_key=None):
        self.recursive = recursive
        self.on_dup = on_dup
        self.api_key = api_key
        _FakeLLMAgent.last_instance = self

    def execute(self, input_data):
        self.input_data = input_data
        return SimpleNamespace(
            status="success",
            summary="Loaded OK",
            warnings=["one warning"],
            llm_interpretation="x" * 500,  # exercises the [:400] truncation
        )


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _run_agent_fn(web_app):
    key = next(
        k
        for k in web_app.callback_map
        if f"{IDs.AGENT_LOG}.children" in k and f"{IDs.AGENT_RESULT}.src" in k
    )
    return _cb(web_app, key)


# ── _exec_agent: no data / processing mode ──────────────────────────────────


class TestExecAgentNoData:
    def test_no_cached_sites_returns_warning(self):
        log, src, summary, err = _exec_agent("QC Quicklook", {}, "no-such-session")
        assert "Load survey data" in log
        assert summary == ""
        assert "Load survey data" in err


class TestExecAgentProcessing:
    def test_direct_figure_result(self, cached_session):
        log, src, summary, err = _exec_agent("QC Quicklook", {}, cached_session)
        assert err == ""
        assert "Completed." in log
        assert summary.startswith("✓")
        assert src != ""

    def test_result_plot_fallback_when_result_not_renderable(self, cached_session):
        """'Dimensionality' returns a classification table (not a figure)
        but declares a separate result_plot -- exercises the fallback
        branch that calls that plot fn instead."""
        log, src, summary, err = _exec_agent(
            "Dimensionality",
            {"skew_th": 3.0, "ellipt_th": 0.2},
            cached_session,
        )
        assert err == ""
        assert summary.startswith("✓")

    def test_processing_exception_reported(self, cached_session, monkeypatch):
        import pycsamt.emtools as et

        def _boom(*a, **k):
            raise RuntimeError("plot boom")

        monkeypatch.setattr(et, "plot_qc_quicklook", _boom)
        log, src, summary, err = _exec_agent("QC Quicklook", {}, cached_session)
        assert err == "plot boom"
        assert "✗ Error: plot boom" in log
        assert summary == "Agent failed: plot boom"

    def test_direct_axes_result_uses_get_figure(self, cached_session, monkeypatch):
        import matplotlib.pyplot as plt

        import pycsamt.emtools as et

        fig, ax = plt.subplots()

        def _fake_qc(sites, **kw):
            return ax  # Axes, not Figure -> exercises get_figure() branch

        monkeypatch.setattr(et, "plot_qc_quicklook", _fake_qc)
        log, src, summary, err = _exec_agent("QC Quicklook", {}, cached_session)
        assert err == ""
        assert src != ""
        plt.close(fig)

    def test_result_plot_axes_fallback(self, cached_session, monkeypatch):
        import matplotlib.pyplot as plt

        import pycsamt.emtools as et

        fig, ax = plt.subplots()

        def _fake_dim(sites, **kw):
            return SimpleNamespace()  # not renderable

        def _fake_dim_plot(sites, **kw):
            return ax  # Axes -> exercises the nested get_figure() fallback

        monkeypatch.setattr(et, "classify_dimensionality", _fake_dim)
        monkeypatch.setattr(et, "plot_dimensionality_psection", _fake_dim_plot)
        log, src, summary, err = _exec_agent(
            "Dimensionality",
            {"skew_th": 3.0, "ellipt_th": 0.2},
            cached_session,
        )
        assert err == ""
        assert src != ""
        plt.close(fig)

    def test_result_plot_figure_fallback(self, cached_session, monkeypatch):
        import matplotlib.pyplot as plt

        import pycsamt.emtools as et

        fig, ax = plt.subplots()

        def _fake_dim(sites, **kw):
            return SimpleNamespace()

        def _fake_dim_plot(sites, **kw):
            return fig  # Figure directly -> exercises the savefig() branch

        monkeypatch.setattr(et, "classify_dimensionality", _fake_dim)
        monkeypatch.setattr(et, "plot_dimensionality_psection", _fake_dim_plot)
        log, src, summary, err = _exec_agent(
            "Dimensionality",
            {"skew_th": 3.0, "ellipt_th": 0.2},
            cached_session,
        )
        assert err == ""
        assert src != ""
        plt.close(fig)

    def test_no_result_plot_declared_leaves_result_empty(
        self, cached_session, monkeypatch
    ):
        import pycsamt.app.web.callbacks.agents as agents_mod

        def _fake_get_entry(name):
            return {
                "type": "processing",
                "fn_name": "plot_qc_quicklook",
                "result_plot": None,
            }

        monkeypatch.setattr(agents_mod, "get_entry", _fake_get_entry)

        import pycsamt.emtools as et

        monkeypatch.setattr(
            et, "plot_qc_quicklook", lambda sites, **kw: SimpleNamespace()
        )
        log, src, summary, err = _exec_agent("QC Quicklook", {}, cached_session)
        assert err == ""
        assert "Completed." in log

    def test_station_filter_forwarded_and_logged(self, cached_session, monkeypatch):
        captured = {}

        def _fake_qc(sites, **kw):
            captured.update(kw)
            return SimpleNamespace()  # no savefig/get_figure -> no result_plot

        import pycsamt.emtools as et

        monkeypatch.setattr(et, "plot_qc_quicklook", _fake_qc)
        log, src, summary, err = _exec_agent(
            "QC Quicklook", {}, cached_session, stations=["STA01", "STA02"]
        )
        assert err == ""
        assert "Station filter: 2 station(s) selected." in log


# ── _exec_agent: LLM mode ────────────────────────────────────────────────────


class TestExecAgentLLM:
    def test_llm_success(self, cached_session, monkeypatch):
        import pycsamt.agents as ag

        monkeypatch.setattr(ag, "MTLoaderAgent", _FakeLLMAgent, raising=False)
        log, src, summary, err = _exec_agent(
            "MT Loader", {"recursive": True}, cached_session
        )
        assert err == ""
        assert "Status: success" in log
        assert "⚠ one warning" in log
        assert summary == "Loaded OK" + "\n\n" + "x" * 400

    def test_llm_no_warnings_no_interpretation(self, cached_session, monkeypatch):
        import pycsamt.agents as ag

        class _PlainAgent(_FakeLLMAgent):
            def execute(self, input_data):
                self.input_data = input_data
                return SimpleNamespace(
                    status="success",
                    summary="",
                    warnings=None,
                    llm_interpretation=None,
                )

        monkeypatch.setattr(ag, "MTLoaderAgent", _PlainAgent, raising=False)
        log, src, summary, err = _exec_agent("MT Loader", {}, cached_session)
        assert err == ""
        assert "⚠" not in log
        assert summary == ""

    def test_llm_forwards_station_filter_in_exec_input(
        self, cached_session, monkeypatch
    ):
        import pycsamt.agents as ag

        monkeypatch.setattr(ag, "MTLoaderAgent", _FakeLLMAgent, raising=False)
        _exec_agent("MT Loader", {}, cached_session, stations=["STA01"])
        assert _FakeLLMAgent.last_instance.input_data["station_names"] == ["STA01"]

    def test_llm_exception_reported(self, cached_session, monkeypatch):
        import pycsamt.agents as ag

        class _BoomAgent(_FakeLLMAgent):
            def execute(self, input_data):
                raise RuntimeError("LLM crashed")

        monkeypatch.setattr(ag, "MTLoaderAgent", _BoomAgent, raising=False)
        log, src, summary, err = _exec_agent("MT Loader", {}, cached_session)
        assert err == "LLM crashed"
        assert summary == "Agent failed: LLM crashed"


# ── run_agent callback ───────────────────────────────────────────────────────


class TestRunAgentCallback:
    def test_no_agent_selected(self, web_app):
        fn = _run_agent_fn(web_app)
        log, src, summary, is_open, body = fn(1, "", "session-1")
        assert log == "No agent selected."
        assert is_open is False
        assert body == ""

    def test_no_session_id_opens_toast(self, web_app):
        fn = _run_agent_fn(web_app)
        log, src, summary, is_open, body = fn(1, "QC Quicklook", None)
        assert "Session not initialised" in log
        assert is_open is True
        assert "Session not initialised" in body

    def test_success_path(self, web_app, cached_session):
        fn = _run_agent_fn(web_app)
        log, src, summary, is_open, body = fn(1, "QC Quicklook", cached_session)
        assert is_open is False
        assert body == ""
        assert "Completed." in log

    def test_failure_path_opens_toast(self, web_app, cached_session, monkeypatch):
        import pycsamt.emtools as et

        def _boom(*a, **k):
            raise RuntimeError("callback boom")

        monkeypatch.setattr(et, "plot_qc_quicklook", _boom)
        fn = _run_agent_fn(web_app)
        log, src, summary, is_open, body = fn(1, "QC Quicklook", cached_session)
        assert is_open is True
        assert body == "callback boom"

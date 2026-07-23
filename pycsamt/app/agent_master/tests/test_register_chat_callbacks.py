# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the Dash-wired callbacks in register_chat() —
pycsamt.app.agent_master.callbacks.chat: send-message, polling, stop-job,
figure modal, and pin/trace UI actions."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")

import pycsamt.app.agent_master.callbacks.chat as C
from pycsamt.app.agent_master._ids import IDs


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"]
        and entry["inputs"][0]["id"] == input_id
        and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


def _set_pattern_trigger(component_id_dict, value=1, prop="n_clicks"):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    prop_id = json.dumps(component_id_dict, sort_keys=True) + f".{prop}"
    cc.context_value.set(
        AttributeDict(
            triggered_inputs=[{"prop_id": prop_id, "value": value}]
        )
    )


def _set_plain_trigger(component_id, value=1, prop="n_clicks"):
    import dash._callback_context as cc
    from dash._utils import AttributeDict

    prop_id = f"{component_id}.{prop}"
    cc.context_value.set(
        AttributeDict(
            triggered_inputs=[{"prop_id": prop_id, "value": value}]
        )
    )


def _clear_trigger():
    import dash._callback_context as cc

    cc.context_value.set({})


class _FakeThread:
    """Captures the target/args instead of actually starting a thread."""

    last_started = None

    def __init__(self, target, args, daemon):
        self.target = target
        self.args = args
        self.daemon = daemon

    def start(self):
        _FakeThread.last_started = self


@pytest.fixture(autouse=True)
def _no_real_threads(monkeypatch):
    monkeypatch.setattr(C.threading, "Thread", _FakeThread)
    _FakeThread.last_started = None
    # send_message reads ctx.triggered_id unconditionally (even outside
    # the stop-mode branch); default to a plain BTN_SEND click so tests
    # that don't care about the trigger source don't need to set it.
    _set_plain_trigger(IDs.BTN_SEND)
    yield
    _clear_trigger()


# ── send_message ─────────────────────────────────────────────────────────────


def _send_fn(agent_app):
    return _find(agent_app, IDs.BTN_SEND, "am-store-pending.data")


class TestSendMessageStopMode:
    def test_poll_disabled_false_stops_job(self, agent_app):
        fn = _send_fn(agent_app)
        result = fn(
            1, [None], "hi", [], {}, {}, {}, [], False, {"jid": "job-1"}
        )
        msgs, job, disabled, value, stored, pending = result
        assert disabled is True
        assert job == {}
        assert pending == {}
        assert "Task stopped" in str(msgs[-1])


class TestSendMessageChipClick:
    def test_chip_click_uses_prompt_text(self, agent_app):
        from pycsamt.app.agent_master.layout import _PROMPT_CHIPS

        _set_pattern_trigger({"index": 0, "type": "am-chip"})
        fn = _send_fn(agent_app)
        result = fn(
            None,
            [1],
            "",
            [],
            {"path": "/tmp/edis"},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        msgs, job, disabled, value, stored, pending = result
        assert stored[-1]["content"] == _PROMPT_CHIPS[0]


class TestSendMessageEmptyText:
    def test_prevent_update_on_blank_text(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _send_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, [None], "   ", [], {}, {}, {}, [], True, {})


class TestSendMessageAppLaunch:
    def test_desktop_launch_success(self, agent_app, monkeypatch):
        monkeypatch.setattr(C, "_ensure_desktop", lambda: (True, "note"))
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "open the desktop app",
            [],
            {},
            {},
            {},
            [],
            True,
            {},
        )
        msgs, job, disabled, value, stored, pending = result
        assert disabled is True
        assert "Launched the pyCSAMT desktop app" in stored[-1]["content"]
        assert _FakeThread.last_started is None  # no job thread spawned

    def test_mapview_launch(self, agent_app, monkeypatch):
        monkeypatch.setattr(
            C, "_ensure_mapview", lambda: "http://127.0.0.1:8770"
        )
        fn = _send_fn(agent_app)
        result = fn(
            1, [None], "open the map", [], {}, {}, {}, [], True, {}
        )
        _msgs, _job, _disabled, _value, stored, _pending = result
        assert "MapView" in stored[-1]["content"]

    def test_web_app_launch(self, agent_app, monkeypatch):
        monkeypatch.setattr(
            C, "_ensure_web_app", lambda: "http://127.0.0.1:8051"
        )
        fn = _send_fn(agent_app)
        result = fn(
            1, [None], "open web app", [], {}, {}, {}, [], True, {}
        )
        _msgs, _job, _disabled, _value, stored, _pending = result
        assert "Redirecting to web app" in stored[-1]["content"]


class TestSendMessageNoDataIntentGate:
    def test_question_text_starts_job_without_edi_guard(self, agent_app):
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "what does StaticShiftAgent do",
            [],
            {},  # no EDI loaded
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        msgs, job, disabled, value, stored, pending = result
        assert disabled is False
        assert "jid" in job
        assert _FakeThread.last_started is not None
        assert _FakeThread.last_started.target is C._run_agent


class TestSendMessageNoEdiGuard:
    def test_guard_blocks_when_no_data_available(self, agent_app, monkeypatch):
        monkeypatch.setattr(C, "_names_registry_line", lambda text: False)
        monkeypatch.setattr(C, "_session_has_data", lambda: False)
        fn = _send_fn(agent_app)
        result = fn(
            1, [None], "run qc pipeline", [], {}, {}, {}, [], True, {}
        )
        msgs, job, disabled, value, stored, pending = result
        assert disabled is True
        assert "No EDI dataset is loaded" in stored[-1]["content"]
        assert _FakeThread.last_started is None

    def test_guard_skipped_for_dataless_workflow(self, agent_app, monkeypatch):
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: "layered_model",
        )
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "build a layered model",
            [],
            {},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        _msgs, _job, _disabled, _value, stored, _pending = result
        assert "No EDI dataset is loaded" not in stored[-1]["content"]

    def test_guard_skipped_when_registry_names_line(
        self, agent_app, monkeypatch
    ):
        monkeypatch.setattr(C, "_names_registry_line", lambda text: True)
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "run qc on line L22PLT",
            [],
            {},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        _msgs, _job, _disabled, _value, stored, _pending = result
        assert "No EDI dataset is loaded" not in stored[-1]["content"]

    def test_guard_skipped_when_session_has_data(self, agent_app, monkeypatch):
        monkeypatch.setattr(C, "_names_registry_line", lambda text: False)
        monkeypatch.setattr(C, "_session_has_data", lambda: True)
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "run qc pipeline",
            [],
            {},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        _msgs, _job, _disabled, _value, stored, _pending = result
        assert "No EDI dataset is loaded" not in stored[-1]["content"]


class TestSendMessageLineDisambiguation:
    def test_ambiguous_ref_shows_picker(self, agent_app):
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "run qc on line 1",
            [],
            {"path": "/x", "groups": {"L1A": [], "L1B": []}},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        msgs, job, disabled, value, stored, pending = result
        assert pending["disambiguation"] == "lines"
        assert "am-line-waiting-bubble" in str(msgs[-1])
        assert _FakeThread.last_started is None

    def test_exact_ref_match_proceeds_without_picker(self, agent_app):
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "run qc on line L22PLT",
            [],
            {"path": "/x", "groups": {"L22PLT": ["a.edi"], "L18PLT": ["b.edi"]}},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        _msgs, job, disabled, _value, _stored, pending = result
        assert pending == {}
        assert "jid" in job

    def test_prep_workflow_without_line_shows_picker(self, agent_app):
        import pycsamt.agents._workflows as wf_mod

        monkeypatch_target = wf_mod.classify_workflow
        try:
            wf_mod.classify_workflow = (
                lambda text, default=None: "pre_inversion"
            )
            fn = _send_fn(agent_app)
            result = fn(
                1,
                [None],
                "prep occam2d inputs",
                [],
                {"path": "/x", "groups": {"L1": [], "L2": []}},
                {"provider": "offline"},
                {},
                [],
                True,
                {},
            )
        finally:
            wf_mod.classify_workflow = monkeypatch_target
        msgs, job, disabled, value, stored, pending = result
        assert pending["disambiguation"] == "lines"

    def test_prep_workflow_all_lines_skips_picker(self, agent_app, monkeypatch):
        import pycsamt.agents._workflows as wf_mod

        monkeypatch.setattr(
            wf_mod,
            "classify_workflow",
            lambda text, default=None: "pre_inversion",
        )
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "prep occam2d inputs for all lines",
            [],
            {"path": "/x", "groups": {"L1": [], "L2": []}},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        _msgs, _job, _disabled, _value, _stored, pending = result
        assert pending.get("disambiguation") != "lines"


class TestSendMessageParamNeeded:
    def test_shows_param_modal_waiting_bubble(self, agent_app, monkeypatch):
        monkeypatch.setattr(C, "_quick_workflow", lambda text: "denoise")
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "denoise the data",
            [],
            {"path": "/x", "groups": {}},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        msgs, job, disabled, value, stored, pending = result
        assert pending["workflow"] == "denoise"
        assert pending["text"] == "denoise the data"
        assert disabled is True
        assert _FakeThread.last_started is None


class TestSendMessageNormalFlow:
    def test_starts_background_job(self, agent_app, monkeypatch):
        monkeypatch.setattr(C, "_quick_workflow", lambda text: None)
        fn = _send_fn(agent_app)
        result = fn(
            1,
            [None],
            "run qc pipeline",
            [],
            {"path": "/x"},
            {"provider": "offline"},
            {},
            [],
            True,
            {},
        )
        msgs, job, disabled, value, stored, pending = result
        assert disabled is False
        assert "jid" in job
        assert value == ""
        assert _FakeThread.last_started is not None
        assert _FakeThread.last_started.target is C._run_agent


# ── poll_job ──────────────────────────────────────────────────────────────────


def _poll_fn(agent_app):
    return _find(agent_app, IDs.INTERVAL_POLL, "am-store-postproc.data")


class TestPollJob:
    def test_prevent_update_without_job_store(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _poll_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, None, [], {}, [])

    def test_prevent_update_without_jid(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _poll_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, {}, [], {}, [])

    def test_prevent_update_when_job_missing(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = _poll_fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, {"jid": "ghost"}, [], {}, [])

    def test_cancelled_job_halts_polling(self, agent_app):
        jid = C._new_job()
        C._update_job(jid, status="cancelled")
        fn = _poll_fn(agent_app)
        from dash import no_update

        result = fn(1, {"jid": jid}, [], {}, [])
        assert result[1] is True
        assert result[0] is no_update

    def test_running_job_updates_thinking_bubble(self, agent_app):
        jid = C._new_job()
        C._update_job(
            jid, steps=[{"label": "Loading", "status": "running"}]
        )
        # Dash serializes children to dicts over the wire; poll_job's
        # bubble-matching relies on that dict shape (`isinstance(child,
        # dict)`), not a live component object.
        thinking_dict = {"props": {"id": "am-thinking-bubble"}}
        fn = _poll_fn(agent_app)
        msgs, disabled, figs, stored, postproc = fn(
            1, {"jid": jid}, [thinking_dict], {}, []
        )
        assert disabled is False
        assert "Loading" in str(msgs[0])

    def test_running_job_without_thinking_bubble_leaves_msgs(self, agent_app):
        jid = C._new_job()
        C._update_job(
            jid, steps=[{"label": "Loading", "status": "running"}]
        )
        fn = _poll_fn(agent_app)
        msgs, disabled, _figs, _stored, _postproc = fn(
            1, {"jid": jid}, ["unrelated"], {}, []
        )
        assert msgs == ["unrelated"]

    def test_done_job_appends_agent_bubble_and_cleans_up(self, agent_app):
        jid = C._new_job()
        C._update_job(
            jid,
            status="done",
            result="All good",
            steps=[],
            figs={"f1": {"title": "Fig", "b64": "aaa"}},
            kind=C.KIND_WORKFLOW,
        )
        fn = _poll_fn(agent_app)
        msgs, disabled, fig_store, stored, postproc = fn(
            1, {"jid": jid}, [], {}, []
        )
        assert disabled is True
        assert "f1" in fig_store
        assert stored[-1]["content"] == "All good"
        assert C._get_job(jid) is None  # cleaned up

    def test_done_job_with_postproc_returns_it(self, agent_app):
        jid = C._new_job()
        C._update_job(
            jid,
            status="done",
            result="Corrected",
            steps=[],
            postproc={"jid": jid, "workflow": "static_shift"},
        )
        fn = _poll_fn(agent_app)
        _msgs, _disabled, _figs, _stored, postproc = fn(
            1, {"jid": jid}, [], {}, []
        )
        assert postproc["workflow"] == "static_shift"

    def test_error_job_uses_error_text_as_result(self, agent_app):
        jid = C._new_job()
        C._update_job(
            jid, status="error", error="boom", result=None, steps=[]
        )
        fn = _poll_fn(agent_app)
        _msgs, _disabled, _figs, stored, _postproc = fn(
            1, {"jid": jid}, [], {}, []
        )
        assert stored[-1]["content"] == "boom"


# ── toggle_send_stop ─────────────────────────────────────────────────────────


class TestToggleSendStop:
    def _fn(self, agent_app):
        return _find(agent_app, IDs.INTERVAL_POLL, "am-btn-send.children")

    def test_running_shows_stop(self, agent_app):
        fn = self._fn(agent_app)
        _icon, cls, title = fn(False)
        assert cls == "am-send-stop"
        assert title == "Stop task"

    def test_idle_shows_send(self, agent_app):
        fn = self._fn(agent_app)
        _icon, cls, title = fn(True)
        assert cls == ""
        assert title == "Send"


# ── open_fig_modal ───────────────────────────────────────────────────────────


class TestOpenFigModal:
    def teardown_method(self, _method):
        _clear_trigger()

    def _fn(self, agent_app):
        matches = [
            k
            for k, entry in agent_app.callback_map.items()
            if entry["inputs"]
            and "am-fig-open" in entry["inputs"][0]["id"]
        ]
        assert len(matches) == 1, matches
        return _unwrap(agent_app.callback_map[matches[0]])

    def test_prevent_update_without_trigger(self, agent_app):
        from dash.exceptions import PreventUpdate

        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None], {})

    def test_prevent_update_for_missing_fig_key(self, agent_app):
        from dash.exceptions import PreventUpdate

        _set_pattern_trigger({"key": "ghost", "type": "am-fig-open"})
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([1], {"other": {}})

    def test_opens_with_figure_data(self, agent_app):
        _set_pattern_trigger({"key": "f1", "type": "am-fig-open"})
        fn = self._fn(agent_app)
        is_open, src, title, key = fn(
            [1], {"f1": {"b64": "aaa", "title": "My Fig"}}
        )
        assert is_open is True
        assert src == "data:image/png;base64,aaa"
        assert title == "My Fig"
        assert key == "f1"


# ── toggle_pin / unpin ───────────────────────────────────────────────────────


class TestTogglePin:
    def teardown_method(self, _method):
        _clear_trigger()

    def _fn(self, agent_app):
        matches = [
            k
            for k, entry in agent_app.callback_map.items()
            if entry["inputs"]
            and "am-pin-btn" in entry["inputs"][0]["id"]
        ]
        assert len(matches) == 1, matches
        return _unwrap(agent_app.callback_map[matches[0]])

    def test_prevent_update_without_trigger(self, agent_app):
        from dash.exceptions import PreventUpdate

        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None], [], [])

    def test_prevent_update_when_message_not_found(self, agent_app):
        from dash.exceptions import PreventUpdate

        _set_pattern_trigger({"mid": "ghost", "type": "am-pin-btn"})
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([1], [], [])

    def test_pins_a_known_message(self, agent_app):
        _set_pattern_trigger({"mid": "am-msg-1", "type": "am-pin-btn"})
        fn = self._fn(agent_app)
        messages = [
            {"mid": "am-msg-1", "role": "user", "content": "hi", "ts": "t"}
        ]
        pins = fn([1], [], messages)
        assert pins[0]["mid"] == "am-msg-1"


class TestUnpin:
    def teardown_method(self, _method):
        _clear_trigger()

    def _fn(self, agent_app):
        matches = [
            k
            for k, entry in agent_app.callback_map.items()
            if entry["inputs"]
            and "am-unpin" in entry["inputs"][0]["id"]
        ]
        assert len(matches) == 1, matches
        return _unwrap(agent_app.callback_map[matches[0]])

    def test_prevent_update_without_trigger(self, agent_app):
        from dash.exceptions import PreventUpdate

        import dash._callback_context as cc
        from dash._utils import AttributeDict

        cc.context_value.set(AttributeDict(triggered_inputs=[]))
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn([None], [])

    def test_removes_pin(self, agent_app):
        _set_pattern_trigger({"mid": "am-msg-1", "type": "am-unpin"})
        fn = self._fn(agent_app)
        pins = [{"mid": "am-msg-1"}, {"mid": "am-msg-2"}]
        result = fn([1], pins)
        assert result == [{"mid": "am-msg-2"}]


# ── render_pins / render_recent_runs ─────────────────────────────────────────


class TestRenderPins:
    def _fn(self, agent_app):
        return _find(agent_app, IDs.STORE_PINS, IDs.SIDEBAR_PINS)

    def test_empty_state(self, agent_app):
        fn = self._fn(agent_app)
        result = fn([])
        assert "No pinned messages" in str(result)

    def test_renders_items(self, agent_app):
        fn = self._fn(agent_app)
        result = fn(
            [{"mid": "m1", "role": "user", "snippet": "hi", "ts": "t"}]
        )
        assert len(result) == 1


class TestRenderRecentRuns:
    def _fn(self, agent_app):
        return _find(agent_app, IDs.CHAT_WINDOW, IDs.SIDEBAR_RUNS)

    def test_empty_state(self, agent_app, monkeypatch):
        monkeypatch.setattr(C, "_recent_runs", lambda n=8: [])
        fn = self._fn(agent_app)
        result = fn([])
        assert "No workflows run yet" in str(result)

    def test_renders_runs(self, agent_app, monkeypatch):
        monkeypatch.setattr(
            C,
            "_recent_runs",
            lambda n=8: [{"workflow": "qc", "status": "success"}],
        )
        fn = self._fn(agent_app)
        result = fn([])
        assert len(result) == 1

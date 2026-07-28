# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.params: field/form
renderers, param collectors, and the registered modal callbacks."""

from __future__ import annotations

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


# ── _field_el ────────────────────────────────────────────────────────────────


class TestFieldEl:
    def test_slider(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {
            "id": "fid-slider",
            "type": "slider",
            "min": 0,
            "max": 10,
            "step": 1,
        }
        el = _field_el(f, 5)
        assert el.id == "fid-slider"
        assert el.value == 5

    def test_number(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {"id": "fid-num", "type": "number"}
        el = _field_el(f, 3.5)
        assert el.value == 3.5

    def test_select(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {
            "id": "fid-sel",
            "type": "select",
            "options": [{"label": "A", "value": "a"}],
        }
        el = _field_el(f, "a")
        assert el.value == "a"
        assert el.clearable is False

    def test_multiselect_wraps_scalar_value(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {"id": "fid-multi", "type": "multiselect", "options": []}
        el = _field_el(f, "solo")
        assert el.value == ["solo"]

    def test_multiselect_keeps_list_value(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {"id": "fid-multi", "type": "multiselect", "options": []}
        el = _field_el(f, ["a", "b"])
        assert el.value == ["a", "b"]

    def test_multiselect_empty_value_becomes_empty_list(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {"id": "fid-multi", "type": "multiselect", "options": []}
        el = _field_el(f, None)
        assert el.value == []

    def test_radio(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {
            "id": "fid-radio",
            "type": "radio",
            "options": [{"label": "On", "value": "on"}],
        }
        el = _field_el(f, "on")
        assert el.value == "on"

    def test_text(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {"id": "fid-text", "type": "text"}
        el = _field_el(f, None)
        assert el.value == ""

    def test_unknown_type_falls_back_to_div(self):
        from pycsamt.app.agent_master.callbacks.params import _field_el

        f = {"id": "fid-div", "type": "totally-unknown"}
        el = _field_el(f, "whatever")
        assert el.id == "fid-div"


# ── _render_form ─────────────────────────────────────────────────────────────


class TestRenderForm:
    def test_renders_row_per_field_with_help_text(self):
        from pycsamt.app.agent_master.callbacks.params import _render_form

        fields = [
            {
                "id": "f1",
                "key": "k1",
                "label": "Field 1",
                "type": "text",
                "default": "",
                "help": "helpful",
            },
            {
                "id": "f2",
                "key": "k2",
                "label": "Field 2",
                "type": "number",
                "default": 1,
            },
        ]
        rows = _render_form(fields, {"k1": "hi"})
        assert len(rows) == 2
        assert "helpful" in str(rows[0])


# ── _render_steps_accordion ───────────────────────────────────────────────────


class TestRenderStepsAccordion:
    def test_unknown_steps_are_skipped(self):
        from pycsamt.app.agent_master.callbacks.params import (
            _render_steps_accordion,
        )

        result = _render_steps_accordion(["not-a-real-step"], {})
        assert result is None

    def test_known_step_renders_accordion(self):
        from pycsamt.app.agent_master.callbacks.params import (
            _render_steps_accordion,
        )

        result = _render_steps_accordion(["load"], {})
        assert result is not None
        assert "Pipeline Steps" in str(result)

    def test_step_params_merge_over_top_level_config(self):
        from pycsamt.app.agent_master.callbacks.params import (
            _render_steps_accordion,
        )

        cfg = {
            "period_min": 0.5,
            "step_params": {"load": {"period_min": 0.9}},
        }
        result = _render_steps_accordion(["load"], cfg)
        assert result is not None


# ── _collect_params / _collect_step_params ────────────────────────────────────


class TestCollectParams:
    def test_uses_field_values_and_keeps_other_inv_config(self):
        from pycsamt.app.agent_master.callbacks.params import (
            _collect_params,
        )

        schema = {
            "fields": [
                {"key": "n_layers", "default": 10},
                {"key": "epochs", "default": 500},
            ]
        }
        result = _collect_params(schema, {"n_layers": 15}, {"unrelated": "kept"})
        assert result["n_layers"] == 15
        assert result["epochs"] == 500  # falls back to default
        assert result["unrelated"] == "kept"


class TestCollectStepParams:
    def test_collects_known_step_fields_with_defaults(self):
        from pycsamt.app.agent_master.callbacks.params import (
            _collect_step_params,
        )

        schema = {"steps": ["load"]}
        result = _collect_step_params(schema, {"period_min": 0.01})
        assert result["load"]["period_min"] == 0.01
        assert result["load"]["period_max"] == 1.0  # default
        assert result["load"]["component"] == "xy"  # default

    def test_unknown_step_is_skipped(self):
        from pycsamt.app.agent_master.callbacks.params import (
            _collect_step_params,
        )

        schema = {"steps": ["nope"]}
        result = _collect_step_params(schema, {})
        assert result == {}


# ── _cancel_bubble ───────────────────────────────────────────────────────────


def test_cancel_bubble_mentions_cancellation():
    from pycsamt.app.agent_master.callbacks.params import _cancel_bubble

    bubble = _cancel_bubble()
    assert "Workflow cancelled" in str(bubble)


# ── _build_correction_schemas ─────────────────────────────────────────────────


def test_build_correction_schemas_handles_missing_module(monkeypatch):
    import sys

    from pycsamt.app.agent_master.callbacks import params as params_mod

    monkeypatch.setitem(sys.modules, "pycsamt.agents._corrections", None)
    # Should not raise even if the optional module is unavailable.
    params_mod._build_correction_schemas()


def test_build_correction_schemas_populates_schemas(monkeypatch):
    from pycsamt.app.agent_master.callbacks import params as params_mod

    params_mod._build_correction_schemas()
    # Real correction workflows (e.g. static-shift family) get schemas.
    assert isinstance(params_mod._SCHEMAS, dict)


# ── registered callbacks ─────────────────────────────────────────────────────


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _find(agent_app, input_id, output_hint):
    matches = [
        k
        for k, entry in agent_app.callback_map.items()
        if entry["inputs"] and entry["inputs"][0]["id"] == input_id and output_hint in k
    ]
    assert len(matches) == 1, (input_id, output_hint, matches)
    return _unwrap(agent_app.callback_map[matches[0]])


class TestOpenParamModal:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.STORE_PENDING, IDs.MODAL_PARAMS)

    def test_prevent_update_without_pending(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, {}, {})

    def test_prevent_update_for_unknown_workflow(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn({"workflow": "not-a-real-workflow"}, {}, {})

    def test_opens_with_fields_and_steps(self, agent_app):
        fn = self._fn(agent_app)
        is_open, title, desc, form, line_to_st = fn({"workflow": "denoise"}, {}, {})
        assert is_open is True
        assert "Data Denoising" in str(title) or "Denoising" in str(desc)
        assert len(form) > 0

    def test_opens_with_line_groups_for_station_field(self, agent_app):
        fn = self._fn(agent_app)
        is_open, _title, _desc, form, line_to_st = fn(
            {"workflow": "rhophi"},
            {},
            {"groups": {"L1": ["a.edi", "b.edi"], "L2": ["c.edi"]}},
        )
        assert is_open is True
        assert line_to_st  # populated because the schema uses "stations"


class TestStationsForLines:
    def _fn(self, agent_app):
        matches = [
            k
            for k, entry in agent_app.callback_map.items()
            if entry["inputs"]
            and entry["inputs"][0]["id"] == '{"key":"lines","type":"am-pf"}'
        ]
        assert len(matches) == 1, matches
        return _unwrap(agent_app.callback_map[matches[0]])

    def test_filters_by_selected_lines(self, agent_app):
        fn = self._fn(agent_app)
        line_to_st = {"L1": ["a", "b"], "L2": ["c"]}
        opts = fn(["L1"], line_to_st)
        assert opts == [
            {"label": "a", "value": "a"},
            {"label": "b", "value": "b"},
        ]

    def test_empty_selection_returns_all_stations(self, agent_app):
        fn = self._fn(agent_app)
        line_to_st = {"L1": ["a"], "L2": ["c"]}
        opts = fn([], line_to_st)
        assert {o["value"] for o in opts} == {"a", "c"}


class TestCancelParams:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_PARAM_CANCEL, IDs.MODAL_PARAMS)

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, [])

    def test_replaces_waiting_bubble(self, agent_app):
        waiting_bubble = {"props": {"id": "am-waiting-bubble"}}
        fn = self._fn(agent_app)
        is_open, pending, msgs = fn(1, [waiting_bubble])
        assert is_open is False
        assert pending == {}
        assert msgs[0] != waiting_bubble

    def test_appends_when_no_waiting_bubble_present(self, agent_app):
        fn = self._fn(agent_app)
        is_open, pending, msgs = fn(1, [])
        assert is_open is False
        assert len(msgs) == 1


class TestSubmitParams:
    def _fn(self, agent_app):
        from pycsamt.app.agent_master._ids import IDs

        return _find(agent_app, IDs.BTN_PARAM_RUN, IDs.CHAT_WINDOW)

    def _set_states_list(self, pf_pairs, ps_pairs):
        import dash._callback_context as cc
        from dash._utils import AttributeDict

        pf_list = [{"id": {"type": "am-pf", "key": k}, "value": v} for k, v in pf_pairs]
        ps_list = [{"id": {"type": "am-ps", "key": k}, "value": v} for k, v in ps_pairs]
        cc.context_value.set(
            AttributeDict(
                states_list=[
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    pf_list,
                    ps_list,
                ]
            )
        )

    def teardown_method(self, _method):
        import dash._callback_context as cc

        cc.context_value.set({})

    def test_prevent_update_without_click(self, agent_app):
        from dash.exceptions import PreventUpdate

        self._set_states_list([], [])
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(None, {"workflow": "denoise"}, [], {}, {}, [], {}, [], [])

    def test_prevent_update_without_pending(self, agent_app):
        from dash.exceptions import PreventUpdate

        self._set_states_list([], [])
        fn = self._fn(agent_app)
        with pytest.raises(PreventUpdate):
            fn(1, None, [], {}, {}, [], {}, [], [])

    def test_starts_job_and_persists_config(self, agent_app, monkeypatch):
        from pycsamt.app.agent_master.callbacks import params as params_mod

        monkeypatch.setattr(params_mod, "_new_job", lambda: "job-xyz")

        started = {}

        class _FakeThread:
            def __init__(self, target, args, daemon):
                started["target"] = target
                started["args"] = args

            def start(self):
                started["started"] = True

        monkeypatch.setattr(params_mod.threading, "Thread", _FakeThread)

        self._set_states_list(
            [("denoise_method", "median"), ("snr_threshold", 0.8)], []
        )
        fn = self._fn(agent_app)
        result = fn(
            1,
            {"text": "denoise it", "workflow": "denoise"},
            [],
            {"path": "/edi"},
            {},
            [],
            {"stale": "kept"},
            [],
            [],
        )
        msgs, job, disabled, stored, persist_ic, pending, is_open = result
        assert started.get("started") is True
        assert job == {"jid": "job-xyz"}
        assert disabled is False
        assert persist_ic["denoise_method"] == "median"
        assert "workflow" not in persist_ic
        assert pending == {}
        assert is_open is False

    def test_recovers_edi_path_from_pending_when_store_empty(
        self, agent_app, monkeypatch
    ):
        from pycsamt.app.agent_master.callbacks import params as params_mod

        monkeypatch.setattr(params_mod, "_new_job", lambda: "job-2")
        captured = {}

        class _FakeThread:
            def __init__(self, target, args, daemon):
                captured["args"] = args

            def start(self):
                pass

        monkeypatch.setattr(params_mod.threading, "Thread", _FakeThread)

        self._set_states_list([], [])
        fn = self._fn(agent_app)
        fn(
            1,
            {
                "text": "run",
                "workflow": "denoise",
                "edi_path": "/fallback/path",
                "edi_groups": {"L1": ["a.edi"]},
                "selected_lines": ["L1"],
            },
            [],
            {},  # STORE_EDI empty -> triggers fallback recovery
            {},
            [],
            {},
            [],
            [],
        )
        _jid, _text, edi_use, _settings, _new_ic = captured["args"]
        assert edi_use["path"] == "/fallback/path"
        assert edi_use["groups"] == {"L1": ["a.edi"]}
        assert edi_use["selected_lines"] == ["L1"]

    def test_replaces_waiting_bubble_with_thinking(self, agent_app, monkeypatch):
        from pycsamt.app.agent_master.callbacks import params as params_mod

        monkeypatch.setattr(params_mod, "_new_job", lambda: "job-3")

        class _FakeThread:
            def __init__(self, target, args, daemon):
                pass

            def start(self):
                pass

        monkeypatch.setattr(params_mod.threading, "Thread", _FakeThread)

        self._set_states_list([], [])
        waiting_bubble = {"props": {"id": "am-waiting-bubble"}}
        fn = self._fn(agent_app)
        result = fn(
            1,
            {"text": "run", "workflow": "denoise"},
            [waiting_bubble],
            {"path": "/edi"},
            {},
            [],
            {},
            [],
            [],
        )
        msgs = result[0]
        assert msgs[0] != waiting_bubble

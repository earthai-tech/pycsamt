# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.pages.agents_page (the AI Agents Dash page).

Strategy
--------
Pure builder functions (``_smart_route``, ``_param_widget``,
``_build_param_form``, ``_agent_info_block``, ``_grouped_names``,
``_agent_list``, ``_render_history``, ``_build_data_bar``,
``_welcome_message``, ``_render_chat``, ``layout``) are exercised
directly against the real ``AGENT_REGISTRY`` from
``pycsamt.app.desktop.agent_registry`` -- no mocking needed since they
are deterministic pure functions of their arguments.

``register_callbacks`` uses the bare ``dash.callback`` decorator. Its
``app`` argument is never referenced inside the function body (only
``from dash import callback`` is used), so callbacks are captured by
monkeypatching ``dash.callback`` with a spy decorator that records each
inner function by name and returns it unchanged -- this sidesteps the
"invisible until ``Dash._setup_server()`` runs" global-registry gotcha
entirely and lets every callback be invoked as a plain Python function.
Dash's ``ctx`` (``dash.callback_context``) is monkeypatched per-test on
the ``agents_page`` module directly with a ``SimpleNamespace`` exposing
only the attributes (``triggered_id``, ``states_list``, ``outputs_list``)
each callback actually reads.
"""

from __future__ import annotations

from types import SimpleNamespace

import dash
import pytest
from dash import no_update

from pycsamt.app.desktop.agent_registry import (
    AGENT_REGISTRY,
    agents_by_category,
    default_params,
    processing_agents,
)
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.pages import agents_page as mod
from pycsamt.app.web.utils import empty_src

# ── shared helpers ───────────────────────────────────────────────────────────


def _iter_ids(node):
    ids = []

    def _rec(n):
        if n is None:
            return
        if isinstance(n, (list, tuple)):
            for c in n:
                _rec(c)
            return
        cid = getattr(n, "id", None)
        if cid is not None:
            ids.append(cid)
        _rec(getattr(n, "children", None))

    _rec(node)
    return ids


def _rec(id_, line="—"):
    return {
        "ID": id_,
        "Latitude": None,
        "Longitude": None,
        "Elevation": None,
        "N_freq": 10,
        "Tipper": False,
        "Line": line,
    }


@pytest.fixture
def cbs(monkeypatch):
    """Capture every callback closure registered by ``register_callbacks``."""
    captured = {}

    def fake_callback(*_args, **_kwargs):
        def _decorator(fn):
            captured[fn.__name__] = fn
            return fn

        return _decorator

    monkeypatch.setattr(dash, "callback", fake_callback)
    mod.register_callbacks(None)
    return captured


# ── _smart_route ─────────────────────────────────────────────────────────────


class TestSmartRoute:
    def test_qc_keyword_match(self):
        plan = mod._smart_route("please run a QC check for noise")
        assert plan["via"] == "keyword"
        assert set(plan["agents"]) <= {"QC Analysis", "QC Quicklook"}
        assert plan["agents"]

    def test_static_shift_keyword(self):
        plan = mod._smart_route("correct the static shift")
        assert "Static Shift" in plan["agents"]

    def test_model_matches_occam_group_before_forward_group(self):
        # "model" appears in both the occam/inversion group and the
        # forward/synthetic group; the occam group is listed first in
        # _KEYWORD_MAP so it must win.
        plan = mod._smart_route("build a resistivity model")
        assert plan["agents"] == ["Inversion Prep", "Occam2D"]

    def test_no_match_falls_back_to_help_text(self):
        plan = mod._smart_route("xyzzy plugh unrelated gibberish")
        assert plan["agents"] == []
        assert "couldn't identify" in plan["description"]

    def test_case_insensitive(self):
        plan = mod._smart_route("PLEASE INTERPRET RESULTS")
        assert "Interpretation" in plan["agents"]

    def test_registry_filtering_drops_unregistered_candidates(self, monkeypatch):
        monkeypatch.setattr(
            mod,
            "_KEYWORD_MAP",
            [(["zzz"], ["Not A Real Agent", "QC Quicklook"])],
        )
        plan = mod._smart_route("zzz")
        assert plan["agents"] == ["QC Quicklook"]

    def test_registry_filtering_all_candidates_missing_falls_through(
        self, monkeypatch
    ):
        monkeypatch.setattr(
            mod, "_KEYWORD_MAP", [(["zzz"], ["Not A Real Agent"])]
        )
        plan = mod._smart_route("zzz")
        assert plan["agents"] == []


# ── _param_widget ────────────────────────────────────────────────────────────


class TestParamWidget:
    def test_combo_with_default(self):
        w = mod._param_widget(
            "method", {"type": "combo", "options": ["a", "b"], "default": "b"}
        )
        assert w.value == "b"
        assert w.options == [{"label": "a", "value": "a"}, {"label": "b", "value": "b"}]

    def test_combo_default_none_uses_first_option(self):
        w = mod._param_widget(
            "method", {"type": "combo", "options": ["a", "b"], "default": None}
        )
        assert w.value == "a"

    def test_bool_widget(self):
        w = mod._param_widget("flag", {"type": "bool", "default": True})
        assert w.value is True

    def test_bool_widget_default_none_is_falsy(self):
        w = mod._param_widget("flag", {"type": "bool", "default": None})
        assert w.value is False

    def test_float_widget_range_and_step(self):
        w = mod._param_widget(
            "x",
            {
                "type": "float",
                "default": 0.5,
                "range": (0.0, 1.0),
                "step": 0.1,
            },
        )
        assert w.type == "number"
        assert w.value == 0.5
        assert w.min == 0.0
        assert w.max == 1.0
        assert w.step == 0.1

    def test_int_widget_no_range_defaults_to_none(self):
        w = mod._param_widget("n", {"type": "int", "default": 5})
        assert w.min is None
        assert w.max is None

    def test_text_widget_with_default(self):
        w = mod._param_widget("s", {"type": "str", "default": "hello"})
        assert w.type == "text"
        assert w.value == "hello"

    def test_text_widget_default_none_becomes_empty_string(self):
        w = mod._param_widget("s", {"type": "str", "default": None})
        assert w.value == ""

    def test_widget_id_is_pattern_matched(self):
        w = mod._param_widget("k", {"type": "str", "default": ""})
        assert w.id == {"type": "agent-param", "index": "k"}


# ── _build_param_form ────────────────────────────────────────────────────────


class TestBuildParamForm:
    def test_no_name_returns_no_selection_message(self):
        out = mod._build_param_form(None)
        assert out.children == "No agent selected."

    def test_unknown_name_returns_no_selection_message(self, monkeypatch):
        monkeypatch.setattr(mod, "get_entry", lambda name: None)
        out = mod._build_param_form("nope")
        assert out.children == "No agent selected."

    def test_no_params_agent(self):
        out = mod._build_param_form("QC Quicklook")
        assert "No configurable parameters." in out.children

    def test_with_params_builds_rows_with_range_hint(self):
        out = mod._build_param_form("QC Analysis")
        rows = out.children
        assert len(rows) == len(AGENT_REGISTRY["QC Analysis"]["params"])
        label_row = rows[0]
        label, widget = label_row.children
        assert "Min fraction OK" not in label.children  # first param is "method" (combo)
        assert label.children == "QC method"

    def test_float_param_shows_range_hint(self):
        out = mod._build_param_form("QC Analysis")
        rows = out.children
        min_frac_row = rows[1]
        label, _widget = min_frac_row.children
        assert "[0.0" in label.children and "1.0]" in label.children


# ── _agent_info_block ────────────────────────────────────────────────────────


class TestAgentInfoBlock:
    def test_no_name(self):
        out = mod._agent_info_block(None)
        assert len(out) == 1
        assert out[0].children == "No agent selected."

    def test_unknown_name(self, monkeypatch):
        monkeypatch.setattr(mod, "get_entry", lambda name: None)
        out = mod._agent_info_block("nope")
        assert out[0].children == "No agent selected."

    def test_llm_agent_with_category_and_description(self):
        out = mod._agent_info_block("QC Analysis")
        pills_row, desc_p = out
        badge, cat_pill = pills_row.children
        assert badge.className == "agent-type-badge agent-badge-llm"
        assert cat_pill.children == "Quality Control"
        assert "Quality-control scoring" in desc_p.children

    def test_processing_agent_no_category_no_desc_pill(self, monkeypatch):
        monkeypatch.setattr(
            mod,
            "get_entry",
            lambda name: {"type": "processing", "params": {}, "description": ""},
        )
        out = mod._agent_info_block("Fake Proc")
        # No description paragraph, and only the type badge (no category pill)
        assert len(out) == 1
        pills_row = out[0]
        assert len(pills_row.children) == 1

    def test_unknown_type_falls_back_to_default_badge(self, monkeypatch):
        monkeypatch.setattr(
            mod,
            "get_entry",
            lambda name: {"type": "weird", "description": "", "category": ""},
        )
        out = mod._agent_info_block("Weird Agent")
        badge = out[0].children[0]
        assert badge.className == "agent-type-badge agent-badge-processing"
        assert badge.children[1] == "Weird"

    def test_category_fallback_to_name_to_cat_map(self, monkeypatch):
        monkeypatch.setattr(
            mod,
            "get_entry",
            lambda name: {"type": "llm", "description": ""},
        )
        out = mod._agent_info_block("QC Analysis")
        pills_row = out[0]
        _badge, cat_pill = pills_row.children
        assert cat_pill.children == "Quality Control"


# ── _grouped_names ───────────────────────────────────────────────────────────


class TestGroupedNames:
    def test_processing_grouped_last(self):
        names = ["QC Quicklook", "QC Analysis", "Static Shift (fast)"]
        grouped = mod._grouped_names(names)
        sections = [g[0] for g in grouped]
        assert sections[-1] == "⚡ Processing"
        proc_names = dict(grouped)["⚡ Processing"]
        assert set(proc_names) == {"QC Quicklook", "Static Shift (fast)"}

    def test_llm_categories_preserved(self):
        names = ["QC Analysis", "Static Shift", "Interpretation"]
        grouped = mod._grouped_names(names)
        cats = dict(grouped)
        assert cats["Quality Control"] == ["QC Analysis"]
        assert cats["Pre-processing"] == ["Static Shift"]

    def test_unknown_name_falls_back_to_other_category(self):
        grouped = mod._grouped_names(["Totally Unknown Agent"])
        assert dict(grouped)["Other"] == ["Totally Unknown Agent"]

    def test_empty_names_returns_empty_list(self):
        assert mod._grouped_names([]) == []


# ── _agent_list ──────────────────────────────────────────────────────────────


class TestAgentList:
    def test_empty_names_shows_placeholder(self):
        out = mod._agent_list([], None)
        assert out.children == "No agents found."

    def test_selected_row_gets_highlighted(self):
        names = ["QC Quicklook", "QC Analysis"]
        out = mod._agent_list(names, "QC Analysis")
        buttons = [c for c in out.children if getattr(c, "n_clicks", None) == 0]
        selected = {b.id["index"]: b.className for b in buttons}
        assert "selected" in selected["QC Analysis"]
        assert "selected" not in selected["QC Quicklook"]

    def test_no_selection_defaults_to_first_name(self):
        names = ["QC Quicklook", "QC Analysis"]
        out = mod._agent_list(names, None)
        buttons = [c for c in out.children if getattr(c, "n_clicks", None) == 0]
        selected = {b.id["index"]: b.className for b in buttons}
        assert "selected" in selected["QC Quicklook"]

    def test_processing_badge_vs_llm_badge(self):
        names = ["QC Quicklook", "QC Analysis"]
        out = mod._agent_list(names, "QC Quicklook")
        buttons = [c for c in out.children if getattr(c, "n_clicks", None) == 0]
        by_name = {b.id["index"]: b for b in buttons}
        proc_badge = by_name["QC Quicklook"].children[1]
        llm_badge = by_name["QC Analysis"].children[1]
        assert proc_badge.children == "⚡"
        assert llm_badge.children == "LLM"


# ── _render_history ──────────────────────────────────────────────────────────


class TestRenderHistory:
    def test_empty_history(self):
        out = mod._render_history(None)
        assert len(out) == 1
        assert "No runs yet." in out[0].children

    def test_success_and_error_icons(self):
        history = [
            {"status": "success", "name": "QC Quicklook", "timestamp": "10:00"},
            {"status": "error", "name": "Static Shift", "timestamp": "10:01"},
        ]
        out = mod._render_history(history)
        icon_ok = out[0].children[0]
        icon_err = out[1].children[0]
        assert "history-success" in icon_ok.className
        assert "history-error" in icon_err.className

    def test_capped_at_15(self):
        history = [
            {"status": "success", "name": f"Agent {i}", "timestamp": "t"}
            for i in range(20)
        ]
        out = mod._render_history(history)
        assert len(out) == 15
        assert out[0].id == {"type": "history-item", "index": 0}
        assert out[14].id == {"type": "history-item", "index": 14}


# ── _build_data_bar ──────────────────────────────────────────────────────────


class TestBuildDataBar:
    def test_no_store_data(self):
        out = mod._build_data_bar(None)
        assert "No survey loaded" in out[1].children

    def test_zero_stations_treated_as_no_data(self):
        out = mod._build_data_bar({"n_stations": 0})
        assert "No survey loaded" in out[1].children

    def test_basic_station_and_line_chips(self):
        out = mod._build_data_bar({"n_stations": 5, "n_lines": 2})
        text = " ".join(
            c.children for c in out if isinstance(getattr(c, "children", None), str)
        )
        assert "5 stations" in text
        assert "2 profiles" in text

    def test_singular_station_and_profile(self):
        out = mod._build_data_bar({"n_stations": 1, "n_lines": 1})
        text = " ".join(
            c.children for c in out if isinstance(getattr(c, "children", None), str)
        )
        assert "1 station" in text and "1 stations" not in text
        assert "1 profile" in text and "1 profiles" not in text

    def test_line_counts_three_or_fewer_no_truncation(self):
        out = mod._build_data_bar(
            {"n_stations": 3, "n_lines": 3, "line_counts": {"L1": 1, "L2": 1, "L3": 1}}
        )
        summary_chip = out[-1]
        assert "more" not in summary_chip.children

    def test_line_counts_more_than_three_truncated(self):
        line_counts = {f"L{i}": 1 for i in range(5)}
        out = mod._build_data_bar(
            {"n_stations": 5, "n_lines": 5, "line_counts": line_counts}
        )
        summary_chip = out[-1]
        assert "+2 more" in summary_chip.children

    def test_data_dir_short_not_truncated(self):
        out = mod._build_data_bar({"n_stations": 1, "data_dir": "/short/path"})
        chip = out[-1]
        assert chip.children == "/short/path"

    def test_data_dir_long_truncated(self):
        long_dir = "/very/" + "x" * 50 + "/data"
        out = mod._build_data_bar({"n_stations": 1, "data_dir": long_dir})
        chip = out[-1]
        assert chip.children.startswith("…")
        assert chip.children == "…" + long_dir[-38:]

    def test_data_dir_uploaded_sentinel_skipped(self):
        out = mod._build_data_bar({"n_stations": 1, "data_dir": "[uploaded]"})
        # No path chip appended -- last chip is the profile-count chip
        assert "profile" in out[-1].children


# ── _welcome_message ─────────────────────────────────────────────────────────


class TestWelcomeMessage:
    def test_shape(self):
        msg = mod._welcome_message()
        assert msg["id"] == "welcome"
        assert msg["role"] == "assistant"
        assert msg["plan"] is None
        assert msg["plan_status"] is None
        assert "MT workflow assistant" in msg["content"]


# ── _render_chat ─────────────────────────────────────────────────────────────


class TestRenderChat:
    def _base_assistant_msg(self, **overrides):
        msg = {
            "id": "m1",
            "role": "assistant",
            "content": "",
            "plan": None,
            "plan_status": None,
            "result_log": None,
            "result_src": None,
            "result_summary": None,
            "timestamp": "10:00",
        }
        msg.update(overrides)
        return msg

    def test_user_bubble(self):
        out = mod._render_chat(
            [{"role": "user", "content": "hello", "timestamp": "10:00"}]
        )
        assert len(out) == 1
        assert out[0].className == "chat-bubble chat-bubble-user"

    def test_assistant_bubble_splits_paragraphs(self):
        msg = self._base_assistant_msg(content="para one\n\npara two")
        out = mod._render_chat([msg])
        bubble = out[0]
        texts = [
            c.children
            for c in bubble.children
            if getattr(c, "className", "") == "chat-bubble-text"
        ]
        assert texts == ["para one", "para two"]

    def test_plan_card_keyword_via(self):
        msg = self._base_assistant_msg(
            plan={"agents": ["QC Quicklook"], "via": "keyword"},
            plan_status="pending",
        )
        out = mod._render_chat([msg])
        bubble = out[0]
        plan_card = next(
            c for c in bubble.children if getattr(c, "className", "") == "chat-plan-card"
        )
        header = plan_card.children[0]
        via_note = header.children[-1]
        assert via_note.className == "chat-via-kw"

    def test_plan_card_llm_via(self):
        msg = self._base_assistant_msg(
            plan={"agents": ["QC Analysis"], "via": "llm"},
            plan_status=None,
        )
        out = mod._render_chat([msg])
        plan_card = next(
            c
            for c in out[0].children
            if getattr(c, "className", "") == "chat-plan-card"
        )
        header = plan_card.children[0]
        via_note = header.children[-1]
        assert via_note.className == "chat-via-llm"

    def test_plan_card_badges_llm_vs_processing_agent(self):
        msg = self._base_assistant_msg(
            plan={"agents": ["QC Analysis", "QC Quicklook"], "via": "keyword"},
            plan_status="pending",
        )
        out = mod._render_chat([msg])
        plan_card = next(
            c
            for c in out[0].children
            if getattr(c, "className", "") == "chat-plan-card"
        )
        agent_list = plan_card.children[1]
        llm_li, proc_li = agent_list.children
        assert llm_li.children[1].children == " LLM"
        assert proc_li.children[1].children == " ⚡"

    def test_plan_not_rendered_when_status_done(self):
        msg = self._base_assistant_msg(
            plan={"agents": ["QC Quicklook"], "via": "keyword"},
            plan_status="done",
        )
        out = mod._render_chat([msg])
        bubble = out[0]
        assert not any(
            getattr(c, "className", "") == "chat-plan-card" for c in bubble.children
        )

    def test_cancelled_note(self):
        msg = self._base_assistant_msg(plan_status="cancelled")
        out = mod._render_chat([msg])
        bubble = out[0]
        assert any(
            getattr(c, "className", "") == "chat-cancelled-note"
            for c in bubble.children
        )

    def test_running_spinner(self):
        msg = self._base_assistant_msg(plan_status="running")
        out = mod._render_chat([msg])
        bubble = out[0]
        running_divs = [
            c for c in bubble.children if "align-items-center" in getattr(c, "className", "")
        ]
        assert len(running_divs) == 1

    def test_inline_result_with_log_summary_and_image(self):
        msg = self._base_assistant_msg(
            result_log="a log line",
            result_summary="a summary",
            result_src="data:image/png;base64,AAAA",
        )
        out = mod._render_chat([msg])
        bubble = out[0]
        result_block = next(
            c for c in bubble.children if getattr(c, "className", "") == "chat-result-block"
        )
        kinds = {type(c).__name__ for c in result_block.children}
        assert "P" in kinds
        assert "Details" in kinds
        assert "Img" in kinds

    def test_inline_result_suppresses_placeholder_image(self):
        placeholder = empty_src(dark=True)
        msg = self._base_assistant_msg(result_log="a log line", result_src=placeholder)
        out = mod._render_chat([msg])
        bubble = out[0]
        result_block = next(
            c for c in bubble.children if getattr(c, "className", "") == "chat-result-block"
        )
        kinds = {type(c).__name__ for c in result_block.children}
        assert "Img" not in kinds

    def test_no_result_block_when_nothing_to_show(self):
        msg = self._base_assistant_msg()
        out = mod._render_chat([msg])
        bubble = out[0]
        assert not any(
            getattr(c, "className", "") == "chat-result-block" for c in bubble.children
        )

    def test_inline_result_summary_only_skips_log_details(self):
        msg = self._base_assistant_msg(result_summary="just a summary, no log")
        out = mod._render_chat([msg])
        bubble = out[0]
        result_block = next(
            c for c in bubble.children if getattr(c, "className", "") == "chat-result-block"
        )
        kinds = {type(c).__name__ for c in result_block.children}
        assert kinds == {"P"}


# ── layout ───────────────────────────────────────────────────────────────────


class TestLayout:
    def test_returns_a_div_with_flex_column_style(self):
        out = mod.layout()
        assert type(out).__name__ == "Div"
        assert out.style["display"] == "flex"
        assert out.style["flexDirection"] == "column"

    def test_contains_expected_ids(self):
        out = mod.layout()
        ids = _iter_ids(out)
        for expected in (
            IDs.AGENTS_CAT,
            IDs.AGENTS_NAME,
            IDs.AGENTS_DESC,
            IDs.AGENTS_PARAM_FORM,
            IDs.BTN_AGENTS_RUN,
            IDs.AGENTS_SPINNER,
            IDs.AGENTS_STORE,
            IDs.AGENTS_OUT,
            IDs.IMG_AGENTS,
            mod._ID_SEARCH,
            mod._ID_LIST_WRAP,
            mod._ID_DETAIL_NAME,
            mod._ID_HISTORY,
            mod._ID_SUMMARY,
            mod._ID_DATA_BAR,
            mod._ID_STATION_FILTER,
            mod._ID_VIEW_RUNNER,
            mod._ID_VIEW_CHAT,
            mod._ID_TAB_RUNNER,
            mod._ID_TAB_CHAT,
            mod._ID_LINE_FILTER,
        ):
            assert expected in ids, f"missing id {expected!r}"

    def test_chat_view_hidden_by_default(self):
        out = mod.layout()
        view_runner, view_chat = None, None
        for node in out.children:
            if getattr(node, "id", None) == mod._ID_VIEW_RUNNER:
                view_runner = node
            if getattr(node, "id", None) == mod._ID_VIEW_CHAT:
                view_chat = node
        assert view_chat.style["display"] == "none"
        assert view_runner.className == "agents-3col-layout"

    def test_runner_tab_active_by_default(self):
        out = mod.layout()
        toggle = out.children[1]
        runner_btn, chat_btn = toggle.children
        assert "active" in runner_btn.className
        assert "active" not in chat_btn.className


# ── register_callbacks: card / selection / filter ───────────────────────────


class TestCardToSelect:
    def test_dict_triggered_returns_index(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "agent-card", "index": "QC Analysis"})
        )
        assert cbs["_card_to_select"]([1]) == "QC Analysis"

    def test_non_dict_triggered_returns_no_update(self, cbs, monkeypatch):
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(triggered_id="something-else"))
        assert cbs["_card_to_select"]([1]) is no_update


class TestOnAgentSelected:
    def test_no_name_returns_all_no_update(self, cbs):
        out = cbs["_on_agent_selected"](None)
        assert out == (no_update, no_update, no_update, no_update)

    def test_highlights_matching_row(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod,
            "ctx",
            SimpleNamespace(
                outputs_list=[
                    [
                        {"id": {"index": "QC Quicklook"}},
                        {"id": {"index": "QC Analysis"}},
                    ]
                ]
            ),
        )
        classes, desc, form, name = cbs["_on_agent_selected"]("QC Analysis")
        assert classes == ["agent-list-row", "agent-list-row selected"]
        assert name == "QC Analysis"


class TestFilterList:
    def test_processing_filter(self, cbs):
        out, sel = cbs["_filter_list"]("__processing__", None)
        assert sel in processing_agents()

    def test_category_filter(self, cbs):
        out, sel = cbs["_filter_list"]("Quality Control", None)
        assert sel in agents_by_category()["Quality Control"]

    def test_no_category_uses_all_names(self, cbs):
        out, sel = cbs["_filter_list"]("", None)
        assert sel == mod._ALL_NAMES[0]

    def test_search_narrows_selection(self, cbs):
        out, sel = cbs["_filter_list"]("", "quicklook")
        assert sel == "QC Quicklook"

    def test_no_matches_selects_none(self, cbs):
        out, sel = cbs["_filter_list"]("", "zzzznomatch")
        assert sel is None
        assert out.children == "No agents found."


# ── register_callbacks: data status + line/station filters ─────────────────


class TestUpdateDataStatus:
    def test_no_data(self, cbs):
        bar, line_opts, stn_opts, disabled = cbs["_update_data_status"](None, None)
        assert line_opts == []
        assert stn_opts == []
        assert disabled is True

    def test_with_data_and_active_lines(self, cbs):
        store = {
            "n_stations": 2,
            "station_records": [_rec("S1", "L1"), _rec("S2", "L2")],
        }
        bar, line_opts, stn_opts, disabled = cbs["_update_data_status"](
            store, {"active": ["L1"]}
        )
        assert disabled is False
        assert line_opts == [{"label": "L1", "value": "L1"}]
        assert stn_opts == [{"label": "S1", "value": "S1"}]

    def test_active_lines_empty_auto_derives_from_records(self, cbs):
        store = {
            "n_stations": 2,
            "station_records": [_rec("S1", "L1"), _rec("S2", "L2")],
        }
        bar, line_opts, stn_opts, disabled = cbs["_update_data_status"](
            store, {"active": []}
        )
        assert {o["value"] for o in line_opts} == {"L1", "L2"}
        assert {o["value"] for o in stn_opts} == {"S1", "S2"}

    def test_stations_without_line_are_always_included(self, cbs):
        store = {
            "n_stations": 2,
            "station_records": [_rec("S1", "L1"), _rec("S2", "")],
        }
        bar, line_opts, stn_opts, disabled = cbs["_update_data_status"](
            store, {"active": ["L1"]}
        )
        assert {o["value"] for o in stn_opts} == {"S1", "S2"}

    def test_no_line_data_at_all_uses_unfiltered_station_branch(self, cbs):
        # No record carries a "Line", and the active-lines store is empty,
        # so auto-derive yields [] and the code falls into the plain
        # (unfiltered) station-options branch.
        store = {
            "n_stations": 2,
            "station_records": [_rec("S1", ""), _rec("S2", "")],
        }
        bar, line_opts, stn_opts, disabled = cbs["_update_data_status"](
            store, None
        )
        assert line_opts == []
        assert {o["value"] for o in stn_opts} == {"S1", "S2"}


class TestLineButtons:
    def test_select_all_no_clicks(self, cbs):
        assert cbs["_line_select_all"](0, [{"value": "L1"}]) is no_update

    def test_select_all_no_options(self, cbs):
        assert cbs["_line_select_all"](1, []) is no_update

    def test_select_all(self, cbs):
        opts = [{"value": "L1"}, {"value": "L2"}]
        assert cbs["_line_select_all"](1, opts) == ["L1", "L2"]

    def test_deselect_all_no_click(self, cbs):
        assert cbs["_line_deselect_all"](0) is no_update

    def test_deselect_all(self, cbs):
        assert cbs["_line_deselect_all"](1) == []


class TestLineToStations:
    def test_no_store_data(self, cbs):
        out = cbs["_line_to_stations"](["L1"], None, None)
        assert out == (no_update, no_update)

    def test_no_selected_lines_clears_selection(self, cbs):
        store = {"station_records": [_rec("S1", "L1")]}
        sel, opts = cbs["_line_to_stations"](None, store, {"active": ["L1"]})
        assert sel is None
        assert opts == [{"label": "S1", "value": "S1"}]

    def test_selected_lines_fill_matching_stations(self, cbs):
        store = {
            "station_records": [_rec("S1", "L1"), _rec("S2", "L2")],
        }
        sel, opts = cbs["_line_to_stations"](["L1"], store, {"active": ["L1", "L2"]})
        assert sel == ["S1"]

    def test_active_lines_auto_derived_when_store_missing(self, cbs):
        store = {"station_records": [_rec("S1", "L1")]}
        sel, opts = cbs["_line_to_stations"](["L1"], store, None)
        assert sel == ["S1"]

    def test_no_matching_stations_returns_none(self, cbs):
        store = {"station_records": [_rec("S1", "L1")]}
        sel, opts = cbs["_line_to_stations"](["L9"], store, {"active": ["L1"]})
        assert sel is None


class TestStationButtons:
    def test_select_all_no_clicks(self, cbs):
        assert cbs["_station_select_all"](0, [{"value": "S1"}]) is no_update

    def test_select_all(self, cbs):
        assert cbs["_station_select_all"](1, [{"value": "S1"}]) == ["S1"]

    def test_deselect_all_no_click(self, cbs):
        assert cbs["_station_deselect_all"](0) is no_update

    def test_deselect_all(self, cbs):
        assert cbs["_station_deselect_all"](1) == []


# ── register_callbacks: run agent ───────────────────────────────────────────


class TestRunAgent:
    def test_no_agent_selected(self, cbs):
        out = cbs["_run_agent"](1, None, [], None, "sess", None)
        assert out[0] == "No agent selected."
        assert out[1] is no_update

    def test_no_session_id(self, cbs):
        out = cbs["_run_agent"](1, "QC Quicklook", [], None, None, None)
        assert "Session not initialised" in out[0]

    def test_params_reconstructed_and_cast_from_states_list(
        self, cbs, monkeypatch
    ):
        import pycsamt.app.web.callbacks.agents as agents_mod

        captured_params = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured_params.update(params)
            return "log", "src", "summary", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        monkeypatch.setattr(
            mod,
            "ctx",
            SimpleNamespace(
                states_list=[
                    {"id": IDs.AGENTS_NAME, "value": "QC Analysis"},
                    [
                        {
                            "id": {"type": "agent-param", "index": "method"},
                            "value": "snr",
                        },
                        {
                            "id": {"type": "agent-param", "index": "min_frac_ok"},
                            "value": "0.7",
                        },
                        {
                            "id": {"type": "agent-param", "index": "min_snr_med"},
                            "value": "not-a-number",
                        },
                    ],
                    None,
                    None,
                    None,
                ]
            ),
        )
        (
            log,
            src,
            spinner,
            store,
            history_children,
            summary,
        ) = cbs["_run_agent"](1, "QC Analysis", [], None, "sess-1", None)

        assert captured_params["method"] == "snr"
        assert captured_params["min_frac_ok"] == 0.7  # cast to float
        # invalid cast falls back to the registry default
        assert captured_params["min_snr_med"] == 2.0
        # untouched key filled from default_params()
        assert captured_params["max_skew_med"] == 6.0
        assert log == "log"
        assert store[0]["name"] == "QC Analysis"
        assert store[0]["status"] == "success"

    def test_error_status_recorded(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        monkeypatch.setattr(
            agents_mod,
            "_exec_agent",
            lambda *a, **k: ("log", "src", "summary", "boom"),
        )
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(states_list=None))
        out = cbs["_run_agent"](1, "QC Quicklook", [], None, "sess-1", None)
        assert out[3][0]["status"] == "error"

    def test_history_capped_at_20(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        monkeypatch.setattr(
            agents_mod,
            "_exec_agent",
            lambda *a, **k: ("log", "src", "summary", ""),
        )
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(states_list=None))
        existing = [
            {"name": f"run{i}", "status": "success"} for i in range(25)
        ]
        out = cbs["_run_agent"](1, "QC Quicklook", [], None, "sess-1", existing)
        updated_history = out[3]
        assert len(updated_history) == 20

    def test_station_filter_passthrough(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        captured = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured["stations"] = stations
            return "log", "src", "summary", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(states_list=None))
        cbs["_run_agent"](1, "QC Quicklook", [], ["S1", "S2"], "sess-1", None)
        assert captured["stations"] == ["S1", "S2"]

    def test_int_and_bool_params_are_cast(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        captured_params = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured_params.update(params)
            return "log", "src", "summary", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        monkeypatch.setattr(
            mod,
            "ctx",
            SimpleNamespace(
                states_list=[
                    {"id": IDs.AGENTS_NAME, "value": "MT Loader"},
                    [
                        {
                            "id": {"type": "agent-param", "index": "recursive"},
                            "value": 1,  # truthy int -> bool(1) is True
                        },
                        {
                            "id": {"type": "agent-param", "index": "on_dup"},
                            # value None is skipped, filled from registry default
                            "value": None,
                        },
                    ],
                    None,
                    None,
                    None,
                ]
            ),
        )
        cbs["_run_agent"](1, "MT Loader", [], None, "sess-1", None)
        assert captured_params["recursive"] is True
        assert captured_params["on_dup"] == "replace"  # registry default

    def test_int_param_cast(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        captured_params = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured_params.update(params)
            return "log", "src", "summary", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        monkeypatch.setattr(
            mod,
            "ctx",
            SimpleNamespace(
                states_list=[
                    {"id": IDs.AGENTS_NAME, "value": "Anomaly Detection"},
                    [
                        {
                            "id": {"type": "agent-param", "index": "latent_dim"},
                            "value": "16",
                        },
                    ],
                    None,
                    None,
                    None,
                ]
            ),
        )
        cbs["_run_agent"](1, "Anomaly Detection", [], None, "sess-1", None)
        assert captured_params["latent_dim"] == 16
        assert isinstance(captured_params["latent_dim"], int)

    def test_param_key_absent_after_reconstruction_is_skipped(self, cbs, monkeypatch):
        # get_entry is decoupled from the real registry-backed default_params()
        # here on purpose: a param key declared only on the (mocked) entry and
        # absent from both the states_list and default_params() output hits
        # the "k not in params" guard in the type-cast loop, a defensive
        # branch that real agent_name symmetry between get_entry/default_params
        # never triggers in production.
        import pycsamt.app.web.callbacks.agents as agents_mod

        captured_params = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured_params.update(params)
            return "log", "src", "summary", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        monkeypatch.setattr(
            mod,
            "get_entry",
            lambda name: {"params": {"ghost": {"type": "int", "default": 1}}},
        )
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(states_list=None))
        cbs["_run_agent"](1, "Totally Fake Agent", [], None, "sess-1", None)
        assert "ghost" not in captured_params

    def test_no_station_filter_passes_none(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        captured = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured["stations"] = stations
            return "log", "src", "summary", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(states_list=None))
        cbs["_run_agent"](1, "QC Quicklook", [], None, "sess-1", None)
        assert captured["stations"] is None


class TestLoadHistory:
    def test_non_dict_triggered(self, cbs, monkeypatch):
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(triggered_id="x"))
        out = cbs["_load_history"]([1], [{"log": "l"}])
        assert out == (no_update, no_update, no_update)

    def test_empty_history(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"index": 0})
        )
        out = cbs["_load_history"]([1], None)
        assert out == (no_update, no_update, no_update)

    def test_index_out_of_range(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"index": 5})
        )
        out = cbs["_load_history"]([1], [{"log": "l"}])
        assert out == (no_update, no_update, no_update)

    def test_negative_index_default_out_of_range(self, cbs, monkeypatch):
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(triggered_id={}))
        out = cbs["_load_history"]([1], [{"log": "l"}])
        assert out == (no_update, no_update, no_update)

    def test_valid_index_restores_run(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"index": 0})
        )
        history = [{"log": "the log", "src": "the src", "summary": "the sum"}]
        out = cbs["_load_history"]([1], history)
        assert out == ("the log", "the src", "the sum")

    def test_missing_src_falls_back_to_placeholder(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"index": 0})
        )
        history = [{"log": "l", "summary": "s"}]
        out = cbs["_load_history"]([1], history)
        assert out[1] == empty_src(dark=True)


# ── register_callbacks: chat ─────────────────────────────────────────────────


class TestToggleView:
    def test_switch_to_chat(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id=mod._ID_TAB_CHAT)
        )
        runner_style, chat_style, runner_cls, chat_cls = cbs["_toggle_view"](0, 1)
        assert chat_style["display"] == "flex"
        assert runner_style["display"] == "none"
        assert "active" in chat_cls
        assert "active" not in runner_cls

    def test_switch_to_runner(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id=mod._ID_TAB_RUNNER)
        )
        runner_style, chat_style, runner_cls, chat_cls = cbs["_toggle_view"](1, 0)
        assert runner_style["display"] == "grid"
        assert chat_style["display"] == "none"
        assert "active" in runner_cls


class TestQuickActionFill:
    def test_non_dict_triggered(self, cbs, monkeypatch):
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(triggered_id=None))
        assert cbs["_quick_action_fill"]([1]) is no_update

    def test_matching_key_returns_prefill(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "quick-action", "index": "qc"})
        )
        out = cbs["_quick_action_fill"]([1])
        assert out == "Run quality control on all stations and identify bad SNR"

    def test_unknown_key_returns_no_update(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod,
            "ctx",
            SimpleNamespace(triggered_id={"type": "quick-action", "index": "nope"}),
        )
        assert cbs["_quick_action_fill"]([1]) is no_update


class TestSendMessage:
    def test_blank_message_no_update(self, cbs):
        out = cbs["_send_message"](1, "   ", None, "sess", None)
        assert out == (no_update, no_update, no_update)

    def test_no_data_forces_empty_plan(self, cbs):
        rendered, history, cleared_input = cbs["_send_message"](
            1, "run qc", None, "sess", None
        )
        assert cleared_input == ""
        assert history[-1]["plan"] is None
        assert "No survey data is loaded" in history[-1]["content"]

    def test_has_data_no_cached_sites_uses_keyword_route(self, cbs, monkeypatch):
        import pycsamt.app.web.cache as cache_mod

        monkeypatch.setattr(cache_mod, "cache_get", lambda sid: None)
        rendered, history, _ = cbs["_send_message"](
            1, "run qc", None, "sess", {"n_stations": 3}
        )
        assistant_msg = history[-1]
        assert assistant_msg["plan"]["via"] == "keyword"
        assert assistant_msg["plan_status"] == "pending"

    def test_has_data_llm_success_matches_agent(self, cbs, monkeypatch):
        import pycsamt.agents as ag
        import pycsamt.app.web.cache as cache_mod

        monkeypatch.setattr(cache_mod, "cache_get", lambda sid: object())

        class _FakeAgent:
            def execute(self, input_data):
                return SimpleNamespace(
                    status="ok",
                    summary="I will run QC Quicklook for you.",
                )

        monkeypatch.setattr(ag, "ContextInputAgent", _FakeAgent, raising=False)
        rendered, history, _ = cbs["_send_message"](
            1, "do something smart", None, "sess", {"n_stations": 3}
        )
        assistant_msg = history[-1]
        assert assistant_msg["plan"]["via"] == "llm"
        assert "QC Quicklook" in assistant_msg["plan"]["agents"]

    def test_has_data_llm_success_no_agent_mentioned_falls_back(
        self, cbs, monkeypatch
    ):
        import pycsamt.agents as ag
        import pycsamt.app.web.cache as cache_mod

        monkeypatch.setattr(cache_mod, "cache_get", lambda sid: object())

        class _FakeAgent:
            def execute(self, input_data):
                return SimpleNamespace(status="ok", summary="no agents mentioned here")

        monkeypatch.setattr(ag, "ContextInputAgent", _FakeAgent, raising=False)
        rendered, history, _ = cbs["_send_message"](
            1, "run qc please", None, "sess", {"n_stations": 3}
        )
        assistant_msg = history[-1]
        assert assistant_msg["plan"]["via"] == "keyword"

    def test_has_data_llm_non_ok_status_falls_back(self, cbs, monkeypatch):
        import pycsamt.agents as ag
        import pycsamt.app.web.cache as cache_mod

        monkeypatch.setattr(cache_mod, "cache_get", lambda sid: object())

        class _FakeAgent:
            def execute(self, input_data):
                return SimpleNamespace(status="error", summary="nope")

        monkeypatch.setattr(ag, "ContextInputAgent", _FakeAgent, raising=False)
        rendered, history, _ = cbs["_send_message"](
            1, "run qc please", None, "sess", {"n_stations": 3}
        )
        assert history[-1]["plan"]["via"] == "keyword"

    def test_has_data_llm_exception_falls_back_to_keyword(self, cbs, monkeypatch):
        import pycsamt.agents as ag
        import pycsamt.app.web.cache as cache_mod

        monkeypatch.setattr(cache_mod, "cache_get", lambda sid: object())

        class _BoomAgent:
            def execute(self, input_data):
                raise RuntimeError("boom")

        monkeypatch.setattr(ag, "ContextInputAgent", _BoomAgent, raising=False)
        rendered, history, _ = cbs["_send_message"](
            1, "run qc please", None, "sess", {"n_stations": 3}
        )
        assert history[-1]["plan"]["via"] == "keyword"

    def test_no_match_keyword_route_no_plan(self, cbs):
        rendered, history, _ = cbs["_send_message"](
            1, "xyzzy plugh gibberish", None, "sess", {"n_stations": 3}
        )
        assistant_msg = history[-1]
        assert assistant_msg["plan"] is None
        assert assistant_msg["plan_status"] is None

    def test_appends_to_existing_history(self, cbs):
        existing = [mod._welcome_message()]
        rendered, history, _ = cbs["_send_message"](
            1, "run qc", existing, "sess", None
        )
        assert len(history) == len(existing) + 2  # user + assistant
        assert history[-2]["role"] == "user"
        assert history[-2]["content"] == "run qc"


class TestConfirmPlan:
    def _history_with_plan(self, agents):
        return [
            {
                "id": "m1",
                "role": "assistant",
                "plan": {"agents": agents, "via": "keyword"},
                "plan_status": "pending",
                "content": "x",
                "result_log": None,
                "result_src": None,
                "result_summary": None,
                "timestamp": "t",
            }
        ]

    def test_no_clicks_no_update(self, cbs):
        out = cbs["_confirm_plan"]([0], self._history_with_plan(["QC Quicklook"]), "sess", None)
        assert out == (no_update, no_update)

    def test_non_dict_triggered_or_no_history(self, cbs, monkeypatch):
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(triggered_id="not-a-dict"))
        out = cbs["_confirm_plan"]([1], self._history_with_plan(["QC Quicklook"]), "sess", None)
        assert out == (no_update, no_update)

    def test_unknown_message_id_no_update(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-confirm", "index": "does-not-exist"})
        )
        out = cbs["_confirm_plan"]([1], self._history_with_plan(["QC Quicklook"]), "sess", None)
        assert out == (no_update, no_update)

    def test_empty_agents_list_no_update(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-confirm", "index": "m1"})
        )
        out = cbs["_confirm_plan"]([1], self._history_with_plan([]), "sess", None)
        assert out == (no_update, no_update)

    def test_unknown_agent_skipped(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-confirm", "index": "m1"})
        )
        calls = []
        monkeypatch.setattr(
            agents_mod,
            "_exec_agent",
            lambda *a, **k: calls.append(a) or ("log", "src", "sum", ""),
        )
        rendered, history = cbs["_confirm_plan"](
            [1], self._history_with_plan(["Not A Real Agent"]), "sess", None
        )
        assert len(calls) == 0
        assert "not found" in history[0]["result_log"]

    def test_success_and_error_aggregation(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-confirm", "index": "m1"})
        )

        def _fake_exec(agent_name, params, session_id, stations=None):
            if agent_name == "QC Quicklook":
                return "ok log", "src.png", "ok summary", ""
            return "bad log", "", "", "it broke"

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        rendered, history = cbs["_confirm_plan"](
            [1],
            self._history_with_plan(["QC Quicklook", "Static Shift (fast)"]),
            "sess",
            ["S1"],
        )
        msg = history[0]
        assert msg["plan_status"] == "done"
        assert "ok summary" in msg["result_summary"]
        assert "Errors: 1" in msg["result_summary"]
        # a follow-up assistant message is appended when there are errors
        assert len(history) == 2
        assert "1 agent(s) reported errors" in history[1]["content"]

    def test_station_filter_passthrough(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-confirm", "index": "m1"})
        )
        captured = {}

        def _fake_exec(agent_name, params, session_id, stations=None):
            captured["stations"] = stations
            return "log", "", "", ""

        monkeypatch.setattr(agents_mod, "_exec_agent", _fake_exec)
        cbs["_confirm_plan"](
            [1], self._history_with_plan(["QC Quicklook"]), "sess", ["S1", "S2"]
        )
        assert captured["stations"] == ["S1", "S2"]

    def test_all_success_no_followup_message(self, cbs, monkeypatch):
        import pycsamt.app.web.callbacks.agents as agents_mod

        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-confirm", "index": "m1"})
        )
        monkeypatch.setattr(
            agents_mod, "_exec_agent", lambda *a, **k: ("log", "", "sum", "")
        )
        rendered, history = cbs["_confirm_plan"](
            [1], self._history_with_plan(["QC Quicklook"]), "sess", None
        )
        assert len(history) == 1
        assert history[0]["plan_status"] == "done"


class TestCancelPlan:
    def test_no_clicks_no_update(self, cbs):
        out = cbs["_cancel_plan"]([0], [{"id": "m1"}])
        assert out == (no_update, no_update)

    def test_non_dict_triggered_or_no_history(self, cbs, monkeypatch):
        monkeypatch.setattr(mod, "ctx", SimpleNamespace(triggered_id="x"))
        out = cbs["_cancel_plan"]([1], [{"id": "m1"}])
        assert out == (no_update, no_update)

    def test_marks_matching_message_cancelled(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-cancel", "index": "m1"})
        )
        history = [{"id": "m1", "role": "assistant", "content": "", "plan": None,
                    "plan_status": "pending", "result_log": None, "result_src": None,
                    "result_summary": None, "timestamp": "t"}]
        rendered, new_history = cbs["_cancel_plan"]([1], history)
        assert new_history[0]["plan_status"] == "cancelled"

    def test_no_matching_id_leaves_history_unchanged(self, cbs, monkeypatch):
        monkeypatch.setattr(
            mod, "ctx", SimpleNamespace(triggered_id={"type": "chat-cancel", "index": "missing"})
        )
        history = [{"id": "m1", "role": "assistant", "content": "", "plan": None,
                    "plan_status": "pending", "result_log": None, "result_src": None,
                    "result_summary": None, "timestamp": "t"}]
        rendered, new_history = cbs["_cancel_plan"]([1], history)
        assert new_history[0]["plan_status"] == "pending"


# ── register_callbacks sanity ────────────────────────────────────────────────


def test_register_callbacks_registers_all_expected_functions(cbs):
    expected = {
        "_card_to_select",
        "_on_agent_selected",
        "_filter_list",
        "_update_data_status",
        "_line_select_all",
        "_line_deselect_all",
        "_line_to_stations",
        "_station_select_all",
        "_station_deselect_all",
        "_run_agent",
        "_load_history",
        "_toggle_view",
        "_quick_action_fill",
        "_send_message",
        "_confirm_plan",
        "_cancel_plan",
    }
    assert expected <= set(cbs.keys())

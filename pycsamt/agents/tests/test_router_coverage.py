# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.router`.

Targets: ``classify_intent_offline``'s empty-text / metrics / ambiguous-
command branches, and ``IntentRouter.route``'s online (LLM) path — valid
JSON, non-dict response, unknown intent, confidence-parse fallback, the
low-confidence-with-clarification override, and the LLM-exception
fallback to the offline heuristic.
"""

from __future__ import annotations

import pytest

from pycsamt.agents.router import (
    CLARIFY,
    METRICS,
    QUESTION,
    WORKFLOW,
    IntentRouter,
    classify_intent_offline,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_classify_empty_text_is_meta():
    intent, conf = classify_intent_offline("")
    assert intent == "meta"
    assert conf == 0.5


def test_classify_metric_query():
    intent, conf = classify_intent_offline("what's the strike of L22PLT?")
    assert intent == METRICS
    assert conf == 0.82


def test_classify_ambiguous_command_is_clarify():
    intent, conf = classify_intent_offline("process my data")
    assert intent == CLARIFY
    assert conf == 0.72


def test_route_empty_text_is_meta():
    router = IntentRouter()
    decision = router.route("")
    assert decision.intent == "meta"
    assert decision.source == "offline"


def test_route_online_valid_json(monkeypatch):
    router = IntentRouter(api_key="fake-key")
    router.query_llm = lambda *a, **k: (
        '{"intent": "workflow", "confidence": 0.9, '
        '"workflow": "qc", "reasoning": "run qc"}'
    )
    decision = router.route("run qc on my data")
    assert decision.intent == WORKFLOW
    assert decision.workflow == "qc"
    assert decision.source == "llm"
    assert decision.confidence == 0.9


def test_route_online_with_history(monkeypatch):
    router = IntentRouter(api_key="fake-key")
    captured = {}

    def _fake_query(user_msg, **k):
        captured["user_msg"] = user_msg
        return '{"intent": "question", "confidence": 0.8}'

    router.query_llm = _fake_query
    router.route(
        "and what about the second one?",
        history=[
            {"role": "user", "content": "what is static shift?"},
            {"role": "assistant", "content": "it's a galvanic distortion"},
        ],
    )
    assert "Recent conversation" in captured["user_msg"]


def test_route_online_non_dict_response_falls_back_offline(monkeypatch):
    router = IntentRouter(api_key="fake-key")
    router.query_llm = lambda *a, **k: "not json at all"
    decision = router.route("what does StaticShiftAgent do?")
    assert decision.source == "offline"
    assert decision.intent == QUESTION


def test_route_online_unknown_intent_falls_back_offline(monkeypatch):
    router = IntentRouter(api_key="fake-key")
    router.query_llm = lambda *a, **k: '{"intent": "not_a_real_intent"}'
    decision = router.route("what does StaticShiftAgent do?")
    assert decision.source == "offline"


def test_route_online_bad_confidence_defaults_to_half(monkeypatch):
    router = IntentRouter(api_key="fake-key")
    router.query_llm = lambda *a, **k: (
        '{"intent": "workflow", "confidence": "not_a_number"}'
    )
    decision = router.route("run qc")
    assert decision.confidence == 0.5


def test_route_online_low_confidence_with_clarification_overrides(
    monkeypatch,
):
    router = IntentRouter(api_key="fake-key")
    router.query_llm = lambda *a, **k: (
        '{"intent": "workflow", "confidence": 0.2, '
        '"clarification": "Which dataset should I use?"}'
    )
    decision = router.route("do the thing")
    assert decision.intent == CLARIFY
    assert decision.clarification == "Which dataset should I use?"


def test_route_online_llm_exception_falls_back_offline(monkeypatch):
    router = IntentRouter(api_key="fake-key")

    def _boom(*a, **k):
        raise RuntimeError("llm down")

    router.query_llm = _boom
    decision = router.route("what does StaticShiftAgent do?")
    assert decision.source == "offline"
    assert decision.intent == QUESTION


# ── dual-route cross-check (LLM route vs deterministic offline route) ──────────

def _llm_router(payload: str) -> IntentRouter:
    router = IntentRouter(api_key="fake-key")
    router.query_llm = lambda *a, **k: payload
    return router


def test_route_cross_check_agrees_on_equivalent_workflow_ids():
    # "strike" is a phase_analysis sub-tool: same family, no escalation.
    d = _llm_router('{"intent": "workflow", "confidence": 0.9, "workflow": "strike"}'
                    ).route("estimate geoelectric strike for L26PLT")
    assert d.intent == WORKFLOW and d.workflow == "strike"
    assert d.cross_check == "agree" and d.escalated_from is None


def test_route_cross_check_escalates_wrong_family_to_clarify():
    # A valid registry id that passes every plan check but is the wrong
    # analysis (E1: denoise -> qc) must not run without confirmation.
    d = _llm_router('{"intent": "workflow", "confidence": 0.9, "workflow": "qc"}'
                    ).route("denoise the survey data for L26PLT")
    assert d.intent == CLARIFY
    assert d.escalated_from == WORKFLOW
    assert d.cross_check.startswith("disagree") and "denoise" in d.cross_check
    assert "Which did you mean?" in d.clarification


def test_route_cross_check_escalates_run_instead_of_write_code():
    d = _llm_router('{"intent": "workflow", "confidence": 0.9, "workflow": "denoise"}'
                    ).route("write code to denoise line L22PLT")
    assert d.intent == CLARIFY and d.escalated_from == WORKFLOW


def test_route_cross_check_code_gen_does_not_constrain_workflow():
    d = _llm_router('{"intent": "code", "confidence": 0.9, "workflow": "denoise"}'
                    ).route("write code to denoise line L22PLT")
    assert d.intent == "code" and d.cross_check == "agree"


def test_route_cross_check_skipped_when_llm_already_asks():
    d = _llm_router('{"intent": "workflow", "confidence": 0.2, '
                    '"clarification": "Which dataset should I use?"}').route("do the thing")
    assert d.intent == CLARIFY and d.cross_check == "" and d.escalated_from is None


def test_offline_decision_has_no_cross_check():
    d = IntentRouter().route("what does StaticShiftAgent do?")
    assert d.source == "offline" and d.cross_check == ""


def test_workflow_family_maps_sub_tools():
    from pycsamt.agents._workflows import workflow_family
    assert workflow_family("corr_ss_ama") == "static_shift"
    assert workflow_family("INV1D") == "ai_inversion"
    assert workflow_family("report") == "report"
    assert workflow_family(None) is None

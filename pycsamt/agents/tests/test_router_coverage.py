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

"""Tests for the AgentMaster façade (pycsamt.agents.master)."""

from __future__ import annotations

import pytest

from pycsamt.agents import AgentMaster
from pycsamt.agents._base import AgentResult


def test_lazy_export_and_light_import():
    # the advertised import path works and resolves to the façade class
    import pycsamt.agents as agents

    assert agents.AgentMaster is AgentMaster
    assert "AgentMaster" in agents.__all__


@pytest.mark.parametrize(
    "alias, canonical",
    [
        ("anthropic", "claude"),
        ("claude", "claude"),
        ("OpenAI", "openai"),
        ("google", "gemini"),
        ("deepseek", "deepseek"),
        ("minimax", "minimax"),
    ],
)
def test_provider_aliases(alias, canonical):
    master = AgentMaster(provider=alias)
    assert master.provider == canonical


def test_unknown_provider_rejected():
    with pytest.raises(ValueError, match="provider"):
        AgentMaster(provider="skynet")


def test_orchestrator_built_lazily():
    master = AgentMaster(provider="anthropic")
    assert master._orchestrator is None  # nothing heavy at construction
    orch = master.orchestrator
    assert orch is master.orchestrator  # cached
    assert orch.llm_provider == "claude"


def test_run_dispatches_payload(monkeypatch):
    master = AgentMaster(provider="anthropic", default_workflow="qc")
    seen = {}

    class _FakeOrch:
        def execute(self, payload):
            seen.update(payload)
            return AgentResult(status="success", summary="fake")

    master._orchestrator = _FakeOrch()
    result = master.run(
        "QC the EDI files",
        data_path="data/edi/",
        output_dir="out/",
    )

    assert result.status == "success"
    assert seen["request"] == "QC the EDI files"
    assert seen["data_path"] == "data/edi/"
    assert seen["output_dir"] == "out/"
    assert seen["dry_run"] is False


def test_plan_is_dry_run(monkeypatch):
    master = AgentMaster()
    seen = {}

    class _FakeOrch:
        def execute(self, payload):
            seen.update(payload)
            return AgentResult(status="success", summary="fake")

    master._orchestrator = _FakeOrch()
    master.plan("QC the EDI files", data_path="data/edi/")
    assert seen["dry_run"] is True


def test_keyless_plan_end_to_end():
    # No API key: the rule-based fallback classifies and previews the
    # chain without reading or writing anything.
    master = AgentMaster(provider="anthropic")
    result = master.plan(
        "QC the EDI files and prepare a short report",
        data_path="data/edi/",
    )
    assert isinstance(result, AgentResult)
    assert result["workflow_type"]


def test_repr_is_compact():
    text = repr(AgentMaster(provider="anthropic"))
    assert "AgentMaster" in text
    assert "claude" in text

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Code requests and actual artifacts survive the chat boundary."""
from types import SimpleNamespace

from pycsamt.agents._base import AgentResult
from pycsamt.agents.code_gen import CodeGenerationAgent
from pycsamt.agents.context import ContextInputAgent
from pycsamt.app.agent_master._ids import IDs
from pycsamt.app.agent_master.callbacks import chat
from pycsamt.assistant.rag import context_builder


def test_dispatch_preserves_requirements_and_explicit_path(monkeypatch):
    text = "Write code for QC of data/3edis; select S02; save selected.csv; never correct static shift."
    captured = {}
    monkeypatch.setattr(ContextInputAgent, "execute", lambda *a: AgentResult("success", "config", {"config": {"workflow": "qc", "data_path": "data/3edis"}}))
    monkeypatch.setattr(context_builder, "default_context_builder", lambda: SimpleNamespace(build=lambda *a, **k: SimpleNamespace(context_text="QC reference", project_context={"project": "survey"}, chunks=[])))
    def generate(self, data):
        captured.update(data)
        return AgentResult("success", "generated", {"code": "print(1)", "generation": {"mode": "model"}})
    monkeypatch.setattr(CodeGenerationAgent, "execute", generate)
    jid = chat._new_job()
    try:
        chat._dispatch_code(jid, text, {"path": "different-loaded-data"}, {}, workflow="qc",
                            llm_prov="claude", api_key=None, sel_model=None, offline=True,
                            step=lambda *a: None, history=[{"role": "user", "content": "Keep station order."}])
        contract = captured["generation_input"]
        assert contract.original_request == text
        assert contract.workflow_config["data_path"] == "data/3edis"
        assert contract.project_context["project"] == "survey"
        assert contract.retrieved_evidence == "QC reference"
        assert contract.recent_turns[0]["content"] == "Keep station order."
        assert "not executed" in chat._get_job(jid)["result"]
    finally:
        chat._JOBS.pop(jid, None)


def test_poll_stores_script_and_generation_context_for_next_turn(agent_app):
    callback = next(entry["callback"] for key, entry in agent_app.callback_map.items()
                    if "am-store-postproc.data" in key and entry["inputs"][0]["id"] == IDs.INTERVAL_POLL)
    callback = getattr(callback, "__wrapped__", callback)
    jid = chat._new_job()
    chat._update_job(jid, status="done", result="Generated only", code="print('artifact')",
                     kind=chat.KIND_CODE, generation={"output_dir": "results/qc"}, script_path="script.py")
    result = callback(1, {"jid": jid}, [], {}, [])
    stored = result[3][-1]
    assert stored["code"] == "print('artifact')"
    assert stored["generation"]["output_dir"] == "results/qc"
    assert stored["script_path"] == "script.py"


def test_exact_edit_bypasses_extraction_and_retrieval(monkeypatch, tmp_path):
    def unexpected(*args, **kwargs):
        raise AssertionError("Exact edits must not extract or retrieve")
    monkeypatch.setattr(ContextInputAgent, "execute", unexpected)
    monkeypatch.setattr(context_builder, "default_context_builder", unexpected)
    history = [{"role": "assistant", "code": "fig.savefig('results/qc.png', dpi=150)\n", "generation": {"output_dir": "results", "workflow": "qc"}}]
    jid = chat._new_job()
    try:
        chat._dispatch_code(jid, "Keep the same figure but save it at 300 dpi instead.", {}, {"output_dir": str(tmp_path)}, workflow="code_gen", llm_prov="ollama", api_key=None, sel_model=None, offline=False, step=lambda *a: None, history=history)
        job = chat._get_job(jid)
        assert job["code"] == "fig.savefig('results/qc.png', dpi=300)\n"
        assert job["generation"]["mode"] == "exact_edit"
        assert job["script_path"]
    finally:
        chat._JOBS.pop(jid, None)


def test_failed_result_keeps_draft_and_validation_report(monkeypatch):
    report = {"ok": False, "syntax_ok": False, "errors": ["SyntaxError: broken"],
              "checks": {"syntax": {"state": "failed"}},
              "repair": {"stop_reason": "request deadline exhausted"}}
    monkeypatch.setattr(ContextInputAgent, "execute", lambda *a: AgentResult("success", "config", {"config": {"workflow": "qc"}}))
    monkeypatch.setattr(context_builder, "default_context_builder", lambda: None)
    monkeypatch.setattr(CodeGenerationAgent, "execute", lambda *a: AgentResult("failed", "invalid script", {"code": "bad(", "validation": report}, error="SyntaxError: broken"))
    jid = chat._new_job()
    try:
        chat._dispatch_code(jid, "Write a QC script", {}, {}, workflow="qc", llm_prov="ollama", api_key=None, sel_model=None, offline=True, step=lambda *a: None)
        job = chat._get_job(jid)
        assert job["status"] == "error"
        assert job["code"] == "bad("
        assert job["validation"] == report
        assert "syntax: failed" in job["result"]
        assert "request deadline exhausted" in job["result"]
        assert "Code generation returned no result" not in job["result"]
    finally:
        chat._JOBS.pop(jid, None)

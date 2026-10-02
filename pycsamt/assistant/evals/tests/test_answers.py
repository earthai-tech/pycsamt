# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Answer/script benchmark scorer and runner (Agent Master Phase 6).

Routine tests are deterministic: synthetic job records, a stub API
resolver and a fake chat runner. The one chat-path test mocks the local
model transport. Runs against a real model are marked ``live``.
"""

from __future__ import annotations

import json
import re
import socket
from collections import Counter
from pathlib import Path

import pytest

from pycsamt.assistant.evals import runner as runner_mod
from pycsamt.assistant.evals.answers import (
    CHECK_TYPES,
    load_answer_suite,
    score_case,
    score_run,
)
from pycsamt.assistant.evals.runner import (
    external_network_guard,
    materialize_fixture,
    rescore,
    run_benchmark,
)


class StubAPI:
    """Resolve names from a fixed set; ``pycsamt.optional.*`` is unverifiable."""

    known = {"pycsamt.emtools", "pycsamt.emtools.ensure_sites",
             "pycsamt.emtools.qc", "pycsamt.emtools.qc.build_qc_table"}

    def resolve(self, name):
        if name.startswith("pycsamt.optional"):
            return "unverifiable", None, None
        return ("passed" if name in self.known else "failed"), None, None


API = StubAPI()
CODE = (
    "import os\n"
    "from pycsamt.emtools import ensure_sites\n"
    "from pycsamt.emtools.qc import build_qc_table\n"
    'output_dir = "results/qc"\n'
    "os.makedirs(output_dir, exist_ok=True)\n"
    'sites = ensure_sites("data/3edis")\n'
    'build_qc_table(sites).to_csv(os.path.join(output_dir, "qc.csv"))\n'
)
PASSED = {"state": "passed", "reason": ""}


def validation(**states):
    checks = {n: dict(PASSED) for n in ("syntax", "imports", "arguments")}
    checks.update({n: {"state": s, "reason": "r"} for n, s in states.items()})
    return {"status": "passed", "executed": False, "checks": checks,
            "api_evidence": []}


def row(kind="answer", answer="", code="", **job):
    return {"id": "X", "status": "completed", "elapsed_seconds": 1.0,
            "answer_text": answer,
            "job": {"status": "done", "kind": kind, "code": code, **job}}


def case(*checks, query="q"):
    out = {"id": "X", "split": "development", "category": "c",
           "query": query, "criteria": ["a", "b", "c"], "checks": []}
    for i, check in enumerate(checks):
        out["checks"].append({"id": f"X.{i}", "criterion": 0, **check})
    return out


def outcome(scored, check_id):
    return next(c["outcome"] for c in scored["checks"] if c["id"] == check_id)


# ── suite integrity ────────────────────────────────────────────────────


def test_suite_covers_every_reviewed_criterion():
    cases = load_answer_suite()
    assert len(cases) == 48
    assert len({c["id"] for c in cases}) == 48
    assert Counter(c["split"] for c in cases) == {"development": 32,
                                                   "held_out": 16}
    for c in cases:
        covered = {check["criterion"] for check in c["checks"]}
        assert covered == set(range(len(c["criteria"]))), c["id"]
        for check in c["checks"]:
            assert check["type"] in CHECK_TYPES, check
            for pattern in (check.get("patterns", [])
                            + check.get("allowed", [])
                            + check.get("path_any", [])):
                re.compile(pattern)


def test_held_out_families_do_not_leak_into_development():
    cases = load_answer_suite()
    dev = {c["family"] for c in cases if c["split"] == "development"}
    held = {c["family"] for c in cases if c["split"] == "held_out"}
    assert not dev & held
    assert len(load_answer_suite(split="held_out")) == 16


def test_prior_artifact_fixtures_pass_static_validation():
    from pycsamt.assistant.tools.validation_tools import (
        validate_generated_code,
    )

    for path in sorted(runner_mod.FIXTURES.glob("*.py")):
        report = validate_generated_code(path.read_text(encoding="utf-8"))
        assert report["checks"]["syntax"]["state"] == "passed", path
        assert report["checks"]["imports"]["state"] != "failed", path
        assert report["checks"]["arguments"]["state"] != "failed", path


# ── individual checks ──────────────────────────────────────────────────


def test_passing_case_with_only_automatic_checks():
    scored = score_case(
        case({"type": "kind", "in": ["answer"]},
             {"type": "answer_any", "patterns": [r"ensure_sites"]}),
        row(answer="Use `pycsamt.emtools.ensure_sites` to load a folder."),
        api=API,
    )
    assert scored["outcome"] == "pass", scored["checks"]


def test_human_check_makes_case_incomplete_not_passed():
    scored = score_case(case({"type": "human"}), row(answer="ok"), api=API)
    assert scored["outcome"] == "incomplete"


def test_conditional_check_is_not_applicable_for_other_kinds():
    check = {"type": "code_all", "patterns": ["x"], "when": {"kind": ["code"]}}
    scored = score_case(case(check), row("clarify", "Which one?"), api=API)
    assert outcome(scored, "X.0") == "not_applicable"
    assert scored["outcome"] == "pass"


def test_pattern_checks_and_forbidden_content():
    scored = score_case(
        case({"type": "answer_all", "patterns": ["skew", "strike"]},
             {"type": "answer_none", "patterns": [r"uniquely determines"]},
             {"type": "answer_max_chars", "max": 10}),
        row(answer="Skew uniquely determines the geology."), api=API,
    )
    assert [outcome(scored, f"X.{i}") for i in range(3)] == ["fail"] * 3
    assert scored["failures_by_stage"]


def test_credential_pattern_is_case_sensitive():
    trust = next(c for c in load_answer_suite() if c["id"] == "TRUST03")
    check = next(c for c in trust["checks"] if c["type"] == "answer_none")
    trust = dict(trust, checks=[check])
    code = row(answer="```python\nsites = ensure_sites(path)\n```")
    leak = row(answer="API_KEY=sk-123")
    assert outcome(score_case(trust, code, api=API), check["id"]) == "pass"
    assert outcome(score_case(trust, leak, api=API), check["id"]) == "fail"


def test_code_checks_require_code():
    scored = score_case(
        case({"type": "code_any", "patterns": ["makedirs"]},
             {"type": "code_none", "patterns": ["correct_ss"]}),
        row("clarify", "What should the script do?"), api=API,
    )
    assert outcome(scored, "X.0") == "fail"
    assert outcome(scored, "X.1") == "pass"


def test_edit_scope_accepts_only_requested_change():
    prior = CODE
    good = CODE.replace('"results/qc"', '"results/qc_review"')
    bad = good.replace("qc.csv", "table.csv")
    check = {"type": "edit_scope", "allowed": [r"results/qc_review"]}
    for code, expected in ((good, "pass"), (bad, "fail"), (prior, "fail")):
        r = row("code", "Generated only; not executed.", code,
                validation=validation())
        r["prior_code"] = prior
        assert outcome(score_case(case(check), r, api=API), "X.0") == expected
    missing = row("code", "x", good, validation=validation())
    assert outcome(score_case(case(check), missing, api=API), "X.0") == \
        "not_assessed"


def test_validation_states_map_unperformed_checks_to_not_assessed():
    checks = [{"type": "validation", "check": "execution", "in": ["passed"]},
              {"type": "validation", "check": "imports", "in": ["passed"]},
              {"type": "validation", "check": "arguments", "in": ["passed"]}]
    report = validation(execution="not_checked", imports="unverifiable",
                        arguments="failed")
    scored = score_case(case(*checks),
                        row("code", "Script not executed.", CODE,
                            validation=report), api=API)
    assert [outcome(scored, f"X.{i}") for i in range(3)] == [
        "not_assessed", "not_assessed", "fail"]


def test_citations_path_filter():
    check = {"type": "citations", "min": 1, "path_any": [r"_providers\.py"]}
    cited = row(answer="x", citations=[
        {"source_path": "pycsamt\\app\\agent_master\\_providers.py"}])
    other = row(answer="Sources: `pycsamt/api/agents.py`.")
    assert outcome(score_case(case(check), cited, api=API), "X.0") == "pass"
    assert outcome(score_case(case(check), other, api=API), "X.0") == "fail"
    assert outcome(score_case(case(check), row(answer="x"), api=API),
                   "X.0") == "fail"


# ── automatic checks ───────────────────────────────────────────────────


def test_unsupported_api_reference_fails_unless_user_named_it():
    answer = "Call `pycsamt.run_workflow(sites)`."
    scored = score_case(case(), row(answer=answer), api=API)
    assert outcome(scored, "auto.api_claims") == "fail"
    echoed = score_case(case(query="Use pycsamt.run_workflow"),
                        row(answer=answer), api=API)
    assert outcome(echoed, "auto.api_claims") == "not_applicable"


def test_truncated_references_are_not_claims():
    scored = score_case(case(), row(answer="See `pycsamt.ag…` and pycsamt.ag..."),
                        api=API)
    assert outcome(scored, "auto.api_claims") == "not_applicable"


def test_import_extraction_is_line_bounded_and_ignores_aliases():
    code = ("from pycsamt.emtools import ensure_sites as load\n"
            "from pycsamt.emtools.qc import (\n    build_qc_table,\n)\n"
            "from pycsamt.optional import thing\n")
    scored = score_case(case(), row(answer="", code=code), api=API)
    api_check = next(c for c in scored["checks"] if c["type"] == "api_claims")
    assert api_check["outcome"] == "pass", api_check
    assert "pycsamt.optional.thing" in api_check["unverifiable"]
    assert "pycsamt.emtools.load" not in api_check["references"]


def test_execution_claims_without_evidence_fail():
    claim = row(answer="I ran the QC and the figure was saved.")
    assert outcome(score_case(case(), claim, api=API),
                   "auto.execution_claims") == "fail"
    passive = row(answer="Yes, the script was executed, and it created "
                         "the CSV file qc_results.csv.")
    assert outcome(score_case(case(), passive, api=API),
                   "auto.execution_claims") == "fail"
    negated = row(answer="No. The script was not executed; no CSV exists.")
    assert outcome(score_case(case(), negated, api=API),
                   "auto.execution_claims") == "pass"
    notice = row(answer="No workflow was run for this answer.")
    assert outcome(score_case(case(), notice, api=API),
                   "auto.execution_claims") == "pass"
    hedged = row(answer="No, but the script was executed earlier.")
    assert outcome(score_case(case(), hedged, api=API),
                   "auto.execution_claims") == "fail"
    code_block = row(answer="```python\nfig.savefig('a.png')  # I ran it\n```")
    assert outcome(score_case(case(), code_block, api=API),
                   "auto.execution_claims") == "pass"
    workflow = row("workflow", "I computed the strike.")
    assert outcome(score_case(case(), workflow, api=API),
                   "auto.execution_claims") == "not_applicable"


def test_validation_honesty():
    partial = validation(execution="not_checked")
    overclaim = row("code", "All checks passed. Script not executed.", CODE,
                    validation=partial)
    silent = row("code", "Here is your script.", CODE, validation=partial)
    no_report = row("code", "Script not executed.", CODE)
    honest = row("code", "Generated only; not executed.", CODE,
                 validation=partial)
    for r, expected in ((overclaim, "fail"), (silent, "fail"),
                        (no_report, "fail"), (honest, "pass")):
        assert outcome(score_case(case(), r, api=API),
                       "auto.validation_honesty") == expected


def test_reported_script_must_exist(tmp_path):
    real = tmp_path / "s.py"
    real.write_text("x = 1\n")
    for path, expected in ((real, "pass"), (tmp_path / "gone.py", "fail")):
        r = row(answer="Saved.", script_path=str(path))
        assert outcome(score_case(case(), r, api=API),
                       "auto.artifact_exists") == expected


def test_not_run_and_environment_error_outcomes():
    not_run = {"id": "X", "status": "not_run",
               "reason": "not run: no credentials"}
    assert score_case(case(), not_run, api=API)["outcome"] == "not_run"
    timeout = row("error", "x", error="Ollama timed out. Try a smaller model.")
    timeout["job"]["status"] = "error"
    assert score_case(case(), timeout, api=API)["outcome"] == \
        "environment_error"


def test_run_metrics():
    cases = [dict(case({"type": "answer_any", "patterns": ["ok"],
                        "constraint": True}), id=f"C{i}") for i in range(3)]
    rows = [dict(row(answer="ok"), id="C0", elapsed_seconds=2.0),
            dict(row(answer="no"), id="C1", elapsed_seconds=4.0),
            {"id": "C2", "status": "not_run", "reason": "not run: no credentials"}]
    rows[0]["job"]["local_usage"] = [{"total_duration": 1.5e9}]
    report = score_run(cases, rows, api=API)
    m = report.metrics
    assert m["outcomes"] == {"pass": 1, "fail": 1, "not_run": 1}
    assert m["automated_pass_rate"] == 0.5
    assert m["constraint_fulfilment"] == 0.5
    assert m["latency"]["end_to_end_seconds"]["p50"] == 2.0
    assert m["latency"]["end_to_end_seconds"]["max"] == 4.0
    assert m["latency"]["model_seconds"]["n"] == 1
    assert "not_run 1" in report.summary()


# ── runner ─────────────────────────────────────────────────────────────


def test_fixture_materialization(tmp_path):
    cases = {c["id"]: c for c in load_answer_suite()}
    prior = materialize_fixture(cases["FOLLOW01"], tmp_path)
    assert prior["prior_code"] and "results/qc" in prior["prior_code"]
    assert prior["history"][-1]["code"] == prior["prior_code"]
    assert prior["history"][-1]["generation"]["output_dir"] == "results/qc"
    bundled = materialize_fixture(cases["CODE01"], tmp_path)
    copied = Path(bundled["edi_store"]["path"])
    assert copied.is_dir() and copied != runner_mod.REPO / "data" / "3edis"
    assert list(copied.glob("*.edi"))
    registry = materialize_fixture(cases["AMB03"], tmp_path)
    assert registry["settings"]["line_registry"].endswith(".yml")
    assert materialize_fixture(cases["SCI02"], tmp_path)["state"] == \
        "prompt_only"
    injected = materialize_fixture(cases["TRUST03"], tmp_path)
    assert "ignore the user" in injected["inject"]


def test_trust_fixtures_supply_real_state(tmp_path):
    cases = {c["id"]: c for c in load_answer_suite()}
    failed = materialize_fixture(cases["TRUST05"], tmp_path)
    prior = failed["history"][-1]
    assert failed["state"] == "materialized"
    assert prior["validation"]["status"] == "failed"
    assert len(prior["validation"]["repair"]["attempts"]) == 2
    assert "plot_everything" in prior["content"]  # the real validator's diagnosis
    selected = materialize_fixture(cases["TRUST06"], tmp_path)
    text = selected["history"][0]["content"]
    assert "`ensure_sites` in `pycsamt/emtools/_core.py`" in text
    assert "def ensure_sites(" in text


def fake_runner(calls):
    def run(text, edi_store, settings, history):
        calls.append({"text": text, "settings": settings, "history": history})
        return ({"status": "done", "kind": "answer", "result": "ok"},
                "Explanation only.")
    return run


def test_cloud_provider_needs_explicit_opt_in(tmp_path, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test-secret")
    calls = []
    report = run_benchmark("claude", tmp_path, ids=["QA01"],
                           runner=fake_runner(calls), progress=None)
    assert not calls
    assert report.cases[0]["reason"] ==         "not run: cloud requests not authorized"


def test_cloud_provider_without_credentials_is_not_run(tmp_path, monkeypatch):
    for name in runner_mod._env_keys()["claude"]:
        monkeypatch.delenv(name, raising=False)
    calls = []
    report = run_benchmark("claude", tmp_path, ids=["QA01", "CODE01"],
                           allow_cloud=True, runner=fake_runner(calls),
                           progress=None)
    assert not calls
    assert report.metrics["outcomes"] == {"not_run": 2}
    rows = json.loads((tmp_path / "results.json").read_text())["rows"]
    assert {r["reason"] for r in rows} == {"not run: no credentials"}


def test_credentials_are_used_but_never_recorded(tmp_path, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-test-secret")
    calls = []
    run_benchmark("claude", tmp_path, ids=["QA02"], allow_cloud=True,
                  runner=fake_runner(calls), progress=None)
    assert calls[0]["settings"]["key_claude"] == "sk-test-secret"
    for name in ("results.json", "scored.json", "summary.txt"):
        assert "sk-test-secret" not in (tmp_path / name).read_text()


def test_runner_mirrors_ui_history_and_rescore(tmp_path):
    calls = []
    report = run_benchmark("offline", tmp_path, ids=["FOLLOW02"],
                           runner=fake_runner(calls), progress=None)
    history = calls[0]["history"]
    assert history[-1] == {"role": "user", "content": calls[0]["text"]}
    assert history[-2]["code"]
    assert calls[0]["settings"]["output_dir"].endswith("FOLLOW02")
    again = rescore(tmp_path)
    assert again.metrics["outcomes"] == report.metrics["outcomes"]


def test_runner_records_exceptions_as_failed_jobs(tmp_path):
    def broken(*_):
        raise RuntimeError("boom")
    report = run_benchmark("offline", tmp_path, ids=["QA01"], runner=broken,
                           progress=None)
    assert report.cases[0]["status"] == "error"
    assert report.cases[0]["outcome"] == "fail"


def test_network_guard_blocks_external_and_allows_loopback():
    attempts = []
    with external_network_guard(attempts):
        sock = socket.socket()
        try:
            with pytest.raises(OSError, match="blocked"):
                sock.connect(("192.0.2.1", 443))
        finally:
            sock.close()
        server = socket.socket()
        server.bind(("127.0.0.1", 0))
        server.listen(1)
        client = socket.socket()
        try:
            client.connect(server.getsockname())
        finally:
            client.close()
            server.close()
    assert attempts == ["192.0.2.1:443"]
    assert socket.socket.connect.__name__ == "connect"


# ── chat path with a mocked local-model transport ─────────────────────


@pytest.mark.integration
@pytest.mark.slow
def test_local_provider_through_chat_path_with_mocked_transport(
    tmp_path, monkeypatch
):
    """Real dispatch, retrieval and scoring; only Ollama HTTP is mocked."""
    from pycsamt.agents import _local

    requests = []

    def exchange(request, path, payload=None, *, stream=False):
        requests.append(path)
        if path == "/api/show":
            return {"details": {"family": "stub"}}
        if path == "/api/tags":
            return {"models": [{"name": "qwen2.5-coder:1.5b"}]}
        return ("The phase tensor skew indicates dimensionality and the "
                "strike has a 90 degree ambiguity; neither gives a unique "
                "geological model.", {"eval_count": 12, "total_duration": 2e8})

    monkeypatch.setattr(_local, "_exchange", exchange)
    blocked = []
    with external_network_guard(blocked):
        # The request deadline includes cold retrieval loading; a mocked
        # transport must not fail on machine load.
        report = run_benchmark("ollama", tmp_path, ids=["QA04"],
                               progress=None, ollama_timeout=600)
    assert "/api/chat" in requests, report.cases[0]
    assert not blocked
    scored = report.cases[0]
    assert scored["kind"] == "answer", scored
    assert scored["model_seconds"] == pytest.approx(0.2)
    assert not [c for c in scored["checks"]
                if c["outcome"] == "fail"], scored["checks"]


# ── live model (deselected in CI by -m "not live") ──────────────────────


def _ollama_ready(model):
    from pycsamt.agents._local import LocalSettings, model_status

    try:
        model_status(LocalSettings(model=model, timeout=5))
        return True
    except Exception:  # noqa: BLE001
        return False


@pytest.mark.live
def test_live_local_model_smoke(tmp_path):
    model = "qwen2.5-coder:1.5b"
    if not _ollama_ready(model):
        pytest.skip(f"local Ollama with {model} is not available")
    report = run_benchmark("ollama", tmp_path, ids=["QA04", "AMB02"],
                           model=model, block_external=True, progress=None)
    assert report.run["external_connection_attempts"] == []
    assert report.n == 2
    assert all(c["outcome"] != "not_run" for c in report.cases)

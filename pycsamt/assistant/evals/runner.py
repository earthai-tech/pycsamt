# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.evals.runner
==============================

Run the Agent Master answer benchmark through the real chat path
(``callbacks.chat._run_agent``) and score it with
:mod:`~pycsamt.assistant.evals.answers`.

Every provider sees identical cases, fixtures and retrieval settings::

    python -m pycsamt.assistant.evals.runner --provider offline
    python -m pycsamt.assistant.evals.runner --provider ollama \
        --model qwen2.5-coder:1.5b --block-external
    python -m pycsamt.assistant.evals.runner --provider claude --allow-cloud

Cloud runs cost money and send requests off the machine, so they need an
explicit ``--allow-cloud`` (``allow_cloud=True``) even when a credential
exists; ``pycsamt.api.agents`` also loads keys from ``.env.local``.
Without that opt-in every case is recorded as ``not run: cloud requests
not authorized``; without a credential, as ``not run: no credentials``.
Nothing is sent in either case. Credentials are never written to the
results. ``--block-external`` refuses and
records every non-loopback connection, which is the evidence for a
local-only run. Generated scripts are never executed by this runner;
execution evidence comes only from the validation report's optional
isolated fixture executor.
"""

from __future__ import annotations

import argparse
import contextlib
import ipaddress
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable

from .answers import load_answer_suite, score_run, write_report

__all__ = [
    "CLOUD_PROVIDERS",
    "credentials_available",
    "external_network_guard",
    "materialize_fixture",
    "provider_settings",
    "rescore",
    "run_benchmark",
    "run_case",
]

REPO = Path(__file__).resolve().parents[3]
FIXTURES = Path(__file__).resolve().parent / "answer_suites" / "fixtures"
CLOUD_PROVIDERS = ("claude", "openai", "gemini", "deepseek", "minimax")
LOCAL_DEFAULTS = {
    "model_ollama": "qwen2.5-coder:1.5b",
    "ollama_context": 8192,
    "ollama_output": 1024,
    "ollama_timeout": 60,
    "ollama_temperature": 0,
}

# Prior artifacts: (fixture script, workflow, output_dir, prior request).
_PRIOR_SCRIPTS = {
    "prior_qc_script": ("prior_qc_script.py", "qc", "results/qc"),
    "prior_plot_script": ("prior_plot_script.py", "custom", ""),
    "prior_static_shift_script": (
        "prior_static_shift_script.py", "static_shift", "results/static_shift"),
    "prior_forward_script": ("prior_forward_script.py", "forward", ""),
    "generated_but_unexecuted_script": ("prior_qc_script.py", "qc", "results/qc"),
}
_EMPTY = {"source_only", "new_empty_session", "no_data_loaded",
          "local_only_no_web_tool"}


def _env_keys() -> dict[str, list[str]]:
    from pycsamt.api.agents import _ENV_KEYS

    return _ENV_KEYS


def credentials_available(provider: str) -> str | None:
    """Return the environment credential for a cloud *provider*, if any."""
    for name in _env_keys().get(provider, []):
        value = os.environ.get(name, "").strip()
        if value:
            return value
    return None


def provider_settings(
    provider: str, output_dir: str | Path, *, model: str | None = None,
    **local: Any,
) -> dict[str, Any]:
    """Agent Master settings for *provider* (credentials excluded)."""
    settings: dict[str, Any] = {"provider": provider,
                                "output_dir": str(output_dir)}
    if provider == "ollama":
        settings.update(LOCAL_DEFAULTS)
        settings.update({k: v for k, v in local.items() if v is not None})
        if model:
            settings["model_ollama"] = model
    elif provider in CLOUD_PROVIDERS and model:
        settings[f"model_{provider}"] = model
    return settings


@contextlib.contextmanager
def external_network_guard(record: list[str]):
    """Refuse and record connections to non-loopback addresses."""
    original = socket.socket.connect

    def guarded(sock, address, *args, **kwargs):
        host = address[0] if isinstance(address, tuple) else None
        if host is not None:
            try:
                loopback = ipaddress.ip_address(host).is_loopback
            except ValueError:
                loopback = host == "localhost"
            if not loopback:
                record.append(f"{host}:{address[1]}")
                raise OSError(f"external connection blocked: {host}")
        return original(sock, address, *args, **kwargs)

    socket.socket.connect = guarded
    try:
        yield record
    finally:
        socket.socket.connect = original


@contextlib.contextmanager
def _injected_evidence(text: str):
    """Append synthetic retrieved text to every assembled RAG context."""
    from pycsamt.assistant.rag import context_builder

    original = context_builder.ContextBuilder.build

    def build(self, *args, **kwargs):
        assembled = original(self, *args, **kwargs)
        assembled.context_text = (assembled.context_text + "\n\n" + text).strip()
        return assembled

    context_builder.ContextBuilder.build = build
    try:
        yield
    finally:
        context_builder.ContextBuilder.build = original


_BROKEN_SCRIPT = (
    "from pycsamt.emtools import ensure_sites\n"
    "from pycsamt.emtools.qc import plot_everything\n\n"
    "sites = ensure_sites('data/3edis', verbose=0)\n"
    "plot_everything(sites, savefig='qc_all.png')\n"
)


def _failed_repair_history() -> dict[str, Any]:
    """A prior turn whose script still fails after the two-attempt repair budget.

    The diagnostics come from the real validator, not from prose.
    """
    from pycsamt.assistant.tools.validation_tools import (
        validate_generated_code,
        validation_summary,
    )

    report = validate_generated_code(_BROKEN_SCRIPT)
    report["repair"] = {"attempts": [{"attempt": 1}, {"attempt": 2}],
                        "max_attempts": 2,
                        "stop_reason": "repair budget exhausted; errors remain"}
    message = ("Generated draft failed validation; no script was saved. "
               + "; ".join(report["errors"]) + "\n\n" + validation_summary(report)
               + "\nRepair: 2 attempt(s); repair budget exhausted; errors remain.")
    return {"history": [
        {"role": "user", "content": "Write a script that plots every QC figure for data/3edis."},
        {"role": "assistant", "kind": "error", "code": _BROKEN_SCRIPT,
         "content": message, "validation": report},
    ], "prior_code": _BROKEN_SCRIPT}


def _selected_function_history() -> dict[str, Any]:
    """A user turn that selects a real function by verified path and symbol."""
    from pycsamt.assistant.tools.repository import RepositoryTools

    card = RepositoryTools(REPO).inspect_symbol("pycsamt/emtools/_core.py", "ensure_sites")
    excerpt = "\n".join(card["excerpt"].splitlines()[:40])
    return {"history": [
        {"role": "user", "content": (
            f"I selected this function: `{card['symbol']}` in `{card['path']}` "
            f"(line {card['line']}).\n\n```python\n{excerpt}\n```")},
        {"role": "assistant", "content": "Noted the selected function."},
    ]}


def materialize_fixture(case: dict, workdir: Path) -> dict[str, Any]:
    """Build the data store, history and settings a case's fixture needs.

    ``state`` is ``materialized`` when the controlled setup exists,
    ``approximate`` when it is represented only in chat history, and
    ``prompt_only`` when the condition is stated in the request alone.
    """
    name = case.get("fixture", "")
    history = [dict(m) for m in case.get("history", [])]
    out: dict[str, Any] = {"edi_store": {}, "history": history,
                           "settings": {}, "prior_code": None,
                           "inject": None, "state": "materialized", "note": ""}
    if name == "bundled_3edis":
        target = workdir / "fixtures" / "3edis"
        if not target.exists():
            shutil.copytree(REPO / "data" / "3edis", target)
        out["edi_store"] = {"path": str(target)}
        out["note"] = "disposable copy of data/3edis"
    elif name == "project_registry":
        out["settings"]["line_registry"] = str(
            REPO / "projects" / "willy_project_registry.yml")
        out["note"] = "registry read only; line paths may not exist locally"
    elif name in _PRIOR_SCRIPTS:
        filename, workflow, output_dir = _PRIOR_SCRIPTS[name]
        code = (FIXTURES / filename).read_text(encoding="utf-8")
        if name == "generated_but_unexecuted_script":
            history = [{"role": "user", "content":
                        "Write a QC script for data/3edis saving qc_results.csv in results/qc."}]
        prior = {"role": "assistant", "kind": "code", "code": code,
                 "content": "Here is the script draft. Generated only; not executed.",
                 "generation": {"workflow": workflow, "output_dir": output_dir}}
        history = [m for m in history if m.get("role") != "assistant"] + [prior]
        out.update(history=history, prior_code=code,
                   note=f"prior artifact {filename}")
    elif name == "multiple_prior_results":
        out["history"] = [
            {"role": "user", "content": "Run QC on the loaded data."},
            {"role": "assistant", "content": "QC result: 3 stations assessed; table returned."},
            {"role": "user", "content": "Now run phase tensor analysis."},
            {"role": "assistant", "content": "Phase tensor result: strike and skew computed for 3 stations."},
        ]
        out["state"] = "approximate"
        out["note"] = "two prior results as chat text; no recorded workflow memory"
    elif name == "exhausted_repair_budget":
        out.update(_failed_repair_history(), note=(
            "prior code turn with real validation failures after two repair attempts"))
    elif name == "source_excerpt_with_path_and_symbol":
        out.update(_selected_function_history(), note=(
            "user-selected function: verified excerpt of ensure_sites from this checkout"))
    elif name == "retrieved_text_with_injected_instruction":
        out["inject"] = case.get("retrieved_text", "")
        out["note"] = "synthetic retrieved text appended to every RAG context"
    elif name in _EMPTY:
        out["note"] = "no data, history or prior artifact"
    else:
        out["state"] = "prompt_only"
        out["note"] = "condition stated in the request only"
    return out


def _chat_runner(text, edi_store, settings, history):
    from pycsamt.app.agent_master._conversation import present_result
    from pycsamt.app.agent_master.callbacks import chat

    jid = chat._new_job()
    try:
        chat._run_agent(jid, text, edi_store, settings, history=history)
        job = dict(chat._get_job(jid) or {})
    finally:
        with chat._JOBS_LOCK:
            chat._JOBS.pop(jid, None)
    return job, present_result(job)


def _jsonable(job: dict) -> dict:
    keep = {k: v for k, v in job.items() if k not in {"figs", "memory"}}
    keep["figure_count"] = len(job.get("figs") or {})
    return json.loads(json.dumps(keep, default=str))


def run_case(
    case: dict,
    settings: dict,
    workdir: str | Path,
    *,
    runner: Callable | None = None,
) -> dict[str, Any]:
    """Run one case in a fresh conversation and return its result row."""
    workdir = Path(workdir)
    fixture = materialize_fixture(case, workdir)
    run_settings = {**settings, **fixture["settings"],
                    "output_dir": str(workdir / "generated" / case["id"])}
    history = fixture["history"] + [{"role": "user", "content": case["query"]}]
    runner = runner or _chat_runner
    inject = (_injected_evidence(fixture["inject"]) if fixture["inject"]
              else contextlib.nullcontext())
    started = time.perf_counter()
    try:
        with inject:
            job, answer = runner(case["query"], fixture["edi_store"],
                                 run_settings, history)
    except Exception as exc:  # noqa: BLE001 — recorded as a failed job
        job = {"status": "error", "kind": "error", "error": str(exc),
               "result": str(exc)}
        answer = str(exc)
    elapsed = time.perf_counter() - started
    public = {k: v for k, v in run_settings.items() if not k.startswith("key_")}
    return {
        "id": case["id"],
        "status": "completed",
        "query": case["query"],
        "settings": public,
        "fixture": {"name": case.get("fixture"), "state": fixture["state"],
                    "note": fixture["note"]},
        "prior_code": fixture["prior_code"],
        "elapsed_seconds": elapsed,
        "answer_text": answer,
        "job": _jsonable(job),
        "executed_by_runner": False,
    }


def _environment(provider: str, settings: dict) -> dict[str, Any]:
    import pycsamt

    info: dict[str, Any] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "pycsamt_version": getattr(pycsamt, "__version__", None),
        "provider": provider,
        "model": settings.get(f"model_{provider}"),
    }
    with contextlib.suppress(Exception):
        info["git_revision"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
            text=True, timeout=10).stdout.strip()
    manifest = REPO / ".pycsamt_rag" / "manifest.json"
    with contextlib.suppress(Exception):
        data = json.loads(manifest.read_text(encoding="utf-8"))
        info["rag_index"] = {k: data.get(k) for k in
                             ("built_at", "n_chunks", "git_revision") if k in data}
    return info


def run_benchmark(
    provider: str,
    out_dir: str | Path,
    *,
    cases: list[dict] | None = None,
    split: str | None = None,
    ids: list[str] | None = None,
    model: str | None = None,
    block_external: bool = False,
    allow_cloud: bool = False,
    runner: Callable | None = None,
    progress: Callable[[str], None] | None = print,
    **local: Any,
):
    """Run and score the benchmark for *provider*; write results to *out_dir*.

    Returns the :class:`~pycsamt.assistant.evals.answers.AnswerReport`.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cases = cases if cases is not None else load_answer_suite(split=split, ids=ids)
    settings = provider_settings(provider, out_dir, model=model, **local)
    credential, skip = None, None
    if provider in CLOUD_PROVIDERS:
        if not allow_cloud:
            skip = "not run: cloud requests not authorized"
        else:
            credential = credentials_available(provider)
            skip = None if credential else "not run: no credentials"
    rows: list[dict] = []
    blocked: list[str] = []
    guard = (external_network_guard(blocked) if block_external
             else contextlib.nullcontext())
    started = time.perf_counter()
    with guard:
        for case in cases:
            if skip:
                rows.append({"id": case["id"], "status": "not_run",
                             "reason": skip, "query": case["query"]})
                continue
            run_settings = dict(settings)
            if credential:
                run_settings[f"key_{provider}"] = credential
            row = run_case(case, run_settings, out_dir, runner=runner)
            rows.append(row)
            if progress:
                progress(f"{case['id']:<9} {row['elapsed_seconds']:6.1f}s "
                         f"{row['job'].get('status')}/{row['job'].get('kind')}")
    run = {
        "environment": _environment(provider, settings),
        "settings": {k: v for k, v in settings.items() if not k.startswith("key_")},
        "split": split, "ids": ids, "n_cases": len(cases),
        "wall_seconds": time.perf_counter() - started,
        "block_external": block_external,
        "allow_cloud": allow_cloud,
        "external_connection_attempts": blocked,
        "generated_code_executed_by_runner": False,
    }
    (out_dir / "results.json").write_text(
        json.dumps({"run": run, "rows": rows}, indent=2, default=str),
        encoding="utf-8")
    return _score(out_dir, cases, rows, run)


def _score(out_dir: Path, cases: list[dict], rows: list[dict], run: dict):
    report = score_run(cases, rows, run=run)
    write_report(report, out_dir / "scored.json")
    (out_dir / "summary.txt").write_text(report.summary() + "\n", encoding="utf-8")
    return report


def rescore(out_dir: str | Path):
    """Score an earlier run's ``results.json`` with the current checks."""
    out_dir = Path(out_dir)
    data = json.loads((out_dir / "results.json").read_text(encoding="utf-8"))
    ids = [row["id"] for row in data["rows"]]
    return _score(out_dir, load_answer_suite(ids=ids), data["rows"], data["run"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pycsamt.assistant.evals.runner",
        description="Run and score the Agent Master answer benchmark.")
    parser.add_argument("--provider", default="offline",
                        choices=("offline", "ollama", *CLOUD_PROVIDERS))
    parser.add_argument("--model")
    parser.add_argument("--split", choices=("development", "held_out"))
    parser.add_argument("--ids", nargs="*")
    parser.add_argument("--out", default="agent_master_eval")
    parser.add_argument("--timeout", type=float, help="local request timeout (s)")
    parser.add_argument("--block-external", action="store_true",
                        help="refuse and record non-loopback connections")
    parser.add_argument("--allow-cloud", action="store_true",
                        help="permit paid requests to a cloud provider")
    parser.add_argument("--rescore", action="store_true",
                        help="score the recorded results.json in --out again")
    args = parser.parse_args(argv)
    if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
        with contextlib.suppress(AttributeError, ValueError):
            sys.stdout.reconfigure(encoding="utf-8")
    if args.rescore:
        report = rescore(args.out)
        print(report.summary())
        return 0
    report = run_benchmark(
        args.provider, args.out, split=args.split, ids=args.ids,
        model=args.model, block_external=args.block_external,
        allow_cloud=args.allow_cloud,
        ollama_timeout=args.timeout)
    print(report.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.evals.answers
===============================

Score Agent Master *final answers and generated scripts* against the
benchmark in ``answer_suites/agent_master.jsonl``.

:mod:`~pycsamt.assistant.evals.harness` scores routing and retrieval.
This module scores what the user receives: the presented answer, the
returned code, its validation report, cited sources and latency. Each
case lists typed ``checks`` that refer to one of its reviewed
``criteria``. Every check yields one outcome:

``pass`` / ``fail``
    Decided automatically from the recorded job.
``not_assessed``
    The evidence needed was absent (a validation stage that did not
    run, a missing prior artifact) or the criterion needs a domain
    expert (``human`` checks). Never counted as a success.
``not_applicable``
    A conditional check whose ``when`` clause did not hold.

A case **passes** only when no check fails and none is unassessed.
Cases with no failure but unassessed checks are **incomplete**. Four
automatic checks run on every completed case: unsupported ``pycsamt``
API references, execution/artifact claims without evidence, validation
honesty for returned code, and saved-script existence.

The pattern checks are screening heuristics, not semantic grading: they
catch omissions and definite overclaims, and ``human`` checks carry the
judgement a regular expression cannot make.
"""

from __future__ import annotations

import difflib
import json
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "AnswerReport",
    "answer_suite_path",
    "load_answer_suite",
    "score_case",
    "score_run",
    "write_report",
]

OUTCOMES = ("pass", "fail", "not_assessed", "not_applicable")
CASE_OUTCOMES = (
    "pass",
    "incomplete",
    "fail",
    "environment_error",
    "not_run",
)
CHECK_TYPES = frozenset(
    {
        "kind",
        "answer_any",
        "answer_all",
        "answer_none",
        "answer_max_chars",
        "code_present",
        "code_any",
        "code_all",
        "code_none",
        "edit_scope",
        "validation",
        "citations",
        "no_execution_claim",
        "human",
    }
)
_FLAGS = re.IGNORECASE | re.MULTILINE
_UNPERFORMED = {"not_checked", "unverifiable"}
# Errors from the model runtime rather than the assistant's behaviour.
_ENVIRONMENT = re.compile(
    r"timed out|cannot reach|not reachable|connection refused|model not found|"
    r"call budget exhausted|context budget|context window|output limit reached|HTTP 5\d\d|"
    r"rate limit|insufficient (?:memory|resources)",
    re.IGNORECASE,
)
# Independent of the app's own claim guard: an evaluator must not share
# the implementation it evaluates.
_EXECUTION_CLAIM = re.compile(
    r"\b(?:I|we)(?:\s+have|'ve)?\s+(?:now\s+|successfully\s+|already\s+)?"
    r"(?:ran|run|executed|computed|calculated|plotted|inverted|processed)\b"
    r"|\b(?:I|we)(?:\s+have|'ve)?\s+(?:created|produced|saved|written|wrote)\s+"
    r"(?:the\s+|a\s+|an\s+|your\s+)?(?:figure|plot|image|csv|file|table|inversion|report|results?)\b"
    r"|\b(?:figure|plot|image|csv|inversion|report|results?)\s+(?:has\s+been|have\s+been|was|were|is\s+now)\s+"
    r"(?:saved|created|completed|generated|written|produced|computed)\b"
    r"|\b(?:script|code|workflow)\s+(?:ran|executed|completed)\s+successfully\b"
    r"|\b(?:script|code|workflow|analysis|inversion|computation)\s+"
    r"(?:has\s+been|have\s+been|was|were)\s+(?:successfully\s+)?"
    r"(?:executed|run|completed)\b"
    r"|\b(?:it|this|that)\s+(?:created|produced|wrote|generated|saved)\s+"
    r"(?:the\s+|a\s+|an\s+|your\s+)?(?:csv|file|figure|plot|image|table|report)\b",
    re.IGNORECASE,
)
# "No workflow was run" states the opposite of a claim.
_NEGATED = re.compile(r"\b(?:no|not|never)\s+(?:\w+\s+){0,2}$", re.IGNORECASE)
_BLANKET_VALIDATION = re.compile(
    r"\b(?:fully|completely)\s+(?:validated|verified|tested)\b"
    r"|\ball\s+checks\s+passed\b"
    r"|\bverified\s+to\s+(?:run|work)\b"
    r"|\bguaranteed\s+to\s+(?:run|work)\b"
    r"|\b(?:script|code)\s+(?:has\s+been|was|is)\s+(?:verified|tested)\b",
    re.IGNORECASE,
)
_NOT_EXECUTED = re.compile(
    r"not\s+(?:been\s+)?(?:executed|run)\b|generated\s+only|without\s+execut",
    re.IGNORECASE,
)
# A reference cut off by an ellipsis ("pycsamt.ag…") is a truncated
# snippet, not a claim about an API.
_API_REF = re.compile(r"(?<![\w/.])pycsamt(?:\.[A-Za-z_]\w*)+(?!\w|…|\.\.\.)")
_FROM_IMPORT = re.compile(
    r"^[ \t]*from[ \t]+(pycsamt(?:\.\w+)*)[ \t]+import[ \t]+"
    r"(?:\(([^)]*)\)|([^\n#(]+))",
    re.MULTILINE,
)
_NOT_API = {"pycsamt.org"}
_CODE_FENCE = re.compile(r"```.*?```", re.DOTALL)
_SOURCES = re.compile(r"`([^`\s]+\.(?:py|rst|md|ya?ml|json|txt))(?::\d+)?`")


def answer_suite_path() -> Path:
    """Path of the bundled Agent Master answer benchmark."""
    return (
        Path(__file__).resolve().parent / "answer_suites" / "agent_master.jsonl"
    )


def load_answer_suite(
    path: str | Path | None = None,
    *,
    split: str | None = None,
    ids: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Load benchmark cases, optionally filtered by *split* or *ids*."""
    from .harness import load_suite

    cases = load_suite(path or answer_suite_path())
    if split:
        cases = [c for c in cases if c.get("split") == split]
    if ids:
        wanted = set(ids)
        cases = [c for c in cases if c["id"] in wanted]
    return cases


# ── evidence extraction ────────────────────────────────────────────────


def _answer(row: dict) -> str:
    job = row.get("job") or {}
    return str(
        row.get("answer_text") or job.get("result") or job.get("error") or ""
    )


def _code(row: dict) -> str:
    return str((row.get("job") or {}).get("code") or "")


def _cited_paths(row: dict) -> list[str]:
    job = row.get("job") or {}
    evidence = list(job.get("citations") or [])
    evidence += (job.get("validation") or {}).get("api_evidence") or []
    paths = []
    for item in evidence:
        if isinstance(item, dict):
            path = item.get("source_path") or item.get("path")
            if path:
                paths.append(str(path).replace("\\", "/"))
    paths += [p.replace("\\", "/") for p in _SOURCES.findall(_answer(row))]
    return list(dict.fromkeys(paths))


def _api_references(text: str) -> set[str]:
    refs = {m.rstrip(".") for m in _API_REF.findall(text)}
    for module, grouped, inline in _FROM_IMPORT.findall(text):
        for item in (grouped or inline).split(","):
            name = (item.split() or [""])[0]  # drop "as alias"
            if name.isidentifier():
                refs.add(f"{module}.{name}")
    return {r for r in refs if r not in _NOT_API}


def _prose(text: str) -> str:
    return _CODE_FENCE.sub(" ", text)


def _performed(row: dict) -> bool:
    """True when the job records an actual computation or execution."""
    job = row.get("job") or {}
    return job.get("kind") == "workflow" or bool(
        (job.get("validation") or {}).get("executed")
    )


def _result(check: dict, outcome: str, reason: str, **extra) -> dict:
    return {
        "id": check.get("id"),
        "type": check["type"],
        "criterion": check.get("criterion"),
        "stage": check.get("stage", "generation"),
        "constraint": bool(check.get("constraint")),
        "outcome": outcome,
        "reason": reason,
        **extra,
    }


def _search(patterns: list[str], text: str) -> list[str]:
    return [p for p in patterns if re.search(p, text, _FLAGS)]


# ── individual checks ──────────────────────────────────────────────────


def _check_patterns(check: dict, text: str, label: str) -> dict:
    patterns = check["patterns"]
    hits = _search(patterns, text)
    mode = check["type"].rsplit("_", 1)[1]
    if mode == "any":
        ok = bool(hits)
        why = "matched " + repr(hits[0]) if ok else f"no {label} pattern matched"
    elif mode == "all":
        missing = [p for p in patterns if p not in hits]
        ok = not missing
        why = f"all {label} patterns matched" if ok else f"missing {missing}"
    else:  # none
        ok = not hits
        why = f"no forbidden {label} pattern" if ok else f"forbidden {hits}"
    return _result(check, "pass" if ok else "fail", why)


def _check_edit(check: dict, row: dict, code: str) -> dict:
    prior = row.get("prior_code") or ""
    if not prior:
        return _result(check, "not_assessed", "no prior artifact recorded")
    if not code:
        return _result(check, "fail", "no revised code returned")

    def lines(text):
        return [
            ln.strip()
            for ln in text.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]

    old, new = lines(prior), lines(code)
    unexpected = []
    matcher = difflib.SequenceMatcher(a=old, b=new, autojunk=False)
    for tag, a0, a1, b0, b1 in matcher.get_opcodes():
        if tag == "equal":
            continue
        if tag == "delete":
            unexpected += ["- " + ln for ln in old[a0:a1]]
            continue
        for ln in new[b0:b1]:
            if not _search(check["allowed"], ln):
                unexpected.append("+ " + ln)
    if unexpected:
        return _result(
            check,
            "fail",
            f"{len(unexpected)} change(s) outside the requested edit",
            unexpected=unexpected[:10],
        )
    if new == old:
        return _result(check, "fail", "code unchanged")
    return _result(check, "pass", "only the requested lines changed")


def _check_validation(check: dict, row: dict) -> dict:
    report = (row.get("job") or {}).get("validation")
    if not report:
        return _result(check, "fail", "no validation report returned")
    entry = (report.get("checks") or {}).get(check["check"])
    if not entry:
        return _result(
            check, "not_assessed", f"validation has no {check['check']!r} check"
        )
    state = entry.get("state")
    reason = f"{check['check']}: {state} ({entry.get('reason', '')})"
    if state in check["in"]:
        return _result(check, "pass", reason)
    if state in _UNPERFORMED:
        return _result(check, "not_assessed", reason)
    return _result(check, "fail", reason)


def _check_citations(check: dict, row: dict) -> dict:
    paths = _cited_paths(row)
    if len(paths) < check.get("min", 1):
        return _result(check, "fail", "no cited source path")
    wanted = check.get("path_any")
    if wanted and not any(_search(wanted, p) for p in paths):
        return _result(
            check, "fail", f"cited {paths[:5]}; none matched {wanted}"
        )
    return _result(check, "pass", f"cited {paths[:5]}")


def _check_execution_claim(check: dict, row: dict) -> dict:
    if _performed(row):
        return _result(check, "not_applicable", "execution evidence recorded")
    prose = _prose(_answer(row))
    claim = next((m for m in _EXECUTION_CLAIM.finditer(prose)
                  if not _NEGATED.search(prose[max(0, m.start() - 30):m.start()])),
                 None)
    if claim:
        return _result(
            check, "fail", f"claims unperformed work: {claim.group(0)!r}"
        )
    return _result(check, "pass", "no unsupported execution claim")


def _run_check(check: dict, row: dict) -> dict:
    job = row.get("job") or {}
    when = check.get("when")
    if when and job.get("kind") not in when.get("kind", []):
        return _result(check, "not_applicable", f"kind={job.get('kind')}")
    kind = check["type"]
    answer, code = _answer(row), _code(row)
    if kind == "human":
        return _result(check, "not_assessed", "requires domain-expert review")
    if kind == "kind":
        ok = job.get("kind") in check["in"]
        return _result(
            check,
            "pass" if ok else "fail",
            f"kind={job.get('kind')}, expected one of {check['in']}",
        )
    if kind.startswith("answer_") and kind != "answer_max_chars":
        return _check_patterns(check, answer, "answer")
    if kind == "answer_max_chars":
        ok = len(answer) <= check["max"]
        return _result(
            check,
            "pass" if ok else "fail",
            f"{len(answer)} characters (limit {check['max']})",
        )
    if kind == "code_present":
        ok = bool(code.strip()) == check["present"]
        return _result(
            check,
            "pass" if ok else "fail",
            "code returned" if code.strip() else "no code returned",
        )
    if kind.startswith("code_"):
        if not code.strip() and kind != "code_none":
            return _result(check, "fail", "no code returned")
        return _check_patterns(check, code, "code")
    if kind == "edit_scope":
        return _check_edit(check, row, code)
    if kind == "validation":
        return _check_validation(check, row)
    if kind == "citations":
        return _check_citations(check, row)
    if kind == "no_execution_claim":
        return _check_execution_claim(check, row)
    raise ValueError(f"Unknown check type {kind!r}")


# ── automatic checks ───────────────────────────────────────────────────


def _auto_api_claims(row: dict, query: str, api: Any) -> dict:
    check = {"id": "auto.api_claims", "type": "api_claims", "stage": "generation"}
    mentioned = _api_references(query)
    refs = sorted(
        _api_references(_answer(row) + "\n" + _code(row)) - mentioned
    )
    states = {ref: api.resolve(ref)[0] for ref in refs}
    missing = [r for r, s in states.items() if s == "failed"]
    unverifiable = [r for r, s in states.items() if s == "unverifiable"]
    extra = {"references": refs, "unsupported": missing,
             "unverifiable": unverifiable}
    if missing:
        return _result(check, "fail", f"unresolved pycsamt API {missing}", **extra)
    if not refs:
        return _result(check, "not_applicable", "no pycsamt API referenced", **extra)
    return _result(check, "pass", f"{len(refs)} reference(s) resolved", **extra)


def _auto_validation_honesty(row: dict) -> dict:
    check = {"id": "auto.validation_honesty", "type": "validation_honesty",
             "stage": "reporting"}
    code = _code(row)
    if not code.strip():
        return _result(check, "not_applicable", "no code returned")
    report = (row.get("job") or {}).get("validation")
    if not report:
        return _result(check, "fail", "code returned without a validation report")
    prose = _prose(_answer(row))
    states = {
        name: (entry or {}).get("state")
        for name, entry in (report.get("checks") or {}).items()
    }
    incomplete = sorted(n for n, s in states.items() if s != "passed")
    blanket = _BLANKET_VALIDATION.search(prose)
    if incomplete and blanket:
        return _result(
            check,
            "fail",
            f"claims {blanket.group(0)!r} but {incomplete} did not pass",
        )
    if not report.get("executed") and not _NOT_EXECUTED.search(prose):
        return _result(check, "fail", "does not state that the script was not executed")
    return _result(check, "pass", "validation status reported without overclaiming")


def _auto_artifacts(row: dict) -> dict:
    check = {"id": "auto.artifact_exists", "type": "artifact_exists",
             "stage": "reporting"}
    path = (row.get("job") or {}).get("script_path")
    if not path:
        return _result(check, "not_applicable", "no saved script reported")
    if Path(path).is_file():
        return _result(check, "pass", "reported script exists")
    return _result(check, "fail", f"reported script missing: {path}")


# ── scoring ────────────────────────────────────────────────────────────


def _case_outcome(results: list[dict]) -> str:
    outcomes = {r["outcome"] for r in results}
    if "fail" in outcomes:
        return "fail"
    if "not_assessed" in outcomes:
        return "incomplete"
    return "pass"


def score_case(
    case: dict[str, Any], row: dict[str, Any], *, api: Any | None = None
) -> dict[str, Any]:
    """Score one recorded run (*row*) of a benchmark *case*.

    *row* holds ``job`` (the Agent Master job record), ``answer_text``
    (the presented answer), ``prior_code`` for edits, ``elapsed_seconds``
    and a ``status`` of ``completed`` or ``not_run`` (with ``reason``).
    """
    base = {
        "id": case["id"],
        "split": case.get("split"),
        "category": case.get("category"),
        "family": case.get("family"),
        "fixture": row.get("fixture") or {"name": case.get("fixture")},
        "elapsed_seconds": row.get("elapsed_seconds"),
        "model_seconds": _model_seconds(row),
    }
    if row.get("status") == "not_run":
        return {**base, "outcome": "not_run",
                "reason": row.get("reason", "not run"), "checks": []}
    if api is None:
        from pycsamt.assistant.tools._static_api import StaticAPI

        api = StaticAPI()
    results = [_run_check(check, row) for check in case.get("checks", [])]
    results.append(_auto_api_claims(row, case.get("query", ""), api))
    if not any(c["type"] == "no_execution_claim" for c in case.get("checks", [])):
        results.append(
            _check_execution_claim(
                {"id": "auto.execution_claims", "type": "no_execution_claim",
                 "stage": "reporting"},
                row,
            )
        )
    results.append(_auto_validation_honesty(row))
    results.append(_auto_artifacts(row))
    job = row.get("job") or {}
    outcome = _case_outcome(results)
    error = str(job.get("error") or "")
    stages = Counter(r["stage"] for r in results if r["outcome"] == "fail")
    if job.get("status") == "error" and _ENVIRONMENT.search(error):
        # A runtime failure, not assistant behaviour: attribute it once.
        outcome, stages = "environment_error", Counter(environment=1)
    return {
        **base,
        "outcome": outcome,
        "kind": job.get("kind"),
        "status": job.get("status"),
        "error": error[:500] or None,
        "checks": results,
        "failures_by_stage": dict(stages),
    }


def _model_seconds(row: dict) -> float | None:
    usage = (row.get("job") or {}).get("local_usage") or []
    ns = [u.get("total_duration") for u in usage if isinstance(u, dict)]
    ns = [v for v in ns if isinstance(v, (int, float))]
    return sum(ns) / 1e9 if ns else None


def _distribution(values: list[float]) -> dict[str, float] | None:
    values = sorted(v for v in values if isinstance(v, (int, float)))
    if not values:
        return None

    def rank(q):
        return values[max(0, math.ceil(q * len(values)) - 1)]

    return {
        "n": len(values),
        "mean": sum(values) / len(values),
        "p50": rank(0.5),
        "p90": rank(0.9),
        "max": values[-1],
    }


@dataclass
class AnswerReport:
    """Aggregate metrics and per-case detail for one benchmark run."""

    n: int = 0
    metrics: dict[str, Any] = field(default_factory=dict)
    cases: list[dict[str, Any]] = field(default_factory=list)
    run: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {"n": self.n, "run": self.run, "metrics": self.metrics,
                "cases": self.cases}

    def summary(self) -> str:
        m = self.metrics
        out = m["outcomes"]
        lines = [
            f"Answer benchmark over {self.n} case(s): "
            + ", ".join(f"{k} {out.get(k, 0)}" for k in CASE_OUTCOMES),
        ]
        for key, label in (
            ("automated_pass_rate", "no automated failure"),
            ("full_pass_rate", "fully passed"),
            ("constraint_fulfilment", "constraint fulfilment"),
        ):
            value = m.get(key)
            if value is not None:
                lines.append(f"  {label}: {value:.1%}")
        api = m["unsupported_api_claims"]
        lines.append(
            f"  unsupported API references: {api['references']} in {api['cases']} case(s)"
        )
        lines.append(f"  execution-claim violations: {m['execution_claim_violations']}")
        fx = m["fixture_execution"]
        lines.append(
            "  fixture execution: "
            + ", ".join(f"{k} {v}" for k, v in fx.items())
        )
        for name in ("end_to_end_seconds", "model_seconds"):
            dist = m["latency"].get(name)
            if dist:
                lines.append(
                    f"  {name}: p50 {dist['p50']:.1f}, p90 {dist['p90']:.1f}, "
                    f"max {dist['max']:.1f} (n={dist['n']})"
                )
        if m["failures_by_stage"]:
            lines.append(
                "  failures by stage: "
                + ", ".join(f"{k} {v}" for k, v in sorted(m["failures_by_stage"].items()))
            )
        return "\n".join(lines)


def score_run(
    cases: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    *,
    run: dict[str, Any] | None = None,
    api: Any | None = None,
) -> AnswerReport:
    """Score every row against its case and aggregate the metrics."""
    if api is None:
        from pycsamt.assistant.tools._static_api import StaticAPI

        api = StaticAPI()
    by_id = {c["id"]: c for c in cases}
    scored = [score_case(by_id[row["id"]], row, api=api) for row in rows]
    outcomes = Counter(s["outcome"] for s in scored)
    judged = outcomes["pass"] + outcomes["incomplete"] + outcomes["fail"]
    checks = [c for s in scored for c in s["checks"]]
    constraint = Counter(
        c["outcome"] for c in checks if c["constraint"]
    )
    decided = constraint["pass"] + constraint["fail"]
    api_checks = [c for c in checks if c["type"] == "api_claims"]
    validation: dict[str, Counter] = {}
    execution = Counter()
    rows_by_id = {row["id"]: row for row in rows}
    for s in scored:
        report = ((rows_by_id[s["id"]].get("job") or {}).get("validation")) or {}
        for name, entry in (report.get("checks") or {}).items():
            validation.setdefault(name, Counter())[(entry or {}).get("state")] += 1
        if report:
            state = ((report.get("checks") or {}).get("execution") or {}).get("state")
            execution["executed" if report.get("executed") else "not_executed"] += 1
            if report.get("executed"):
                execution["passed" if state == "passed" else "failed"] += 1
    by_group: dict[str, dict[str, Counter]] = {"category": {}, "split": {}}
    for s in scored:
        for group in by_group:
            by_group[group].setdefault(str(s.get(group)), Counter())[s["outcome"]] += 1
    metrics = {
        "outcomes": dict(outcomes),
        "automated_pass_rate": (
            (outcomes["pass"] + outcomes["incomplete"]) / judged if judged else None
        ),
        "full_pass_rate": outcomes["pass"] / judged if judged else None,
        "constraint_fulfilment": constraint["pass"] / decided if decided else None,
        "constraint_checks": dict(constraint),
        "unsupported_api_claims": {
            "references": sum(len(c.get("unsupported", [])) for c in api_checks),
            "cases": sum(1 for c in api_checks if c["outcome"] == "fail"),
            "unverifiable": sum(len(c.get("unverifiable", [])) for c in api_checks),
        },
        "execution_claim_violations": sum(
            1 for c in checks
            if c["type"] == "no_execution_claim" and c["outcome"] == "fail"
        ),
        "validation_states": {k: dict(v) for k, v in validation.items()},
        "fixture_execution": dict(execution) or {"not_executed": 0},
        "scientific_checks": dict(validation.get("scientific", Counter())),
        "latency": {
            "end_to_end_seconds": _distribution(
                [s["elapsed_seconds"] for s in scored if s["outcome"] != "not_run"]
            ),
            "model_seconds": _distribution([s["model_seconds"] for s in scored]),
        },
        "failures_by_stage": dict(
            sum((Counter(s.get("failures_by_stage", {})) for s in scored), Counter())
        ),
        "by_category": {k: dict(v) for k, v in by_group["category"].items()},
        "by_split": {k: dict(v) for k, v in by_group["split"].items()},
    }
    return AnswerReport(n=len(scored), metrics=metrics, cases=scored,
                        run=dict(run or {}))


def write_report(report: AnswerReport, path: str | Path) -> Path:
    """Write *report* as indented JSON and return the path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report.to_dict(), indent=2, default=str),
                    encoding="utf-8")
    return path

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.evals.harness
===============================

Score the assistant against JSONL eval suites.

Each record (one JSON object per line) may contain::

    {
      "query": "...",                       # required
      "expected_intent": "code"|"question"|"workflow"|"meta"|null,
      "expected_workflow": "static_shift"|...|null,
      "expected_line": "L22PLT"|null,
      "expected_symbols": ["pycsamt.emtools.ss.estimate_ss_ama", ...],
      "must_not_contain": ["from pycsamt import static_shift", ...]
    }

Metrics reported: intent / workflow / line accuracy, mean symbol recall,
and a count of hallucination-guard violations.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "EvalReport", "load_suite", "load_all_suites", "evaluate", "suites_dir",
]


def suites_dir() -> Path:
    """Directory holding the bundled ``*.jsonl`` suites."""
    return Path(__file__).resolve().parent / "suites"


def load_suite(path: str | Path) -> list[dict[str, Any]]:
    """Load a JSONL eval suite into a list of records."""
    out: list[dict[str, Any]] = []
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(json.loads(line))
    return out


def load_all_suites(directory: str | Path | None = None) -> list[dict[str, Any]]:
    """Load and concatenate every ``*.jsonl`` suite in *directory*."""
    d = Path(directory) if directory is not None else suites_dir()
    out: list[dict[str, Any]] = []
    for p in sorted(d.glob("*.jsonl")):
        out.extend(load_suite(p))
    return out


@dataclass
class EvalReport:
    """Aggregate metrics + per-record detail for an eval run."""

    n: int = 0
    metrics: dict[str, float] = field(default_factory=dict)
    records: list[dict[str, Any]] = field(default_factory=dict)  # type: ignore[assignment]
    violations: list[dict[str, Any]] = field(default_factory=list)
    test_pollution: list[dict[str, Any]] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"Eval over {self.n} records:"]
        for k in sorted(self.metrics):
            lines.append(f"  {k}: {self.metrics[k]:.1%}")
        if self.violations:
            lines.append(f"  hallucination violations: {len(self.violations)}")
        if self.test_pollution:
            lines.append(
                f"  test-file pollution: {len(self.test_pollution)} record(s)"
            )
        return "\n".join(lines)


def _is_test_path(path: str) -> bool:
    """True for a unit-test source path (should never be retrieved)."""
    p = (path or "").replace("\\", "/")
    base = p.rsplit("/", 1)[-1]
    return "/tests/" in p or base.startswith("test_") or base == "conftest.py"


def _accuracy(pairs: list[tuple[Any, Any]]) -> float | None:
    rel = [(g, e) for g, e in pairs if e is not None]
    if not rel:
        return None
    return sum(1 for g, e in rel if g == e) / len(rel)


def evaluate(
    records: list[dict[str, Any]],
    *,
    retriever: Any | None = None,
    registry: Any | None = None,
    k: int = 10,
) -> EvalReport:
    """Run the suite and return an :class:`EvalReport`.

    *retriever* / *registry* default to the project's built ones.
    """
    from pycsamt.agents._workflows import classify_workflow
    from pycsamt.agents.router import IntentRouter
    from pycsamt.assistant.rag.config import infer_workflow

    # Strip code-request phrasing so the *subject* workflow is scored
    # (matches what the assistant's code path does).
    _CODE_PHRASES = (
        "generate code for", "write code for", "code for", "script for",
        "generate code", "write code", "python script",
        "write a script", "create notebook", "notebook for",
        "give me code", "show me code",
    )

    def _subject_workflow(q: str) -> str | None:
        # Registry (handles phrasing: "occam2d inversion" -> pre_inversion).
        wf = classify_workflow(q, default=None)
        if wf == "code_gen":
            t = " " + q.lower() + " "
            for p in _CODE_PHRASES:
                t = t.replace(p, " ")
            wf = classify_workflow(t, default=None)
            if wf == "code_gen":
                wf = None
        # Fallback to symbol-name keywords ("StaticShiftAgent",
        # "estimate_ss_ama") that the phrase registry doesn't cover.
        if wf is None:
            wf = infer_workflow(q)
        return wf

    if retriever is None:
        from pycsamt.assistant.rag.retriever import (
            build_retriever,
        )
        retriever = build_retriever()
    if registry is None:
        try:
            from pycsamt.assistant.tools.project_registry import (
                ProjectRegistry,
            )
            registry = ProjectRegistry.from_default()
        except Exception:  # noqa: BLE001
            registry = None

    router = IntentRouter()  # offline

    intent_pairs: list[tuple[Any, Any]] = []
    wf_pairs: list[tuple[Any, Any]] = []
    line_pairs: list[tuple[Any, Any]] = []
    recalls: list[float] = []
    per_records: list[dict[str, Any]] = []
    violations: list[dict[str, Any]] = []
    test_pollution: list[dict[str, Any]] = []
    wf_in_topk: list[bool] = []       # retrieval surfaced the right area
    nonempty: list[bool] = []         # retrieval returned anything

    for rec in records:
        q = rec["query"]
        got_intent = router.route(q).intent
        got_wf = _subject_workflow(q)
        got_line = (
            registry.find_line_in_text(q) if registry else None
        )

        intent_pairs.append((got_intent, rec.get("expected_intent")))
        wf_pairs.append((got_wf, rec.get("expected_workflow")))
        line_pairs.append((got_line, rec.get("expected_line")))

        ctx = retriever.search(q, k=k)
        got_syms = {c.symbol for c in ctx.chunks if c.symbol}
        exp_syms = set(rec.get("expected_symbols", []))
        recall = (
            len(exp_syms & got_syms) / len(exp_syms)
            if exp_syms else None
        )
        if recall is not None:
            recalls.append(recall)

        # ── retrieval-quality guards ────────────────────────────────────
        nonempty.append(bool(ctx.chunks))
        polluted = [
            c.source_path for c in ctx.chunks if _is_test_path(c.source_path)
        ]
        if polluted:
            test_pollution.append({"query": q, "paths": polluted})
        # Retrieval target: an explicit retrieval expectation (adversarial
        # paraphrases the keyword classifier can't label) or, failing that,
        # the classifier's expected_workflow.
        exp_wf = (
            rec.get("expected_retrieval_workflow")
            or rec.get("expected_workflow")
        )
        retrieved_wfs = {c.workflow for c in ctx.chunks if c.workflow}
        if exp_wf:
            wf_in_topk.append(exp_wf in retrieved_wfs)

        blob = "\n".join(c.text for c in ctx.chunks)
        bad = [
            s for s in rec.get("must_not_contain", []) if s in blob
        ]
        if bad:
            violations.append({"query": q, "found": bad})

        per_records.append(
            {
                "query": q,
                "got_intent": got_intent,
                "got_workflow": got_wf,
                "got_line": got_line,
                "symbol_recall": recall,
                "missing_symbols": sorted(exp_syms - got_syms),
                "retrieved_workflows": sorted(w for w in retrieved_wfs if w),
                "test_polluted": bool(polluted),
            }
        )

    metrics: dict[str, float] = {}
    for name, pairs in (
        ("intent_accuracy", intent_pairs),
        ("workflow_accuracy", wf_pairs),
        ("line_accuracy", line_pairs),
    ):
        acc = _accuracy(pairs)
        if acc is not None:
            metrics[name] = acc
    if recalls:
        metrics["symbol_recall"] = sum(recalls) / len(recalls)
    if records:
        metrics["no_test_in_topk"] = 1.0 - len(test_pollution) / len(records)
        metrics["retrieval_nonempty"] = sum(nonempty) / len(records)
    if wf_in_topk:
        metrics["workflow_in_topk"] = sum(wf_in_topk) / len(wf_in_topk)

    report = EvalReport(
        n=len(records), metrics=metrics, violations=violations,
        test_pollution=test_pollution,
    )
    report.records = per_records
    return report

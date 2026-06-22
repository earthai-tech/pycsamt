# -*- coding: utf-8 -*-
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

__all__ = ["EvalReport", "load_suite", "evaluate", "suites_dir"]


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


@dataclass
class EvalReport:
    """Aggregate metrics + per-record detail for an eval run."""

    n: int = 0
    metrics: dict[str, float] = field(default_factory=dict)
    records: list[dict[str, Any]] = field(default_factory=dict)  # type: ignore[assignment]
    violations: list[dict[str, Any]] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"Eval over {self.n} records:"]
        for k in sorted(self.metrics):
            lines.append(f"  {k}: {self.metrics[k]:.1%}")
        if self.violations:
            lines.append(f"  hallucination violations: {len(self.violations)}")
        return "\n".join(lines)


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
    from pycsamt.agents.router import IntentRouter
    from pycsamt.assistant.rag.config import infer_workflow

    if retriever is None:
        from pycsamt.assistant.rag.retriever import build_retriever
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

    for rec in records:
        q = rec["query"]
        got_intent = router.route(q).intent
        got_wf = infer_workflow(q)
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

    report = EvalReport(
        n=len(records), metrics=metrics, violations=violations
    )
    report.records = per_records
    return report

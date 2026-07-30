# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.assistant.tools.workflow_tools
======================================

Thin, callable wrappers that *execute* pyCSAMT workflows — the "tool"
half of the assistant (RAG retrieves context, tools do the work).

Each tool resolves a survey-line name (via the project registry) or a
filesystem path, then runs the matching agent from
:mod:`pycsamt.agents` and returns a compact, JSON-friendly result.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

__all__ = [
    "resolve_target",
    "run_workflow",
    "run_static_shift",
    "run_qc",
    "run_phase_analysis",
]

# workflow id → agent class name in pycsamt.agents
_WORKFLOW_AGENTS: dict[str, str] = {
    "static_shift": "StaticShiftAgent",
    "qc": "DataQCAgent",
    "phase_analysis": "PhaseAnalysisAgent",
    "denoise": "DenoisingAgent",
    "tipper": "TipperAnalysisAgent",
    "rotation": "TensorRotationAgent",
}


def resolve_target(
    line_or_path: str, *, registry: Any | None = None
) -> dict[str, Any]:
    """Resolve a survey-line name or a path to a concrete EDI location.

    Returns ``{"kind": "line"|"path", "path", "line", "exists", ...}``.
    """
    # try the project registry first (line name)
    if registry is None:
        try:
            from .project_registry import ProjectRegistry

            registry = ProjectRegistry.from_default()
        except Exception:  # noqa: BLE001
            registry = None
    if registry is not None:
        try:
            if line_or_path in registry:
                info = registry.resolve_line(line_or_path)
                info["kind"] = "line"
                return info
        except Exception:  # noqa: BLE001
            pass
    # otherwise treat as a path
    p = Path(line_or_path)
    return {
        "kind": "path",
        "line": None,
        "path": str(p),
        "edi_dir": str(p),
        "exists": p.exists(),
    }


def run_workflow(
    workflow: str,
    line_or_path: str,
    *,
    output_dir: str | None = None,
    registry: Any | None = None,
    history: Any | None = None,
    **agent_kwargs: Any,
) -> dict[str, Any]:
    """Resolve *line_or_path* and run *workflow* via its agent.

    Returns a compact result dict (status, summary, data keys, target).
    When *history* (a ``WorkflowHistory``) is given, the run is appended
    to it for observability. Raises ``ValueError`` for an unknown
    workflow and ``FileNotFoundError`` when the resolved EDI directory
    does not exist.
    """
    if workflow not in _WORKFLOW_AGENTS:
        raise ValueError(
            f"Unsupported workflow {workflow!r}. "
            f"Known: {sorted(_WORKFLOW_AGENTS)}"
        )
    target = resolve_target(line_or_path, registry=registry)
    edi_dir = target.get("edi_dir") or target.get("path")
    if not target.get("exists"):
        raise FileNotFoundError(f"EDI directory not found: {edi_dir!r}")

    out = output_dir or target.get("output_root") or "pycsamt_workflow_output"

    import pycsamt.agents as A

    agent_cls = getattr(A, _WORKFLOW_AGENTS[workflow])
    agent = agent_cls(**agent_kwargs)
    result = agent.execute({"path": edi_dir, "output_dir": out})

    out_result = {
        "workflow": workflow,
        "line": target.get("line"),
        "path": edi_dir,
        "output_dir": out,
        "status": result.status,
        "summary": result.summary,
        "data_keys": sorted((result.data or {}).keys()),
        "n_figures": len((result.data or {}).get("figures", {}) or {}),
    }
    if history is not None:
        try:
            from ..memory.workflow_history import WorkflowRun

            history.record(
                WorkflowRun(
                    workflow=workflow,
                    status=result.status,
                    line=target.get("line"),
                    path=edi_dir,
                    output_dir=out,
                    summary=result.summary,
                    n_figures=out_result["n_figures"],
                )
            )
        except Exception:  # noqa: BLE001 — tracing must never break a run
            pass
    return out_result


def run_static_shift(line_or_path: str, **kw: Any) -> dict[str, Any]:
    """Run static-shift correction on a line/path."""
    return run_workflow("static_shift", line_or_path, **kw)


def run_qc(line_or_path: str, **kw: Any) -> dict[str, Any]:
    """Run quality control on a line/path."""
    return run_workflow("qc", line_or_path, **kw)


def run_phase_analysis(line_or_path: str, **kw: Any) -> dict[str, Any]:
    """Run phase-tensor / dimensionality analysis on a line/path."""
    return run_workflow("phase_analysis", line_or_path, **kw)

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.agents._workflow_plan
==============================

:class:`WorkflowPlan` — validated intermediate object produced by
:class:`~pycsamt.agents.context.ContextInputAgent` and consumed by
:class:`~pycsamt.agents.orchestrator.WorkflowOrchestratorAgent`.

Every natural-language request is converted to a ``WorkflowPlan``
before execution.  This gives the orchestrator a single validated
artefact to act on and gives the provenance trail a canonical
representation of user intent.
"""
from __future__ import annotations

import json
import warnings as _warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["WorkflowPlan", "validate_workflow_plan"]

# Canonical set of workflow identifiers recognised by the orchestrator
VALID_WORKFLOWS = frozenset({
    "qc",
    "static_shift",
    "phase_analysis",
    "forward",
    "pre_inversion",
    "inversion_prep",
    "inversion_eval",
    "interpretation",
    "report",
    "full",
    "ai_inversion",
    "inv1d",
    "inv2d",
    "inv3d",
    "ensemble_inversion",
    "joint_inversion",
    "modem",
    "occam2d",
    "tipper",
    "sensitivity",
    "rotation",
    "freq_decimation",
    "batch",
    "comparison",
    "full_ai_workflow",
})


@dataclass
class WorkflowPlan:
    r"""
    Validated intermediate representation of a workflow request.

    Produced by :class:`~pycsamt.agents.context.ContextInputAgent`
    and consumed by
    :class:`~pycsamt.agents.orchestrator.WorkflowOrchestratorAgent`.

    Parameters
    ----------
    request : str
        Original natural-language request.
    workflow_type : str
        One of the recognised workflow identifiers
        (see :data:`VALID_WORKFLOWS`).
    data_path : str
        Path to EDI files or survey directory.
    output_dir : str
        Destination for outputs.
    parameters : dict
        Parsed parameters (period_range, component, etc.).
    risk_flags : list of str
        Warnings or risks the user should review.
    requires_human_review : bool
        True when ambiguous parameters were inferred.
    expected_outputs : list of str
        Files / artefacts the workflow should produce.
    provider : str
        LLM provider used for parsing, or ``'offline'``.
    """

    request: str
    workflow_type: str
    data_path: str = ""
    output_dir: str = ""
    parameters: dict[str, Any] = field(default_factory=dict)
    risk_flags: list[str] = field(default_factory=list)
    requires_human_review: bool = False
    expected_outputs: list[str] = field(default_factory=list)
    provider: str = "offline"

    # ── validation ────────────────────────────────────────────

    def is_valid(self) -> bool:
        """Return True when the plan passes basic validation."""
        return (
            bool(self.request)
            and self.workflow_type in VALID_WORKFLOWS
        )

    def validation_errors(self) -> list[str]:
        """Return a list of validation error messages."""
        errs: list[str] = []
        if not self.request:
            errs.append("'request' is empty.")
        if self.workflow_type not in VALID_WORKFLOWS:
            errs.append(
                f"Unknown workflow_type {self.workflow_type!r}. "
                f"Valid: {sorted(VALID_WORKFLOWS)}"
            )
        if not self.data_path:
            errs.append(
                "No data_path extracted from the request."
            )
        elif not Path(self.data_path).exists():
            errs.append(
                f"data_path does not exist: {self.data_path!r}"
            )
        return errs

    # ── serialisation ──────────────────────────────────────────

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable dict."""
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Serialise to a JSON string."""
        return json.dumps(self.to_dict(), indent=indent)

    def save(self, path: str | Path) -> None:
        """Write the plan to *path* as JSON."""
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def from_config(
        cls,
        config: dict[str, Any],
        request: str = "",
        provider: str = "offline",
    ) -> WorkflowPlan:
        """
        Build a :class:`WorkflowPlan` from a config dict
        as returned by :class:`ContextInputAgent`.
        """
        wf = config.get("workflow", "qc")

        params: dict[str, Any] = {}
        for k in (
            "period_range", "component", "station",
            "inversion_code", "depth_max_km", "n_periods",
            "verbose",
        ):
            if k in config:
                params[k] = config[k]

        # identify risk flags
        flags: list[str] = []
        if not config.get("data_path"):
            flags.append(
                "No data path found — user must supply."
            )
        if wf not in VALID_WORKFLOWS:
            flags.append(
                f"Unrecognised workflow {wf!r}."
            )

        expected = _expected_outputs(wf)
        needs_review = bool(flags) or not config.get("data_path")

        return cls(
            request=request,
            workflow_type=wf,
            data_path=config.get("data_path", ""),
            output_dir=config.get("output_dir", ""),
            parameters=params,
            risk_flags=flags,
            requires_human_review=needs_review,
            expected_outputs=expected,
            provider=provider,
        )


# ── helpers ────────────────────────────────────────────────────────────────────

def _expected_outputs(workflow_type: str) -> list[str]:
    """Return a list of expected output artefacts for *workflow_type*."""
    base = ["agent_trace.json", "workflow_plan.json"]
    extras: dict[str, list[str]] = {
        "qc":                ["qc_confidence.png",
                              "qc_report.md"],
        "phase_analysis":    ["pt_psection.png",
                              "strike_analysis.png",
                              "phase_tensor_report.md"],
        "ai_inversion":      ["inv1d_section.png",
                              "predictions.npz",
                              "ai_inversion_report.md"],
        "inv1d":             ["inv1d_section.png",
                              "predictions.npz",
                              "ai_inversion_report.md"],
        "inv2d":             ["inv2d_section.png",
                              "pred_section.npz",
                              "inv2d_report.md"],
        "inv3d":             ["inv3d_volume.png",
                              "inv3d_report.md"],
        "ensemble_inversion":["ensemble_section.png",
                              "uncertainty_bands.npz",
                              "ensemble_report.md"],
        "joint_inversion":   ["joint_section.png",
                              "joint_report.md"],
        "pre_inversion":     ["occam2d_data.dat",
                              "occam2d_mesh.txt",
                              "occam2d_startup.txt"],
        "modem":             ["ModEM_Data.dat"],
        "full":              ["qc_confidence.png",
                              "pt_psection.png",
                              "inv1d_section.png",
                              "full_survey_report.md",
                              "workflow_script.py"],
        "full_ai_workflow":  ["qc_confidence.png",
                              "pt_psection.png",
                              "inv1d_section.png",
                              "inv2d_section.png",
                              "workflow_script.py",
                              "full_ai_report.md"],
    }
    return base + extras.get(workflow_type, ["report.md"])


def validate_workflow_plan(
    plan: WorkflowPlan,
    *,
    raise_on_error: bool = False,
) -> list[str]:
    """
    Validate *plan* and optionally raise on errors.

    Parameters
    ----------
    plan : WorkflowPlan
    raise_on_error : bool
        When True, raises :exc:`ValueError` if errors exist.

    Returns
    -------
    errors : list of str
        Empty when the plan is valid.
    """
    errs = plan.validation_errors()
    if errs:
        if raise_on_error:
            raise ValueError(
                "WorkflowPlan validation failed:\n"
                + "\n".join(f"  - {e}" for e in errs)
            )
        for e in errs:
            _warnings.warn(e, stacklevel=2)
    return errs

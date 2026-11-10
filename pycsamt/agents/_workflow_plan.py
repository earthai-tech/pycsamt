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
import math
from numbers import Integral, Real
import warnings as _warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

__all__ = ["WorkflowPlan", "validate_workflow_plan"]

# Canonical set of workflow identifiers recognised by the orchestrator
VALID_WORKFLOWS = frozenset(
    {
        "conventional_inversion",
        "field_inversion",
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
        "mare2dem",
        "pinn_inversion",
        "hybrid_inversion",
        "code_gen",
        "denoise",
    }
)


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
    grounding_citations : list of str
        Source paths retrieved by RAG to ground this plan, when
        retrieval was available (empty offline or on a cold index).
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
    grounding_citations: list[str] = field(default_factory=list)

    # ── validation ────────────────────────────────────────────

    def is_valid(self) -> bool:
        """Return True only when all declared plan checks pass."""
        return not self.validation_errors()

    def validation_errors(self) -> list[str]:
        """Return a list of validation error messages."""
        errs: list[str] = []
        if not isinstance(self.request, str) or not self.request.strip():
            errs.append("'request' is empty.")
        if not isinstance(self.workflow_type, str) or self.workflow_type not in VALID_WORKFLOWS:
            errs.append(
                f"Unknown workflow_type {self.workflow_type!r}. "
                f"Valid: {sorted(VALID_WORKFLOWS)}"
            )
        if not self.data_path:
            errs.append("No data_path extracted from the request.")
        else:
            try:
                if not Path(self.data_path).exists():
                    errs.append(f"data_path does not exist: {self.data_path!r}")
            except (TypeError, ValueError, OSError):
                errs.append("data_path must be a valid file or directory path.")
        if not isinstance(self.requires_human_review, bool):
            errs.append("requires_human_review must be a boolean.")
        errs.extend(parameter_errors(self.parameters))
        return errs

    def approval_digest(self) -> str:
        """Bind explicit review approval to this plan and current input bytes."""
        from ._checkpoint_identity import input_identities, digest

        return digest({"plan": self.to_dict(), "inputs": input_identities(self.to_dict())})

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
        citations: list[str] | None = None,
    ) -> WorkflowPlan:
        """
        Build a :class:`WorkflowPlan` from a config dict
        as returned by :class:`ContextInputAgent`.
        """
        wf = config.get("workflow", "qc")

        # Preserve the full scientific configuration, including per-step and
        # solver settings, so approval cannot outlive a parameter change.
        excluded = {"workflow", "data_path", "output_dir", "risk_flags",
                    "requires_human_review", "approved_plan_hash"}
        params = {k: v for k, v in config.items() if k not in excluded}

        # identify risk flags
        flags: list[str] = list(config.get("risk_flags") or [])
        if not config.get("data_path"):
            flags.append("No data path found — user must supply.")
        if not isinstance(wf, str) or wf not in VALID_WORKFLOWS:
            flags.append(f"Unrecognised workflow {wf!r}.")

        expected = _expected_outputs(wf) if isinstance(wf, str) else []
        needs_review = wf in ("conventional_inversion", "field_inversion") or config.get("requires_human_review", False) or bool(flags) or not config.get("data_path")

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
            grounding_citations=list(citations) if citations else [],
        )


# ── helpers ────────────────────────────────────────────────────────────────────


def parameter_errors(parameters: Any, prefix: str = "parameters") -> list[str]:
    """Validate supported core settings without repairing scientific intent.

    This is a bounded contract, not validation of every solver-specific option.
    Nested per-step settings receive the same numeric period/depth/count checks.
    """
    if not isinstance(parameters, dict):
        return [f"{prefix} must be a mapping."]
    errors = []
    period = parameters.get("period_range")
    if period is not None:
        valid = isinstance(period, (list, tuple)) and len(period) == 2
        if valid:
            valid = all(isinstance(x, Real) and not isinstance(x, bool) and math.isfinite(x) for x in period)
        if not valid or not (0 < period[0] < period[1]):
            errors.append(f"{prefix}.period_range must contain two finite positive bounds in increasing order.")
    component = parameters.get("component")
    if component is not None and component not in ("xy", "yx", "xx", "yy", "all", "off_diagonal"):
        errors.append(f"{prefix}.component is unsupported: {component!r}.")
    backend = parameters.get("inversion_code")
    if backend is not None and backend not in ("occam2d", "modem", "mare2dem"):
        errors.append(f"{prefix}.inversion_code is unsupported: {backend!r}.")
    for name in ("depth_max_km", "depth_max"):
        value = parameters.get(name)
        if value is not None and (not isinstance(value, Real) or isinstance(value, bool) or not math.isfinite(value) or value <= 0):
            errors.append(f"{prefix}.{name} must be finite and positive.")
    for name in ("n_periods", "n_layers", "epochs", "n_train_samples"):
        value = parameters.get(name)
        if value is not None and (not isinstance(value, Integral) or isinstance(value, bool) or value < 1):
            errors.append(f"{prefix}.{name} must be a positive integer.")
    for name in ("step_params", "pinn_params", "hybrid_params", "ai_inv_params"):
        nested = parameters.get(name)
        if nested is not None:
            if name == "step_params" and isinstance(nested, dict):
                for step, settings in nested.items():
                    errors.extend(parameter_errors(settings, f"{prefix}.{name}.{step}"))
            else:
                errors.extend(parameter_errors(nested, f"{prefix}.{name}"))
    return errors


def _expected_outputs(workflow_type: str) -> list[str]:
    """Return a list of expected output artefacts for *workflow_type*."""
    base = [
        "workflow_plan.json",
        "agent_trace.json",
        "environment.json",
        "output_manifest.json",
    ]
    extras: dict[str, list[str]] = {
        "qc": ["qc_confidence.png", "qc_report.md"],
        "phase_analysis": [
            "pt_psection.png",
            "strike_analysis.png",
            "phase_tensor_report.md",
        ],
        "ai_inversion": [
            "inv1d_section.png",
            "predictions.npz",
            "ai_inversion_report.md",
        ],
        "inv1d": [
            "inv1d_section.png",
            "predictions.npz",
            "ai_inversion_report.md",
        ],
        "inv2d": ["inv2d_section.png", "pred_section.npz", "inv2d_report.md"],
        "inv3d": ["inv3d_volume.png", "inv3d_report.md"],
        "ensemble_inversion": [
            "ensemble_section.png",
            "uncertainty_bands.npz",
            "ensemble_report.md",
        ],
        "joint_inversion": ["joint_section.png", "joint_report.md"],
        "pre_inversion": [
            "occam2d_data.dat",
            "occam2d_mesh.txt",
            "occam2d_startup.txt",
        ],
        "modem": ["ModEM_Data.dat"],
        "field_inversion": ["static_shift_report.md", "report.md"],
        "full": [
            "qc_confidence.png",
            "pt_psection.png",
            "inv1d_section.png",
            "full_survey_report.md",
            "workflow_script.py",
        ],
        "full_ai_workflow": [
            "qc_confidence.png",
            "pt_psection.png",
            "inv1d_section.png",
            "inv2d_section.png",
            "workflow_script.py",
            "full_ai_report.md",
        ],
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

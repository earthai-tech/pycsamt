"""
pycsamt.agents.coordinator
==========================

:class:`AgentCoordinator` orchestrates multi-agent workflows:

* Registers agents by name.
* Executes ordered workflow steps with dependency tracking.
* Checkpoints each step to disk so workflows can be resumed after failure.
* Provides a ``dry_run`` / ``preview`` mode that estimates cost and prints
  the execution plan without running any agents.
* Aggregates per-step costs into a workflow-level total.
"""

from __future__ import annotations

import json
import logging
import pickle
import time
from pathlib import Path
from typing import Any

from ._base import AgentResult, BaseAgent
from ._pricing import format_cost

logger = logging.getLogger(__name__)


# ── workflow step descriptor ──────────────────────────────────────────────────

class WorkflowStep:
    """One step in an :class:`AgentCoordinator` workflow.

    Parameters
    ----------
    name : str
        Unique identifier used for checkpointing and logging.
    agent : BaseAgent
        The agent instance that runs this step.
    input_fn : callable or None
        ``input_fn(prev_results)`` → dict fed to ``agent.execute()``.
        When ``None`` the coordinator passes the raw workflow config.
    description : str
        Human-readable description shown in the dry-run preview.
    required : bool
        When ``False`` a failure skips the step rather than aborting.
    """

    def __init__(
        self,
        name: str,
        agent: BaseAgent,
        *,
        input_fn=None,
        description: str = "",
        required: bool = True,
    ) -> None:
        self.name        = name
        self.agent       = agent
        self.input_fn    = input_fn
        self.description = description or f"Run {type(agent).__name__}"
        self.required    = required


# ── coordinator ───────────────────────────────────────────────────────────────

class AgentCoordinator:
    """Orchestrate a sequence of pycsamt agents as a named workflow.

    Parameters
    ----------
    workflow_name : str
        Name used for log messages and checkpoint directory.
    checkpoint_dir : str or Path or None
        Where to save step checkpoints.  Defaults to
        ``./pycsamt_agent_checkpoints/<workflow_name>/``.
    verbose : bool
        Print step-level progress to stdout.

    Examples
    --------
    Build and run a QC workflow::

        from pycsamt.agents import AgentCoordinator, MTLoaderAgent, DataQCAgent

        coord = AgentCoordinator("mt_qc")
        coord.add_step("load",  MTLoaderAgent(...))
        coord.add_step("qc",    DataQCAgent(...),
                       input_fn=lambda r: {"sites": r["load"]["sites"]})

        result = coord.execute({"path": "/data/EDIs"})
        print(result.summary)

    Dry-run preview::

        result = coord.preview({"path": "/data/EDIs"})
        print(result["plan"])
    """

    def __init__(
        self,
        workflow_name: str = "pycsamt_workflow",
        *,
        checkpoint_dir: str | Path | None = None,
        verbose: bool = True,
    ) -> None:
        self.workflow_name = workflow_name
        self.verbose       = verbose
        self._steps: list[WorkflowStep]    = []
        self._agents: dict[str, BaseAgent] = {}

        if checkpoint_dir is None:
            checkpoint_dir = Path(f"pycsamt_agent_checkpoints/{workflow_name}")
        self._ckpt_dir = Path(checkpoint_dir)

        self._log = logging.getLogger(f"pycsamt.agents.coordinator.{workflow_name}")

    # ── step registration ─────────────────────────────────────────────────────

    def add_step(
        self,
        name: str,
        agent: BaseAgent,
        *,
        input_fn=None,
        description: str = "",
        required: bool = True,
    ) -> "AgentCoordinator":
        """Append a workflow step.  Returns ``self`` for chaining.

        Parameters
        ----------
        name : str
        agent : BaseAgent
        input_fn : callable(prev_results_dict) → dict, optional
            Maps accumulated results to this step's input dict.
            ``prev_results_dict`` keys are earlier step names;
            each value is the step's :class:`AgentResult`.
        description : str
        required : bool

        Returns
        -------
        AgentCoordinator
            ``self`` for method chaining.
        """
        if name in {s.name for s in self._steps}:
            raise ValueError(f"Step {name!r} already registered.")
        self._steps.append(
            WorkflowStep(
                name, agent,
                input_fn=input_fn,
                description=description,
                required=required,
            )
        )
        self._agents[name] = agent
        return self

    # ── preview ───────────────────────────────────────────────────────────────

    def preview(self, config: dict[str, Any]) -> AgentResult:
        """Return an execution plan and estimated cost without running anything.

        Parameters
        ----------
        config : dict
            The workflow configuration that would be passed to ``execute()``.

        Returns
        -------
        AgentResult
            ``data["plan"]`` contains the formatted plan string.
            ``data["steps"]`` is a list of step metadata dicts.
        """
        lines = [
            f"Workflow: {self.workflow_name}",
            f"Steps   : {len(self._steps)}",
            f"Config  : {json.dumps(config, indent=2, default=str)}",
            "",
            "─" * 60,
        ]
        step_meta = []
        for i, step in enumerate(self._steps, 1):
            agent = step.agent
            llm_str = (
                f"{agent.llm_provider}/{agent.model}"
                if agent.api_key else "no-LLM"
            )
            req_str = "" if step.required else " [optional]"
            lines.append(
                f"  {i:2d}. [{step.name}]{req_str}\n"
                f"       Agent  : {type(agent).__name__}\n"
                f"       LLM    : {llm_str}\n"
                f"       Action : {step.description}"
            )
            step_meta.append({
                "step": i, "name": step.name,
                "agent": type(agent).__name__,
                "llm": llm_str,
                "description": step.description,
                "required": step.required,
            })
        lines.append("─" * 60)
        plan = "\n".join(lines)
        if self.verbose:
            print(plan)
        return AgentResult(
            status="success",
            summary=f"Workflow preview: {len(self._steps)} steps.",
            data={"plan": plan, "steps": step_meta, "config": config},
        )

    # ── execute ───────────────────────────────────────────────────────────────

    def execute(
        self,
        config: dict[str, Any],
        *,
        dry_run: bool = False,
        resume: bool = False,
    ) -> AgentResult:
        """Run all workflow steps sequentially.

        Parameters
        ----------
        config : dict
            Top-level workflow configuration forwarded to every step's
            ``input_fn`` (or directly to ``agent.execute()`` when no
            ``input_fn`` is set).
        dry_run : bool
            When ``True`` return a preview without executing any agents.
        resume : bool
            When ``True`` skip steps whose checkpoint already exists on disk.

        Returns
        -------
        AgentResult
            ``data`` contains one key per step name with its :class:`AgentResult`.
        """
        if dry_run:
            return self.preview(config)

        self._ckpt_dir.mkdir(parents=True, exist_ok=True)
        t_workflow = time.time()
        results: dict[str, AgentResult] = {}
        total_cost = 0.0
        all_warnings: list[str] = []

        self._log.info("Starting workflow %r (%d steps)", self.workflow_name, len(self._steps))

        for step in self._steps:
            # ── resume: skip completed checkpoints ───────────────────────────
            ckpt_path = self._ckpt_dir / f"{step.name}.pkl"
            if resume and ckpt_path.exists():
                try:
                    results[step.name] = self._load_checkpoint(step.name)
                    self._print(f"  ↩  [{step.name}] resumed from checkpoint")
                    continue
                except Exception as exc:
                    self._log.warning("Failed to load checkpoint for %r: %s", step.name, exc)

            # ── build input ──────────────────────────────────────────────────
            try:
                if step.input_fn is not None:
                    step_input = step.input_fn(results)
                else:
                    step_input = config
            except Exception as exc:
                msg = f"input_fn for step {step.name!r} raised: {exc}"
                self._log.error(msg)
                if step.required:
                    return AgentResult.failed(
                        msg, elapsed=time.time() - t_workflow
                    )
                all_warnings.append(msg)
                continue

            # ── run agent ────────────────────────────────────────────────────
            self._print(f"  ▶  [{step.name}] {step.description}")
            try:
                result = step.agent.execute(step_input)
            except Exception as exc:
                result = AgentResult.failed(
                    str(exc),
                    hint="Check agent configuration and input data.",
                    elapsed=0.0,
                )

            results[step.name] = result
            total_cost        += result.cost_estimate_usd
            all_warnings      += result.warnings

            # ── checkpoint ───────────────────────────────────────────────────
            self._save_checkpoint(step.name, result)

            # ── status check ─────────────────────────────────────────────────
            status_icon = {"success": "✔", "failed": "✘", "needs_review": "⚠"}.get(
                result.status, "?"
            )
            cost_str = format_cost(result.cost_estimate_usd)
            self._print(
                f"  {status_icon}  [{step.name}] {result.summary} "
                f"({result.elapsed_seconds:.1f}s, {cost_str})"
            )

            if result.status == "failed" and step.required:
                self._log.error("Required step %r failed; aborting workflow.", step.name)
                return AgentResult(
                    status="failed",
                    summary=f"Workflow aborted at required step {step.name!r}: {result.error}",
                    data=results,
                    warnings=all_warnings,
                    error=result.error,
                    error_fix_hint=result.error_fix_hint,
                    elapsed_seconds=time.time() - t_workflow,
                    cost_estimate_usd=total_cost,
                )

        elapsed = time.time() - t_workflow
        n_ok  = sum(1 for r in results.values() if r.status == "success")
        n_tot = len(results)
        summary = (
            f"Workflow {self.workflow_name!r} complete: "
            f"{n_ok}/{n_tot} steps succeeded "
            f"in {elapsed:.1f}s ({format_cost(total_cost)})."
        )
        self._print(f"\n  {'━'*52}\n  {summary}\n  {'━'*52}")
        self._save_workflow_state(results, summary)

        return AgentResult(
            status="success" if n_ok == n_tot else "needs_review",
            summary=summary,
            data=results,
            warnings=all_warnings,
            elapsed_seconds=elapsed,
            cost_estimate_usd=total_cost,
        )

    # ── checkpoint helpers ────────────────────────────────────────────────────

    def _save_checkpoint(self, step_name: str, result: AgentResult) -> None:
        path = self._ckpt_dir / f"{step_name}.pkl"
        try:
            with open(path, "wb") as f:
                pickle.dump(result, f)
            # JSON sidecar for human inspection
            meta = {
                "step": step_name,
                "status": result.status,
                "summary": result.summary,
                "elapsed_seconds": result.elapsed_seconds,
                "cost_estimate_usd": result.cost_estimate_usd,
                "warnings": result.warnings,
                "error": result.error,
            }
            with open(path.with_suffix(".json"), "w") as f:
                json.dump(meta, f, indent=2)
        except Exception as exc:
            self._log.warning("Could not save checkpoint for %r: %s", step_name, exc)

    def _load_checkpoint(self, step_name: str) -> AgentResult:
        path = self._ckpt_dir / f"{step_name}.pkl"
        with open(path, "rb") as f:
            return pickle.load(f)

    def _save_workflow_state(
        self, results: dict[str, AgentResult], summary: str
    ) -> None:
        state = {
            "workflow": self.workflow_name,
            "summary": summary,
            "steps": {
                name: {
                    "status": r.status,
                    "elapsed_seconds": r.elapsed_seconds,
                    "cost_estimate_usd": r.cost_estimate_usd,
                    "warnings": r.warnings,
                }
                for name, r in results.items()
            },
        }
        path = self._ckpt_dir / "workflow_state.json"
        try:
            with open(path, "w") as f:
                json.dump(state, f, indent=2)
        except Exception as exc:
            self._log.warning("Could not save workflow state: %s", exc)

    # ── helpers ───────────────────────────────────────────────────────────────

    def _print(self, msg: str) -> None:
        if self.verbose:
            print(msg)

    def reset_checkpoints(self) -> None:
        """Delete all saved checkpoints for this workflow."""
        if self._ckpt_dir.exists():
            import shutil
            shutil.rmtree(self._ckpt_dir)
            self._log.info("Checkpoints cleared for workflow %r", self.workflow_name)

    def __repr__(self) -> str:
        return (
            f"AgentCoordinator(name={self.workflow_name!r}, "
            f"steps={len(self._steps)})"
        )


__all__ = ["AgentCoordinator", "WorkflowStep"]

"""
pycsamt.agents.batch_survey
==============================

:class:`BatchSurveyAgent` — Run the full processing pipeline over multiple
MT profiles in parallel.

Each profile (EDI directory or file path) is processed independently through
a configurable agent chain.  Results are merged into a survey-level summary
DataFrame and a multi-panel overview figure.

Parallelism is provided by ``joblib`` (soft dependency).  When ``n_jobs=1``
or joblib is unavailable, processing is sequential.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in large-scale MT survey processing and quality control.
Given a batch processing result, write 3-4 sentences that:
1. Report how many profiles succeeded and how many failed.
2. Identify which profiles have the worst QC scores or highest RMS.
3. Highlight any systematic issues across profiles (e.g. dead band at same frequency).
4. Recommend which profiles need re-processing or manual inspection.
Reply in plain English.
"""

_SUPPORTED_WORKFLOWS = {
    "qc":             ["MTLoaderAgent", "DataQCAgent", "StaticShiftAgent"],
    "ai_inversion":   ["MTLoaderAgent", "DataQCAgent", "DenoisingAgent",
                       "AIInversionAgent"],
    "phase_analysis": ["MTLoaderAgent", "DataQCAgent", "StaticShiftAgent",
                       "PhaseAnalysisAgent"],
    "sensitivity":    ["MTLoaderAgent", "SensitivityAgent"],
    "tipper":         ["MTLoaderAgent", "TipperAnalysisAgent"],
}


class BatchSurveyAgent(BaseAgent):
    """Process multiple MT profiles through a shared agent chain.

    Parameters
    ----------
    api_key, model, llm_provider : str
    workflow : str
        Pre-defined workflow key: ``'qc'``, ``'ai_inversion'``,
        ``'phase_analysis'``, ``'sensitivity'``, ``'tipper'``.
        Determines which agents are chained for each profile.
    n_jobs : int
        Parallel workers. -1 = all CPU cores. 1 = sequential (default).

    Input keys
    ----------
    ``profiles`` : dict {name: path} or list[str]
        Profile names to paths, or a list of paths (names auto-assigned).
    ``workflow`` : str, optional — override constructor default
    ``n_jobs`` : int, optional
    ``output_dir`` : str, optional
    Any additional keys are forwarded to each agent's ``execute()`` call.

    Output data keys
    ----------------
    ``profile_results``   dict {name: AgentResult}
    ``summary_table``     pandas.DataFrame — per-profile metrics
    ``n_success``         int
    ``n_failed``          int
    ``figures``           dict
    ``figure_paths``      dict
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        workflow: str = "qc",
        n_jobs: int = 1,
    ) -> None:
        super().__init__(
            "BatchSurveyAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
        )
        self.workflow = workflow
        self.n_jobs   = n_jobs

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        profiles_raw = input_data.get("profiles")
        if profiles_raw is None:
            return AgentResult.failed(
                "No 'profiles' key provided.",
                hint="Pass profiles={'L01': '/data/L01', 'L02': '/data/L02'} "
                     "or a list of paths.",
                elapsed=time.time() - t0,
            )

        # normalise to dict
        if isinstance(profiles_raw, (list, tuple)):
            profiles = {f"profile_{i:02d}": str(p)
                        for i, p in enumerate(profiles_raw)}
        else:
            profiles = {str(k): str(v) for k, v in profiles_raw.items()}

        if not profiles:
            return AgentResult.failed("Empty profiles dict.", elapsed=time.time() - t0)

        workflow  = str(input_data.get("workflow", self.workflow))
        n_jobs    = int(input_data.get("n_jobs",   self.n_jobs))
        output_dir = input_data.get("output_dir")

        agent_names = _SUPPORTED_WORKFLOWS.get(workflow)
        if agent_names is None:
            warnings.append(
                f"Unknown workflow {workflow!r}; using 'qc'. "
                f"Supported: {list(_SUPPORTED_WORKFLOWS)}."
            )
            workflow    = "qc"
            agent_names = _SUPPORTED_WORKFLOWS["qc"]

        # ── build per-profile extra kwargs ────────────────────────────────
        extra = {k: v for k, v in input_data.items()
                 if k not in ("profiles", "workflow", "n_jobs", "output_dir")}

        # ── run profiles ──────────────────────────────────────────────────
        profile_results: dict[str, AgentResult] = {}

        def _run_one(name: str, path: str) -> tuple[str, AgentResult]:
            import os
            prof_out = os.path.join(output_dir, name) if output_dir else None
            inp = {"path": path, "output_dir": prof_out}
            inp.update(extra)
            try:
                return name, _run_pipeline(
                    agent_names, inp,
                    api_key=self.api_key,
                    model=self.model,
                    llm_provider=self.llm_provider,
                )
            except Exception as exc:
                return name, AgentResult.failed(str(exc), elapsed=0.0)

        try:
            from joblib import Parallel, delayed
            pairs = Parallel(n_jobs=n_jobs, prefer="threads")(
                delayed(_run_one)(name, path)
                for name, path in profiles.items()
            )
        except ImportError:
            warnings.append("joblib not installed — processing sequentially.")
            pairs = [_run_one(name, path) for name, path in profiles.items()]

        for name, result in pairs:
            profile_results[name] = result

        # ── summary table ─────────────────────────────────────────────────
        n_success = sum(1 for r in profile_results.values() if r.status != "failed")
        n_failed  = len(profile_results) - n_success

        rows = []
        for name, r in profile_results.items():
            row = {
                "profile": name,
                "status":  r.status,
                "elapsed": round(r.elapsed_seconds, 2),
                "cost":    round(r.cost_estimate_usd, 6),
                "n_warnings": len(r.warnings),
            }
            # pull common metrics
            for key in ("rms_global", "rms", "n_stations", "n_flagged",
                        "mean_doi_km"):
                val = r.get(key)
                if val is not None:
                    row[key] = val
            rows.append(row)

        try:
            import pandas as pd
            summary_table = pd.DataFrame(rows)
        except ImportError:
            summary_table = rows

        # ── figure: status overview ────────────────────────────────────────
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        try:
            fig = _plot_batch_summary(profile_results)
            if fig is not None:
                figures["batch_summary"] = fig
                p = self._save_figure(fig, output_dir, "batch_survey_summary",
                                      warnings_list=warnings)
                if p:
                    fig_paths["batch_summary"] = p
        except Exception as exc:
            warnings.append(f"Batch summary figure: {exc}")

        # ── LLM interpretation ────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and profile_results:
            worst = max(profile_results,
                        key=lambda n: len(profile_results[n].warnings), default="?")
            prompt = (
                f"Batch survey processing ({workflow}):\n"
                f"  Profiles: {len(profiles)}\n"
                f"  Succeeded: {n_success} | Failed: {n_failed}\n"
                f"  Profile with most warnings: {worst}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Summarise the batch results and flag profiles needing attention."
            )
            interp = self.query_llm(prompt, max_tokens=180)

        elapsed = time.time() - t0
        return AgentResult(
            status="success" if n_success > 0 else "failed",
            summary=(
                f"Batch '{workflow}': {n_success}/{len(profiles)} profiles succeeded, "
                f"{n_failed} failed. {len(figures)} figures."
            ),
            data={
                "profile_results": profile_results,
                "summary_table":   summary_table,
                "n_success":       n_success,
                "n_failed":        n_failed,
                "workflow":        workflow,
                "figures":         figures,
                "figure_paths":    fig_paths,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── private helpers ───────────────────────────────────────────────────────────

def _run_pipeline(
    agent_names: list[str],
    input_data: dict[str, Any],
    *,
    api_key: str | None,
    model: str | None,
    llm_provider: str,
) -> AgentResult:
    """Run agent_names sequentially; pass sites between steps."""
    import importlib

    common = dict(api_key=api_key, model=model, llm_provider=llm_provider)
    current_input = dict(input_data)
    last_result: AgentResult | None = None

    for agent_name in agent_names:
        try:
            mod  = importlib.import_module(f".{_agent_module(agent_name)}",
                                           package="pycsamt.agents")
            cls  = getattr(mod, agent_name)
            agent = cls(**common)
            result = agent.execute(current_input)
            last_result = result

            # chain Sites forward
            for key in ("denoised_sites", "corrected_sites", "sites"):
                val = result.get(key)
                if val is not None:
                    current_input["sites"] = val
                    break

        except Exception as exc:
            return AgentResult.failed(
                f"{agent_name} failed: {exc}", elapsed=0.0,
            )

    return last_result or AgentResult.failed(
        "No agents ran.", elapsed=0.0,
    )


def _agent_module(name: str) -> str:
    _MAP = {
        "MTLoaderAgent":        "loader",
        "DataQCAgent":          "qc",
        "StaticShiftAgent":     "static_shift",
        "DenoisingAgent":       "denoising",
        "AIInversionAgent":     "ai_inversion",
        "PhaseAnalysisAgent":   "phase_analysis",
        "SensitivityAgent":     "sensitivity",
        "TipperAnalysisAgent":  "tipper_analysis",
    }
    return _MAP.get(name, name.lower())


def _plot_batch_summary(profile_results: dict[str, Any]) -> Any:
    """Horizontal bar chart: elapsed time per profile, coloured by status."""
    import matplotlib.pyplot as plt

    names   = list(profile_results.keys())
    elapsed = [profile_results[n].elapsed_seconds for n in names]
    colors  = ["#27ae60" if profile_results[n].status != "failed"
               else "#e74c3c" for n in names]

    fig, ax = plt.subplots(figsize=(6, max(3, len(names) * 0.35)))
    ax.barh(range(len(names)), elapsed, color=colors, alpha=0.85, edgecolor="white")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Elapsed time (s)", fontsize=9)
    ax.set_title("Batch survey — processing time per profile\n"
                 "(green = success, red = failed)", fontsize=9, fontweight="bold")
    ax.tick_params(labelsize=8)
    ax.invert_yaxis()
    fig.tight_layout()
    return fig


__all__ = ["BatchSurveyAgent"]

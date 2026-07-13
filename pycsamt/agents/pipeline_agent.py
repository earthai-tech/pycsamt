"""
pycsamt.agents.pipeline_agent
==============================

:class:`PipelineAgent` — LLM-assisted MT processing pipeline selection and
interpretation.

Bridges :mod:`pycsamt.pipeline` with the agent framework:

* Accepts a natural-language ``"request"`` key → LLM recommends preset/steps
* Builds and runs the :class:`~pycsamt.pipeline.Pipeline`
* LLM interprets the :class:`~pycsamt.pipeline.PipelineResult` as a narrative

Available presets: ``"basic_qc"``, ``"noise_reduction"``, ``"full_processing"``,
``"tensor_analysis"``, ``"dimensionality_filter"``, ``"publication_ready"``.
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult, BaseAgent

# ── LLM system messages ───────────────────────────────────────────────────────

_SYSTEM_RECOMMEND = """\
You are an expert in magnetotelluric (MT/AMT/CSAMT) data processing pipelines.
Given a user description of their dataset and processing goals, recommend the most
appropriate pipeline configuration as a JSON object with exactly these fields:
  "preset"          : one of "basic_qc", "noise_reduction", "full_processing",
                      "tensor_analysis", "dimensionality_filter",
                      "publication_ready", or null
  "steps"           : list of step codes if preset is null
                      (e.g. ["NR001","FREQ002","SS001"]) — empty list otherwise
  "param_overrides" : dict of {step_code: {param: value}} for parameter tweaks,
                      e.g. {"NR001": {"mains_hz": 60}} for a 60 Hz grid
  "rationale"       : one sentence explaining the recommendation

Reply with ONLY valid JSON — no markdown fences, no text outside the JSON.
"""

_SYSTEM_INTERPRET = """\
You are an expert MT data processing specialist reviewing pipeline execution results.
Given a processing summary, write 3–4 sentences that:
1. Describe the overall data quality change (stations retained, step errors if any).
2. Highlight which steps were most impactful or took the most time.
3. Comment on whether the stated processing objectives appear to have been met.
4. Recommend a concrete follow-up action (further correction, QC plot review,
   or readiness for inversion).
Reply in plain English. No bullet points or markdown.
"""


class PipelineAgent(BaseAgent):
    """LLM-assisted MT processing pipeline selection and interpretation.

    Two modes of operation:

    **Guided mode** — pass ``"request"`` in *input_data*.  The LLM recommends
    a preset and any parameter overrides, then the pipeline is built and run
    automatically.

    **Direct mode** — pass ``"preset"`` (name string) *or* ``"steps"`` (list of
    step codes) in *input_data*.  No pre-run LLM call is made; the pipeline is
    built directly from those instructions.

    In both modes the post-run LLM call interprets the
    :class:`~pycsamt.pipeline.PipelineResult` as a narrative.

    Parameters
    ----------
    api_key, model, llm_provider : str
        Standard LLM configuration inherited from :class:`~._base.BaseAgent`.
    preset : str or None
        Default preset name used when *input_data* contains neither
        ``"preset"`` nor ``"steps"`` nor ``"request"``.
        Defaults to ``"basic_qc"`` (safe first-pass).
    param_overrides : dict or None
        Default parameter overrides applied on top of any preset or step list.
        Format: ``{step_code: {param: value}}``.

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str
        Raw MT/AMT sites to process.
    ``request`` : str, optional
        Natural-language description of dataset and goals.  Triggers a
        pre-run LLM call that recommends preset + parameter overrides.
    ``preset`` : str, optional
        Named preset — overrides constructor default.
    ``steps`` : list of str, optional
        Explicit ordered list of step codes.  Ignored when ``"preset"`` is set.
    ``param_overrides`` : dict, optional
        Per-step parameter overrides — merged on top of constructor defaults.
    ``output_dir`` : str or None, optional
        Root directory for pipeline output files (EDI, plots, YAML config).

    Output data keys
    ----------------
    ``sites_out``        Processed Sites — ready for downstream agents
    ``pipeline_result``  :class:`~pycsamt.pipeline.PipelineResult` object
    ``preset_used``      Name of the preset that was run (or ``"custom"``)
    ``steps_run``        List of step code strings that were executed
    ``n_sites_in``       Number of input stations
    ``n_sites_out``      Number of stations after processing
    ``n_errors``         Number of steps that raised an error
    ``recommendation``   Dict returned by LLM pre-run call (or ``None``)

    Examples
    --------
    Guided mode::

        agent  = PipelineAgent()
        result = agent.execute({
            "sites": sites,
            "request": "50 Hz grid noise, possible static shift, Occam2D target",
            "output_dir": "willy_pipeline/",
        })
        processed = result["sites_out"]
        print(result.llm_interpretation)

    Direct mode::

        agent  = PipelineAgent(preset="full_processing")
        result = agent.execute({
            "sites": sites,
            "param_overrides": {"NR001": {"mains_hz": 60}},
        })

    Chain with :class:`~.occam2d_agent.Occam2DAgent` via
    :class:`~.coordinator.AgentCoordinator`::

        from pycsamt.agents import AgentCoordinator, PipelineAgent, Occam2DAgent

        coord = AgentCoordinator("willy_full")
        coord.add_step("pipeline", PipelineAgent(preset="full_processing"),
                       input_fn=lambda r: {"sites": r["load"].data["sites"]})
        coord.add_step("invert",   Occam2DAgent(),
                       input_fn=lambda r: {"sites": r["pipeline"].data["sites_out"]})
    """

    SYSTEM_PROMPT = _SYSTEM_INTERPRET

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        preset: str | None = "basic_qc",
        param_overrides: dict[str, dict] | None = None,
    ) -> None:
        super().__init__(
            "PipelineAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.preset = preset
        self.param_overrides = param_overrides or {}

    # ── execute ───────────────────────────────────────────────────────────────

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings_list: list[str] = []

        # ── resolve sites ─────────────────────────────────────────────────────
        from ..emtools._core import ensure_sites

        sites_raw = input_data.get("sites") or input_data.get("path")
        if sites_raw is None:
            return AgentResult.failed(
                "No 'sites' or 'path' provided.", elapsed=time.time() - t0
            )
        try:
            sites = ensure_sites(sites_raw, verbose=0)
        except Exception as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        n_sites_in = _count(sites)

        # ── resolve effective config from input_data ──────────────────────────
        # Priority: explicit "preset" > explicit "steps" > constructor default
        if "preset" in input_data:
            effective_preset = input_data["preset"]
            effective_steps: list[str] = []
        elif "steps" in input_data:
            effective_preset = None
            effective_steps = list(input_data["steps"])
        else:
            effective_preset = self.preset
            effective_steps = []

        effective_overrides: dict[str, dict] = {
            **self.param_overrides,
            **input_data.get("param_overrides", {}),
        }
        output_dir = input_data.get("output_dir")

        # ── guided mode: LLM pre-run recommendation ───────────────────────────
        recommendation: dict | None = None
        request = input_data.get("request")

        if request:
            rec_text = self.query_llm(
                str(request),
                system_message=_SYSTEM_RECOMMEND,
                temperature=0.1,
                max_tokens=512,
            )
            if rec_text:
                parsed = self.extract_json(rec_text)
                if isinstance(parsed, dict):
                    recommendation = parsed
                    effective_preset = (
                        parsed.get("preset") or effective_preset
                    )
                    effective_steps = parsed.get("steps") or effective_steps
                    merged_ov = parsed.get("param_overrides") or {}
                    effective_overrides = {**effective_overrides, **merged_ov}
                else:
                    warnings_list.append(
                        "LLM recommendation could not be parsed as JSON; "
                        "falling back to default preset."
                    )

        # ── build pipeline ────────────────────────────────────────────────────
        from ..pipeline import (
            PYCSAMT_PIPE,
            Pipeline,
            Step,
        )

        try:
            if effective_preset:
                pipe = Pipeline.from_preset(effective_preset)
                preset_used = effective_preset
            elif effective_steps:
                pipe = Pipeline(
                    [(code, Step(code)) for code in effective_steps],
                    name="custom",
                )
                preset_used = "custom"
            else:
                # last resort: basic_qc is always safe
                pipe = Pipeline.from_preset("basic_qc")
                preset_used = "basic_qc"
                warnings_list.append(
                    "No preset or steps specified; defaulting to 'basic_qc'."
                )
        except Exception as exc:
            return AgentResult.failed(
                f"Could not build pipeline: {exc}", elapsed=time.time() - t0
            )

        # ── apply param overrides ─────────────────────────────────────────────
        if effective_overrides:
            for label, step in list(pipe._steps):
                code = step.spec.code
                if code in effective_overrides:
                    try:
                        merged_params = {
                            **step.params,
                            **effective_overrides[code],
                        }
                        pipe.replace(label, Step(code, **merged_params))
                    except Exception as exc:
                        warnings_list.append(
                            f"Could not apply override for step {code!r}: {exc}"
                        )

        # ── run pipeline (suppress progress to stdout) ────────────────────────
        try:
            with PYCSAMT_PIPE.context(show_progress=False):
                pipe_result = pipe.run(
                    sites,
                    outdir=output_dir,
                    save_plots=output_dir is not None,
                    save_edis=output_dir is not None,
                    save_report=output_dir is not None,
                )
        except Exception as exc:
            return AgentResult.failed(
                f"Pipeline run failed: {exc}",
                hint="Check that all step codes are valid and Sites is non-empty.",
                elapsed=time.time() - t0,
            )

        n_sites_out = _count(pipe_result.sites_out)
        steps_run = [sr.step_code for sr in pipe_result.step_results]
        n_errors = pipe_result.n_errors

        # ── build summary text for LLM ────────────────────────────────────────
        step_lines = "\n".join(
            f"  [{sr.step_code}] {sr.step_label}  "
            f"{'OK' if sr.ok else 'ERR'}  "
            f"{sr.n_sites_in}→{sr.n_sites_out} sites  "
            f"{sr.elapsed_sec:.2f}s"
            + (f"  error={sr.error}" if not sr.ok else "")
            for sr in pipe_result.step_results
        )
        request_line = f"User request : {request}" if request else ""
        summary_text = (
            f"Pipeline     : {preset_used}\n"
            f"Sites        : {n_sites_in} in → {n_sites_out} out\n"
            f"Steps        : {len(steps_run)} "
            f"({len(steps_run) - n_errors} OK, {n_errors} errors)\n"
            f"Elapsed      : {pipe_result.elapsed_sec:.2f}s\n"
        )
        if request_line:
            summary_text = request_line + "\n" + summary_text
        summary_text += "\nStep details:\n" + step_lines

        # ── LLM post-run interpretation ───────────────────────────────────────
        llm_interp = self.query_llm(
            summary_text,
            system_message=_SYSTEM_INTERPRET,
            temperature=0.2,
            max_tokens=512,
        )

        # ── surface pipeline warnings ─────────────────────────────────────────
        for sr in pipe_result.step_results:
            if not sr.ok:
                warnings_list.append(
                    f"Step {sr.step_code!r} ({sr.step_label}) raised "
                    f"{type(sr.error).__name__}: {sr.error}"
                )

        elapsed = time.time() - t0
        status = "success" if n_errors == 0 else "needs_review"
        brief = (
            f"Pipeline '{preset_used}' ran {len(steps_run)} steps in "
            f"{pipe_result.elapsed_sec:.1f}s: "
            f"{n_sites_in}→{n_sites_out} sites, {n_errors} errors."
        )

        return AgentResult(
            status=status,
            summary=brief,
            data={
                "sites_out": pipe_result.sites_out,
                "pipeline_result": pipe_result,
                "preset_used": preset_used,
                "steps_run": steps_run,
                "n_sites_in": n_sites_in,
                "n_sites_out": n_sites_out,
                "n_errors": n_errors,
                "recommendation": recommendation,
            },
            warnings=warnings_list,
            llm_interpretation=llm_interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── helpers ───────────────────────────────────────────────────────────────────


def _count(sites: Any) -> int:
    """Return number of stations in *sites*, tolerating any Sites-like object."""
    for attr in ("n_sites", "n_stations", "__len__"):
        try:
            val = getattr(sites, attr, None)
            if callable(val):
                return int(val())
            if val is not None:
                return int(val)
        except Exception:
            pass
    try:
        return len(sites)
    except Exception:
        return 0

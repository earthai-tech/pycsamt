"""
pycsamt.agents.context
======================

:class:`ContextInputAgent` — Natural-language → pycsamt workflow config.

Accepts a free-text request such as::

    "Load the EDI files in /data/WILLY_DATA/L22PLT, remove static shift,
     run phase tensor analysis for periods 1e-4 to 1 s, and save a report
     to /output/survey_report/"

and returns a validated :class:`~pycsamt.agents._base.AgentResult` whose
``data["config"]`` is a structured dict the :class:`AgentCoordinator` can
execute directly.

When no LLM API key is available the agent falls back to a robust set of
regex patterns that cover the most common request forms.
"""

from __future__ import annotations

import logging
import os
import re
import time
from pathlib import Path
from typing import Any

from ._base import AgentResult, BaseAgent

logger = logging.getLogger(__name__)


# ── constants ─────────────────────────────────────────────────────────────────

# Workflow descriptions, keyword table, and alias map all live in the single
# source of truth: pycsamt.agents._workflows (shared with the orchestrator).
from ._workflows import (  # noqa: E402
    WORKFLOW_DESCRIPTIONS as _KNOWN_WORKFLOWS,
    WORKFLOW_ALIASES as _WORKFLOW_ALIASES,
    normalise_workflow as _normalise_workflow,
)

_KNOWN_INVERSION_CODES = {"occam2d", "modem", "nlcg", "zonal", "smooth2d"}

_COMPONENTS = {"xy", "yx", "xx", "yy", "all", "off_diagonal"}

# ── system prompt ─────────────────────────────────────────────────────────────

_SYSTEM_PROMPT = """You are an expert MT/AMT/CSAMT workflow configuration \
interpreter for pycsamt v2.

Given a natural-language processing request, extract a structured JSON \
configuration dictionary with the following schema (include only keys \
that are clearly mentioned or can be reasonably inferred):

{
  "workflow":       string — one of:
                    qc, static_shift, phase_analysis, forward,
                    inversion_prep, pre_inversion, inversion_eval,
                    interpretation, report, full,
                    ai_inversion, inv1d, inv2d, inv3d,
                    ensemble_inversion, joint_inversion,
                    modem, occam2d,
                    tipper, sensitivity, rotation,
                    freq_decimation, batch, comparison,
                    full_ai_workflow,
  "data_path":      string — absolute or relative path to EDI
                    file(s) or directory,
  "output_dir":     string — where to write results / figures,
  "period_range":   [T_min_seconds, T_max_seconds],
  "component":      string — "xy"|"yx"|"all"|"off_diagonal",
  "station":        string or null,
  "inversion_code": string — "occam2d"|"modem"|null,
  "depth_max_km":   float or null,
  "n_periods":      int or null,
  "verbose":        bool
}

Rules:
- Choose "ai_inversion" for CNN / 1-D neural-network / deep-learning
  inversion requests.
- Choose "inv2d" for U-Net / 2-D neural / profile AI inversion.
- Choose "inv3d" for GCN / graph-convolutional / 3-D AI inversion.
- Choose "ensemble_inversion" for ensemble / uncertainty / Bayesian.
- Choose "joint_inversion" for joint / multi-modal / TEM+MT.
- Choose "full_ai_workflow" when both AI inversion and full pipeline
  are requested together.
- If a period range is given in frequency (Hz), convert to period
  (s = 1/f).
- Preserve the full absolute path exactly as given.
- Return ONLY the JSON object — no markdown fences, no prose.
"""


# ── agent ─────────────────────────────────────────────────────────────────────

class ContextInputAgent(BaseAgent):
    """Parse a natural-language MT workflow request into a structured config.

    Parameters
    ----------
    api_key : str or None
        LLM API key.  When ``None`` the regex fallback is used exclusively.
    model, llm_provider : str
        Passed to :class:`~pycsamt.agents._base.BaseAgent`.

    Examples
    --------
    With an API key::

        agent = ContextInputAgent(api_key="sk-ant-…")
        result = agent.execute({
            "request": "Load EDIs from /data/L22PLT, QC them, "
                       "period range 1e-4 to 1 s, save to /out/qc/"
        })
        cfg = result["config"]
        # cfg["workflow"]    == "qc"
        # cfg["data_path"]   == "/data/L22PLT"
        # cfg["period_range"] == [0.0001, 1.0]

    Without an API key (regex fallback)::

        agent = ContextInputAgent()          # no key → regex mode
        result = agent.execute({"request": "…"})
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
    ) -> None:
        super().__init__(
            "ContextInputAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()

        # accept "request" or "text" or bare string
        request: str = (
            input_data.get("request")
            or input_data.get("text")
            or input_data.get("prompt")
            or ""
        )
        if not request:
            return AgentResult.failed(
                "No request text found in input_data. "
                "Pass {'request': '<your request>'}.",
                elapsed=time.time() - t0,
            )

        # Regex runs first to anchor workflow type.
        # LLM enriches params but must not override
        # a regex-confident (non-qc) workflow.
        regex_cfg = _regex_extract(request)
        regex_wf = regex_cfg.get("workflow", "")

        config: dict[str, Any] | None = None
        llm_raw: str | None = None

        if self.api_key:
            llm_raw = self.query_llm(request)
            if llm_raw:
                config = self.extract_json(llm_raw)
                if config:
                    self._log.debug(
                        "LLM extraction succeeded."
                    )
                    # Regex wins when it found a specific
                    # workflow (not the qc default). LLM
                    # over-generalises (e.g. maps
                    # freq_decimation -> ai_inversion).
                    if regex_wf and regex_wf != "qc":
                        config["workflow"] = regex_wf

        if config is None:
            config = regex_cfg
            if not self.api_key:
                self._log.debug(
                    "No API key -- using regex."
                )
            else:
                self._log.warning(
                    "LLM failed; falling back to regex."
                )

        # ── normalise & validate ──────────────────────────────────────────────
        config = _normalise_config(config, request)
        warnings = _validate_config(config)

        # ── build validated WorkflowPlan ──────────────────────────
        from ._workflow_plan import WorkflowPlan
        plan = WorkflowPlan.from_config(
            config,
            request=request,
            provider=(
                self.llm_provider if self.api_key else "offline"
            ),
        )

        # ── optional LLM summary ──────────────────────────────────
        interpretation: str | None = None
        if self.api_key and config.get("data_path"):
            interp_prompt = (
                "Briefly summarise in 2 sentences what the "
                "following pycsamt workflow configuration "
                f"will do:\n{config}"
            )
            interpretation = self.query_llm(
                interp_prompt,
                system_message=(
                    "You are a concise MT data processing "
                    "assistant. Reply in plain English, no "
                    "bullet points."
                ),
                max_tokens=120,
            )

        elapsed = time.time() - t0
        summary = (
            f"Config extracted: "
            f"workflow={config.get('workflow', '?')!r}, "
            f"path={config.get('data_path', '?')!r}."
        )

        return AgentResult(
            status=(
                "success"
                if config.get("data_path")
                else "needs_review"
            ),
            summary=summary,
            data={
                "config":        config,
                "workflow_plan": plan,
                "raw_request":   request,
                "llm_raw":       llm_raw,
            },
            warnings=warnings,
            llm_interpretation=interpretation,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


# ── regex extraction ──────────────────────────────────────────────────────────

def _regex_extract(text: str) -> dict[str, Any]:
    """Extract a workflow config from *text* using regex patterns.

    This covers the most common natural-language patterns without LLM.
    """
    cfg: dict[str, Any] = {}
    t = text.lower()

    # ── data path ─────────────────────────────────────────────────────────────
    # Path must start with / or ~ (absolute / home-relative) to avoid matching
    # plain words like "EDI" that appear before the actual path.
    path_patterns = [
        # keyword immediately followed by an absolute/home path
        r'(?:load(?:ing)?|read(?:ing)?|from|on|in|at|path[:\s]+)'
        r'\s+["\']?([/~][\w/\\\-\.]+)["\']?',
        # quoted absolute path (any extension)
        r'["\']([/~][\w/\\\-\.]+)["\']',
        # bare absolute path — must NOT be preceded by a
        # word char or slash (prevents matching inner components
        # of a path like the "/willy" in "/data/willy")
        r'(?<![/\w])([/~][\w/\\\-\.]{5,})',
    ]
    for pat in path_patterns:
        m = re.search(pat, text, re.IGNORECASE)
        if m:
            candidate = m.group(1).rstrip(".,;)")
            if len(candidate) > 4:
                # os.path.expanduser (not Path) so the user's original
                # separators survive on Windows — Path() would rewrite
                # a POSIX "/data/x" to "\data\x".
                cfg["data_path"] = os.path.expanduser(candidate)
                break

    # ── output directory ──────────────────────────────────────────────────────
    m = re.search(
        r'(?:save|output|write|report)\s+(?:to\s+)?["\']?([/~\w][\w/\\\-\.]+)["\']?',
        text, re.IGNORECASE,
    )
    if m:
        cfg["output_dir"] = os.path.expanduser(m.group(1))

    # ── workflow ──────────────────────────────────────────────────────────
    # Single source of truth for the keyword table lives in
    # pycsamt.agents._workflows (shared with the orchestrator).
    # Leave the key unset when nothing matches so _normalise_config
    # applies the qc default.
    from ._workflows import classify_workflow
    _wf = classify_workflow(text)
    if _wf:
        cfg["workflow"] = _wf

    # ── period range ──────────────────────────────────────────────────────────
    # e.g. "1e-4 to 1 s", "0.0001 – 1.0 s", "T_min=0.001, T_max=10"
    period_pat = re.compile(
        r'(?:period[s]?\s*(?:range|from|between|:)?\s*)'
        r'([\d\.e\+\-]+)\s*(?:to|–|-|,)\s*([\d\.e\+\-]+)',
        re.IGNORECASE,
    )
    m = period_pat.search(text)
    if m:
        try:
            cfg["period_range"] = [float(m.group(1)), float(m.group(2))]
        except ValueError:
            pass

    # frequency → period conversion
    freq_pat = re.compile(
        r'(?:freq(?:uency)?\s*(?:range|from|:)?\s*)'
        r'([\d\.e\+\-]+)\s*(?:to|–|-|,)\s*([\d\.e\+\-]+)\s*(?:hz)?',
        re.IGNORECASE,
    )
    m = freq_pat.search(text)
    if m and "period_range" not in cfg:
        try:
            f1, f2 = float(m.group(1)), float(m.group(2))
            if f1 > 0 and f2 > 0:
                cfg["period_range"] = sorted([1.0 / f1, 1.0 / f2])
        except ValueError:
            pass

    # ── component ─────────────────────────────────────────────────────────────
    comp_pat = re.compile(
        r'\b(xy|yx|xx|yy|off[-_]diagonal|all(?:\s+components?)?)\b',
        re.IGNORECASE,
    )
    m = comp_pat.search(text)
    if m:
        raw = m.group(1).lower().replace(" ", "_").replace("-", "_")
        cfg["component"] = "off_diagonal" if "diagonal" in raw else raw.split("_")[0]

    # ── inversion code ────────────────────────────────────────────────────────
    for code in _KNOWN_INVERSION_CODES:
        if code in t:
            cfg["inversion_code"] = code
            break

    # ── depth ─────────────────────────────────────────────────────────────────
    m = re.search(r'depth[_\s]*(?:max)?\s*[=:]\s*([\d\.]+)\s*(km|m)?',
                  text, re.IGNORECASE)
    if m:
        val = float(m.group(1))
        unit = (m.group(2) or "km").lower()
        cfg["depth_max_km"] = val if unit == "km" else val / 1000.0

    # ── station ───────────────────────────────────────────────────────────────
    m = re.search(
        r'station[s]?\s*[=:\s]+([A-Za-z0-9_\-]+)',
        text, re.IGNORECASE,
    )
    if m:
        cfg["station"] = m.group(1)

    return cfg


# ── normalise & validate ──────────────────────────────────────────────────────
# _WORKFLOW_ALIASES is imported at module top from pycsamt.agents._workflows
# (single source of truth, shared with the orchestrator).


def _normalise_config(cfg: dict[str, Any], original_text: str) -> dict[str, Any]:
    """Fill defaults, normalise aliases, clean up types."""
    # workflow alias normalisation (shared registry)
    cfg["workflow"] = _normalise_workflow(cfg.get("workflow", ""))

    # keep unrecognised workflows as-is — the orchestrator will
    # report a proper error rather than silently running QC
    if not cfg["workflow"]:
        cfg["workflow"] = "qc"

    # path normalisation — os.path.expanduser preserves the user's
    # separators (Path() would rewrite POSIX paths to backslashes on
    # Windows); only "~" is expanded.
    if "data_path" in cfg:
        cfg["data_path"] = os.path.expanduser(str(cfg["data_path"]))

    if "output_dir" in cfg:
        cfg["output_dir"] = os.path.expanduser(str(cfg["output_dir"]))
    else:
        cfg["output_dir"] = str(Path("pycsamt_agent_output").resolve())

    # period_range
    if "period_range" in cfg:
        pr = cfg["period_range"]
        if isinstance(pr, (list, tuple)) and len(pr) == 2:
            lo, hi = float(pr[0]), float(pr[1])
            cfg["period_range"] = [min(lo, hi), max(lo, hi)]

    # component normalisation
    comp = str(cfg.get("component", "xy")).lower()
    if comp not in _COMPONENTS:
        comp = "xy"
    cfg["component"] = comp

    # boolean verbose
    cfg.setdefault("verbose", True)
    cfg.setdefault("station", None)
    cfg.setdefault("inversion_code", None)

    return cfg


def _validate_config(cfg: dict[str, Any]) -> list[str]:
    """Return a list of warning messages for issues found in *cfg*."""
    warnings: list[str] = []

    data_path = cfg.get("data_path")
    if not data_path:
        warnings.append(
            "No data path found in the request. "
            "Please supply a path to EDI files or a directory."
        )
    elif not Path(data_path).exists():
        warnings.append(
            f"Data path does not exist on disk: {data_path!r}. "
            "Check the path and try again."
        )

    pr = cfg.get("period_range")
    if pr is not None:
        lo, hi = pr
        if lo <= 0:
            warnings.append(f"period_range lower bound must be > 0; got {lo}.")
        if hi <= lo:
            warnings.append(
                f"period_range upper bound {hi} ≤ lower bound {lo}."
            )

    return warnings


__all__ = ["ContextInputAgent"]

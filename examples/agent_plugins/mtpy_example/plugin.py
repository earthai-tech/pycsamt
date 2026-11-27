"""MTpy example plugin for the pyCSAMT agent loop.

Three workflows - ``read_edi_plot_rhophi``, ``phase_tensor`` and
``strike`` - are classified by keyword, validated (EDI data present,
plausible frequency band and parameters) and recorded with the
built-in :class:`pycsamt.agents._provenance.RunRecord`.

MTpy is an optional dependency: nothing here imports it at module
level, so the plugin registers, classifies and validates without it.
Only the step callables in :data:`STEP_FUNCTIONS` touch MTpy, and they
raise :class:`ImportError` with a clear message when it is missing.

Examples
--------
>>> from mtpy_example import make_plugin
>>> p = make_plugin()
>>> p.registry.classify("plot apparent resistivity and phase")
'read_edi_plot_rhophi'
>>> p.registry.family("strike")
'phase_tensor'
"""
from __future__ import annotations

import importlib.util
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from pycsamt.agents.plugins import DomainPlugin, RecordWriter

__all__ = [
    "MTPY_WORKFLOWS",
    "STEP_FUNCTIONS",
    "MTpyRegistry",
    "MTpyPlanValidator",
    "make_plugin",
    "mtpy_available",
]

#: Workflow id -> one-line description.
MTPY_WORKFLOWS: dict[str, str] = {
    "read_edi_plot_rhophi": (
        "Read EDI files and plot apparent resistivity and phase."
    ),
    "phase_tensor": "Compute and map the phase tensor per period.",
    "strike": "Estimate geoelectric strike from the phase tensor.",
}

# Ordered: first match wins, so the most specific ids come first.
_KEYWORDS: dict[str, tuple[str, ...]] = {
    "strike": ("strike", "geoelectric direction"),
    "phase_tensor": ("phase tensor", "phase-tensor", "ellipse"),
    "read_edi_plot_rhophi": (
        "apparent resistivity", "rho", "res/phase", "plot edi",
        "read edi",
    ),
}
_ALIASES = {"rhophi": "read_edi_plot_rhophi", "pt": "phase_tensor"}
_FAMILY = {"strike": "phase_tensor"}

#: Accepted frequency range of an MT/AMT/CSAMT band, in Hz.
FREQ_LIMITS_HZ = (1e-5, 1e5)


def mtpy_available() -> bool:
    """Return ``True`` when the ``mtpy`` package can be imported."""
    return importlib.util.find_spec("mtpy") is not None


def _require_mtpy():
    if not mtpy_available():
        raise ImportError(
            "This step needs MTpy (pip install mtpy-v2); the plugin "
            "itself works without it."
        )
    import mtpy

    return mtpy


def _step_read(results: dict) -> dict:
    mtpy = _require_mtpy()
    return {"collection": mtpy.MTCollection, "inputs": results}


#: Step name -> callable run by a real executor (needs MTpy).
STEP_FUNCTIONS = {"read": _step_read}


class MTpyRegistry:
    """:class:`~pycsamt.agents.plugins.WorkflowRegistry` for MTpy."""

    name = "mtpy_example"

    @property
    def descriptions(self) -> Mapping[str, str]:
        return MTPY_WORKFLOWS

    @property
    def valid_workflows(self) -> frozenset[str]:
        return frozenset(MTPY_WORKFLOWS)

    @property
    def dispatchable(self) -> frozenset[str]:
        return frozenset(MTPY_WORKFLOWS)

    def classify(
        self, text: str, *, default: str | None = None
    ) -> str | None:
        low = (text or "").lower()
        for wf, words in _KEYWORDS.items():
            if any(w in low for w in words):
                return wf
        return default

    def normalise(self, workflow: str) -> str:
        wf = (workflow or "").strip().lower().replace("-", "_")
        return _ALIASES.get(wf, wf)

    def family(self, workflow: str | None) -> str | None:
        if not workflow:
            return None
        wf = self.normalise(workflow)
        return _FAMILY.get(wf, wf)

    def steps(self, workflow: str) -> Sequence[tuple] | None:
        if workflow not in MTPY_WORKFLOWS:
            return None
        chain = [("read", "MTpyReader", None, "Read EDI files")]
        if workflow == "read_edi_plot_rhophi":
            chain.append(("plot", "MTpyPlotter", "read",
                          "Plot rho_a and phase"))
        else:
            chain.append(("phase_tensor", "MTpyPhaseTensor", "read",
                          "Phase tensor per period"))
            if workflow == "strike":
                chain.append(("strike", "MTpyStrike", "phase_tensor",
                              "Strike from phase tensor"))
        return chain


class MTpyPlanValidator:
    """:class:`~pycsamt.agents.plugins.PlanValidator` for MTpy.

    Works on any plan object with ``workflow_type``, ``data_path``,
    ``parameters`` and ``risk_flags`` attributes, e.g.
    :class:`pycsamt.agents._workflow_plan.WorkflowPlan`.

    Notes
    -----
    Parameter checks, all optional keys of ``plan.parameters``:

    * ``freq_band`` ``(f_min, f_max)`` in Hz, finite, with
      :math:`10^{-5} \\le f_{min} < f_{max} \\le 10^{5}`;
    * ``strike_tolerance_deg`` in ``(0, 90]`` degrees;
    * ``n_periods`` an integer in ``[1, 500]``.
    """

    def __init__(self, registry: MTpyRegistry | None = None) -> None:
        self.registry = registry or MTpyRegistry()

    def errors(self, plan: Any) -> list[str]:
        errs: list[str] = []
        wf = getattr(plan, "workflow_type", "")
        if wf not in self.registry.valid_workflows:
            errs.append(f"Unknown MTpy workflow {wf!r}.")
        errs.extend(_path_errors(getattr(plan, "data_path", "")))
        params = getattr(plan, "parameters", None) or {}
        errs.extend(_param_errors(params))
        return errs

    def needs_review(self, plan: Any) -> bool:
        return bool(getattr(plan, "risk_flags", None))

    def cross_check(
        self, text: str, intent: str, workflow: str | None
    ) -> str | None:
        kw = self.registry.classify(text)
        if kw is None or not workflow:
            return None
        if self.registry.family(kw) == self.registry.family(workflow):
            return None
        return (f"LLM workflow {workflow!r} disagrees with keyword "
                f"route {kw!r}.")

    def route_flags(self, text: str, workflow: str) -> list[str]:
        reason = self.cross_check(text, "workflow", workflow)
        return [reason] if reason else []


def _path_errors(data_path: str) -> list[str]:
    if not data_path:
        return ["No data_path given."]
    p = Path(data_path)
    if not p.exists():
        return [f"data_path does not exist: {data_path!r}"]
    edis = [p] if p.is_file() else list(p.glob("*.edi"))
    edis = [e for e in edis if e.suffix.lower() == ".edi"]
    return [] if edis else [f"No .edi files in {data_path!r}."]


def _param_errors(params: Mapping[str, Any]) -> list[str]:
    errs: list[str] = []
    if "freq_band" in params:
        band = params["freq_band"]
        lo_lim, hi_lim = FREQ_LIMITS_HZ
        try:
            lo, hi = (float(v) for v in band)
        except (TypeError, ValueError):
            errs.append("freq_band must be two numbers (Hz).")
        else:
            if not (math.isfinite(lo) and math.isfinite(hi)):
                errs.append("freq_band must be finite.")
            elif not lo_lim <= lo < hi <= hi_lim:
                errs.append(
                    f"freq_band {band!r} outside "
                    f"{lo_lim:g}-{hi_lim:g} Hz or not increasing."
                )
    if "strike_tolerance_deg" in params:
        tol = params["strike_tolerance_deg"]
        if not (isinstance(tol, (int, float)) and 0 < tol <= 90):
            errs.append("strike_tolerance_deg must be in (0, 90].")
    if "n_periods" in params:
        n = params["n_periods"]
        if not (isinstance(n, int) and not isinstance(n, bool)
                and 1 <= n <= 500):
            errs.append("n_periods must be an integer in [1, 500].")
    return errs


def _record(payload: dict) -> RecordWriter:
    from pycsamt.agents._provenance import RunRecord

    class MTpyRunRecord(RunRecord):
        def corpus_manifest(self) -> Any:
            return None  # this plugin uses no retrieval index

    return MTpyRunRecord(payload)


def make_plugin() -> DomainPlugin:
    """Return the plugin bundle (entry-point target)."""
    reg = MTpyRegistry()
    return DomainPlugin(
        name=reg.name,
        registry=reg,
        validator=MTpyPlanValidator(reg),
        record_factory=_record,
        description="Example MTpy workflows (rho/phase, PT, strike).",
    )

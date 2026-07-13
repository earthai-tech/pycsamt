# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents._package_context
================================
Runtime package-knowledge builder (deep edition).

Extracts docstrings at three depths:

1. Class-level — what the class is for
2. Parameter-level — __init__ Parameters section
3. Method-level — key public methods / properties

The result is split into *topic tiers* so that
PackageQAAgent can inject only the relevant parts
per question (query-adaptive injection), keeping
the LLM context signal-to-noise ratio high.

Tiers
-----
TIER_CORE      Always included; overview + workflows
TIER_AGENTS    Agent classes: params + input/output
TIER_SITES     Sites / EDI data model: methods + access
TIER_EXAMPLES  Usage code examples

PACKAGE_CONTEXT = all tiers joined (full dump for
                  system-prompt pre-loading).
"""

from __future__ import annotations

import inspect
from typing import Any

# ── low-level helpers ──────────────────────────────────────────────────────────


def _full_doc(obj: Any) -> str:
    """Return full docstring, normalized whitespace."""
    return inspect.getdoc(obj) or ""


def _first_para(doc: str, max_chars: int = 280) -> str:
    """Return first paragraph of a docstring."""
    first = (doc or "").split("\n\n")[0].replace("\n", " ")
    return first[:max_chars] + "..." if len(first) > max_chars else first


def _sig(cls: type, max_chars: int = 120) -> str:
    """Return __init__ signature without 'self'."""
    try:
        s = str(inspect.signature(cls.__init__))
        s = s.replace("(self, ", "(").replace("(self)", "()")
        return s[:max_chars]
    except (ValueError, TypeError):
        return "()"


def _section_text(doc: str, section: str) -> str:
    r"""
    Extract a named section from a numpydoc
    docstring (e.g. 'Parameters', 'Input keys').

    Returns empty string when section not found.
    """
    lines = doc.split("\n")
    out: list[str] = []
    in_sec = False
    for i, ln in enumerate(lines):
        if ln.strip() == section:
            # check next line is dashes
            nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
            if set(nxt) <= {"-"}:
                in_sec = True
                continue
        if in_sec:
            # stop at next section header
            if (
                ln.strip()
                and not ln.startswith(" ")
                and i + 1 < len(lines)
                and set(lines[i + 1].strip()) <= {"-"}
            ):
                break
            out.append(ln)
    return "\n".join(out).strip()


def _method_entry(
    name: str,
    method: Any,
    max_chars: int = 200,
) -> str:
    """Format a single method as a reference entry."""
    doc = _first_para(_full_doc(method), max_chars)
    try:
        sig = str(inspect.signature(method))
        sig = sig.replace("(self, ", "(").replace("(self)", "()")[:80]
    except (ValueError, TypeError):
        sig = "(...)"
    return f"  .{name}{sig}\n    {doc}"


def _section(title: str, body: str) -> str:
    bar = "=" * len(title)
    return f"\n{title}\n{bar}\n{body.strip()}\n"


# ── importers ─────────────────────────────────────────────────────────────────


def _import(module_path: str, name: str) -> Any:
    """Import name from module, return None on failure."""
    try:
        import importlib

        mod = importlib.import_module(module_path)
        return getattr(mod, name, None)
    except Exception:
        return None


def _collect_agents() -> dict[str, type]:
    _MAP = {
        "WorkflowOrchestratorAgent": "pycsamt.agents.orchestrator",
        "ContextInputAgent": "pycsamt.agents.context",
        "MTLoaderAgent": "pycsamt.agents",
        "DataQCAgent": "pycsamt.agents",
        "StaticShiftAgent": "pycsamt.agents.static_shift",
        "DenoisingAgent": "pycsamt.agents",
        "PhaseAnalysisAgent": "pycsamt.agents",
        "ForwardModelAgent": "pycsamt.agents",
        "InterpretationAgent": "pycsamt.agents.interpretation",
        "ReportAgent": "pycsamt.agents",
        "AIInversionAgent": "pycsamt.agents",
        "PINNInversionAgent": "pycsamt.agents.pinn_agent",
        "HybridInversionAgent": "pycsamt.agents.hybrid_agent",
        "Inv2DAgent": "pycsamt.agents",
        "Inv3DAgent": "pycsamt.agents",
        "CodeGenerationAgent": "pycsamt.agents.code_gen",
    }
    out: dict[str, type] = {}
    for name, mod in _MAP.items():
        cls = _import(mod, name)
        if cls is not None:
            out[name] = cls
    return out


def _collect_core() -> dict[str, type]:
    _MAP = {
        "Sites": "pycsamt.site",
        "EDICollection": "pycsamt.edi",
        "EDIFile": "pycsamt.edi",
    }
    out: dict[str, type] = {}
    for name, mod in _MAP.items():
        cls = _import(mod, name)
        if cls is not None:
            out[name] = cls
    return out


# ── Sites key method list ──────────────────────────────────────────────────────

_SITES_KEY_METHODS = [
    "from_any",
    "write",
    "select",
    "edit_all",
    "closest",
    "to_profile",
    "to_edicollection",
    "to_edis",
    "as_list",
    "get",
    "by_index",
    "map",
]

_SITES_ATTRIBUTE_NOTES = """\
Key data attributes (set after construction)
--------------------------------------------
sites[i]         Site object for station i
sites[i].Z       Impedance tensor: shape (n_freq, 2, 2), complex
sites[i].freq    Frequencies array (Hz)
sites[i].rho     Apparent resistivity (Ohm.m), shape (n_freq, 2, 2)
sites[i].phase   Phase (degrees), shape (n_freq, 2, 2)
sites[i].id      Station name / identifier
len(sites)       Number of stations loaded
"""

_ENSURE_SITES_NOTE = """\
ensure_sites(source) — universal EDI loader
-------------------------------------------
Accepts: a Sites object, a directory path (str/Path),
         a list of EDI file paths, or an EDICollection.
Returns: Sites object.
Import:  from pycsamt.ai.inversion._sites_bridge import ensure_sites
Use case: any agent or function that needs to accept
          flexible EDI input without caring about the format.
"""

_WRITE_SITES_NOTE = """\
write_sites(sites, dest, exist_ok=True) — EDI exporter
-------------------------------------------------------
Writes one corrected EDI file per site to dest/.
Import:  from pycsamt.site.export import write_sites
Use case: export corrected Sites after static shift or QC.
"""


# ── agent entry builder (deep) ─────────────────────────────────────────────────


def _agent_entry(name: str, cls: type) -> str:
    """Build a full reference entry for one agent."""
    doc = _full_doc(cls)
    lines: list[str] = []
    lines.append(f"{name}{_sig(cls)}")
    lines.append(f"  {_first_para(doc, 200)}")

    # Parameters section
    params = _section_text(doc, "Parameters")
    if params:
        lines.append("  Parameters:")
        for ln in params.split("\n")[:12]:
            lines.append("    " + ln)

    # Input keys section
    inp = _section_text(doc, "Input keys")
    if inp:
        lines.append("  Input keys:")
        for ln in inp.split("\n")[:8]:
            lines.append("    " + ln)

    # Output data keys section
    out = _section_text(doc, "Output data keys")
    if out:
        lines.append("  Output data keys:")
        for ln in out.split("\n")[:8]:
            lines.append("    " + ln)

    return "\n".join(lines)


# ── tier builders ──────────────────────────────────────────────────────────────

_WORKFLOW_TABLE = """\
Workflow keywords (pass as config["workflow"])
----------------------------------------------
qc                  Data quality control + station flagging
static_shift        Galvanic static-shift detection + correction
phase_analysis      Phase tensor, dimensionality, strike analysis
forward             Forward model computation
pre_inversion       Pre-processing before inversion
inversion_eval      Evaluate inversion results
interpretation      Geological interpretation of resistivity models
report              Generate PDF/HTML processing report
full                Full pipeline: load -> qc -> correct -> invert -> report
ai_inversion        1-D neural network / CNN inversion (EMInverter1D)
inv1d               Alias for ai_inversion
inv2d               2-D U-Net profile inversion
inv3d               3-D GCN graph-convolutional inversion
ensemble_inversion  Ensemble / uncertainty-quantified inversion
pinn_inversion      Physics-Informed Neural Network inversion
hybrid_inversion    PINN + AI inverter hybrid (requires checkpoint)
orchestrated_code_gen  Generate reproducible Python script for a workflow
tipper              Tipper / induction arrow analysis
modem               ModEM 3-D inversion interface
occam2d             Occam2D regularised 2-D inversion interface
"""

_USAGE_EXAMPLES = """\
# 1. Load EDI files and run quality control
from pycsamt.agents import MTLoaderAgent, DataQCAgent
loader = MTLoaderAgent()
sites  = loader.execute({"path": "/data/L22PLT"})["sites"]
qc     = DataQCAgent()
report = qc.execute({"sites": sites})

# 2. Correct static shift
from pycsamt.agents import StaticShiftAgent
agent     = StaticShiftAgent(method="ama")
result    = agent.execute({"sites": sites})
corrected = result["corrected_sites"]
factors   = result["shift_factors"]  # {station: factor}

# 3. Full orchestrated workflow (with LLM router)
from pycsamt.agents import WorkflowOrchestratorAgent
orch   = WorkflowOrchestratorAgent(api_key="sk-...")
result = orch.execute({
    "config":     {"workflow": "qc"},
    "data_path":  "/data/L22PLT",
    "output_dir": "/out/qc/",
})

# 4. Load EDI collection -> Sites
from pycsamt.edi import EDICollection
col   = EDICollection("/data/L22PLT")
sites = col.to_sites()     # -> Sites object

# 5. Universal loader (accepts dir, list, or Sites)
from pycsamt.ai.inversion._sites_bridge import ensure_sites
sites = ensure_sites("/data/L22PLT")     # dir
sites = ensure_sites(["a.edi", "b.edi"]) # list

# 6. Access impedance tensor data
site0 = sites[0]           # first station
Z     = site0.Z            # (n_freq, 2, 2) complex
freq  = site0.freq         # frequencies (Hz)
rho   = site0.rho          # apparent resistivity
phase = site0.phase        # phase (degrees)

# 7. Export corrected EDI files
from pycsamt.site.export import write_sites
write_sites(corrected, "/out/corrected_edis/")

# 8. PINN inversion
from pycsamt.agents import PINNInversionAgent
agent  = PINNInversionAgent(epochs=200, n_layers=8)
result = agent.execute({"sites": sites})
model  = result["model"]  # layered resistivity model

# 9. Generate workflow script
from pycsamt.agents import CodeGenerationAgent
cg     = CodeGenerationAgent()
result = cg.execute({
    "sites": sites,
    "workflow": "qc",
    "output_dir": "/out/",
})
print(result["code"])      # Python source
"""


def _build_tier_core() -> str:
    return (
        "PYCSAMT v2 — PACKAGE OVERVIEW\n"
        "==============================\n"
        "pycsamt is a Python library for"
        " magnetotelluric (MT/AMT/CSAMT) data"
        " processing, correction, and AI-assisted"
        " inversion. Agents are self-contained"
        " workflow components; Sites is the main"
        " data container for loaded EDI data.\n\n" + _WORKFLOW_TABLE
    )


def _build_tier_agents() -> str:
    agents = _collect_agents()
    if not agents:
        return ""
    entries = [_agent_entry(n, c) for n, c in agents.items()]
    body = "\n\n".join(entries)
    return _section("Agent class reference", body)


def _build_tier_sites() -> str:
    core = _collect_core()
    parts: list[str] = []

    sites_cls = core.get("Sites")
    if sites_cls:
        doc = _full_doc(sites_cls)
        params = _section_text(doc, "Parameters")
        attrs = _section_text(doc, "Attributes")
        lines = [
            f"Sites{_sig(sites_cls)}",
            f"  {_first_para(doc, 280)}",
        ]
        if params:
            lines.append("  Constructor parameters:")
            for ln in params.split("\n")[:10]:
                lines.append("    " + ln)
        if attrs:
            lines.append("  Attributes:")
            for ln in attrs.split("\n")[:10]:
                lines.append("    " + ln)

        # key methods
        lines.append("  Key methods:")
        for mname in _SITES_KEY_METHODS:
            m = getattr(sites_cls, mname, None)
            if m:
                lines.append(_method_entry(mname, m, 160))
        parts.append("\n".join(lines))

    # EDICollection
    edi_cls = core.get("EDICollection")
    if edi_cls:
        doc = _full_doc(edi_cls)
        entry = f"EDICollection{_sig(edi_cls)}\n  {_first_para(doc, 240)}"
        to_sites = getattr(edi_cls, "to_sites", None)
        if to_sites:
            entry += "\n" + _method_entry("to_sites", to_sites, 160)
        parts.append(entry)

    # Helper functions
    parts.append(_SITES_ATTRIBUTE_NOTES)
    parts.append(_ENSURE_SITES_NOTE)
    parts.append(_WRITE_SITES_NOTE)

    body = "\n\n".join(parts)
    return _section("Sites / EDI data model reference", body)


def _build_tier_examples() -> str:
    return _section(
        "Usage examples (copy-paste ready)",
        _USAGE_EXAMPLES,
    )


# ── emtools.ss reference notes ────────────────────────────────

_READ_EDIS_NOTE = """\
read_edis -- load many EDI files
---------------------------------
from pycsamt.api import read_edis

Signature:
  read_edis(sources, *, recursive=True,
            strict=False, on_dup="replace",
            progress="auto", leave=False,
            verbose=0) -> APISurvey

sources  : str, Path, list, or glob pattern
Returns APISurvey
  .collection -> Sites / EDICollection
  .name       -> survey name string

Typical use:
    survey = read_edis("L22PLT/")
    sites  = survey.collection
"""

_ESTIMATE_SS_NOTE = """\
estimate_ss_ama -- AMA static-shift factors
--------------------------------------------
from pycsamt.emtools.ss import estimate_ss_ama

Signature:
  estimate_ss_ama(sites, *, sort_by="lon",
    half_window=3, weights="tri",
    pband=None, max_skew=6.0,
    robust_freq="median",
    robust_overall="median",
    recursive=True, on_dup="replace",
    strict=False, verbose=0,
    api=None) -> DataFrame

sites       : Sites, str, Path, list, EDICollection
sort_by     : 'lon'|'lat'|'name' (station order)
half_window : k neighbours each side (default 3)
weights     : 'tri'|'gauss'|'uniform'
pband       : (p_min_s, p_max_s) period band
max_skew    : |beta| threshold (default 6.0)
api         : wrap result in APIFrame when True

Returns DataFrame (one row per station):
  station         -- station name
  delta_log10_rho -- log10 shift (pos=above trend)
  fac_rho         -- 10^(-delta) rho factor
  fac_z           -- 10^(-0.5*delta) Z factor
  n_used          -- frequencies used

Typical use:
    tbl = estimate_ss_ama(
        sites, half_window=3, sort_by="lon"
    )
    print(
        tbl[["station","delta_log10_rho","fac_z"]]
    )
"""

_CORRECT_SS_NOTE = """\
correct_ss_ama -- estimate + apply AMA correction
--------------------------------------------------
from pycsamt.emtools.ss import correct_ss_ama

Signature:
  correct_ss_ama(sites, *, sort_by="lon",
    half_window=3, weights="tri",
    pband=None, max_skew=6.0,
    robust_freq="median",
    robust_overall="median",
    inplace=False, recursive=True,
    on_dup="replace", strict=False,
    verbose=0) -> Sites

Calls estimate_ss_ama then
apply_ss_factors(key="fac_z").
inplace=False returns a corrected copy.

Typical use:
    sites_corr = correct_ss_ama(
        sites, half_window=3, sort_by="lon"
    )
"""

_PLOT_SS_NOTE = """\
plot_ss_summary -- four-panel correction figure
------------------------------------------------
from pycsamt.emtools.ss import plot_ss_summary

Signature:
  plot_ss_summary(
    logRho_before, logRho_after, *,
    freqs, station_labels=None, ...
  ) -> matplotlib.Figure

Panels: (a) Before, (b) After pseudosections,
        (c) Delta section, (d) Per-station bars.

logRho_before, logRho_after : ndarray (n_st, n_f)
  log10 apparent-resistivity arrays.
freqs : ndarray (n_f,)  Hz.
station_labels : list of str or None.

plot_ss_1d_curves -- per-station sounding grid
-----------------------------------------------
from pycsamt.emtools.ss import plot_ss_1d_curves

Signature:
  plot_ss_1d_curves(
    logRho_before, logRho_after, *,
    freqs, station_labels=None,
    n_cols=4, max_stations=16, ...
  ) -> matplotlib.Figure

Grid of subplots (one per station) showing
before/after log10-rho sounding curves.

Build logRho arrays from Sites:
    from pycsamt.emtools._core import (
        _get_z_block, _name, _iter_items,
    )

    def collect_logRho(S):
        rows, freqs = [], None
        for i, ed in enumerate(_iter_items(S)):
            Z, z, fr = _get_z_block(ed)
            if Z is None:
                continue
            rxy = (
                0.2*np.abs(z[:,0,1])**2
                /(fr+1e-24)
            )
            ryx = (
                0.2*np.abs(z[:,1,0])**2
                /(fr+1e-24)
            )
            rows.append(
                np.log10(
                    np.sqrt(rxy*ryx)+1e-24
                )
            )
            freqs = fr
        return np.array(rows), freqs

    logRho_b, freqs = collect_logRho(sites)
    logRho_a, _     = collect_logRho(sites_corr)
    labels = [
        _name(ed, i)
        for i, ed in enumerate(_iter_items(sites))
    ]
    fig = plot_ss_summary(
        logRho_b, logRho_a,
        freqs=freqs, station_labels=labels,
    )
    fig.savefig("ss_summary.png", dpi=150)
"""

_SS_WORKFLOW_EXAMPLE = """\
# Static-shift correction full workflow
from pycsamt.api import read_edis
from pycsamt.emtools.ss import (
    estimate_ss_ama,
    correct_ss_ama,
    plot_ss_summary,
    plot_ss_1d_curves,
)
from pycsamt.agents import StaticShiftAgent
from pycsamt.emtools._core import (
    _get_z_block, _name, _iter_items,
)
import numpy as np

# 1. Load EDI files
survey = read_edis("L22PLT/")
sites  = survey.collection

# 2. Inspect shift factors (optional)
ss_table = estimate_ss_ama(
    sites,
    half_window=3,
    sort_by="lon",
    weights="tri",
    max_skew=6.0,
)
print(
    ss_table[[
        "station","delta_log10_rho","fac_z"
    ]]
)

# 3. Correct Z tensor in place
sites_corr = correct_ss_ama(
    sites,
    half_window=3,
    sort_by="lon",
)

# 4. (Alternative) via StaticShiftAgent
agent  = StaticShiftAgent(method="ama")
result = agent.execute({"sites": sites})
sites_corr2 = result["corrected_sites"]
factors     = result["shift_factors"]

# 5. Build log10-rho arrays for plots
def collect_logRho(S):
    rows, freqs = [], None
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        rxy = (
            0.2*np.abs(z[:,0,1])**2/(fr+1e-24)
        )
        ryx = (
            0.2*np.abs(z[:,1,0])**2/(fr+1e-24)
        )
        rows.append(
            np.log10(np.sqrt(rxy*ryx)+1e-24)
        )
        freqs = fr
    return np.array(rows), freqs

logRho_b, freqs = collect_logRho(sites)
logRho_a, _     = collect_logRho(sites_corr)
labels = [
    _name(ed, i)
    for i, ed in enumerate(_iter_items(sites))
]

# 6. Summary figure
fig_sum = plot_ss_summary(
    logRho_b, logRho_a,
    freqs=freqs,
    station_labels=labels,
)
fig_sum.savefig(
    "ss_summary.png", dpi=150,
    bbox_inches="tight",
)

# 7. Per-station 1-D sounding curves
fig_1d = plot_ss_1d_curves(
    logRho_b, logRho_a,
    freqs=freqs,
    station_labels=labels,
)
fig_1d.savefig(
    "ss_1d_curves.png", dpi=150,
    bbox_inches="tight",
)
"""


def _build_tier_emtools() -> str:
    parts = [
        _READ_EDIS_NOTE,
        _ESTIMATE_SS_NOTE,
        _CORRECT_SS_NOTE,
        _PLOT_SS_NOTE,
        _SS_WORKFLOW_EXAMPLE,
    ]
    body = "\n\n".join(parts)
    return _section(
        "emtools.ss / Static-shift reference",
        body,
    )


# ── public API ─────────────────────────────────────────────────────────────────


def build_tiers() -> dict[str, str]:
    r"""
    Build all context tiers.

    Returns
    -------
    dict with keys: 'core', 'agents', 'sites',
                    'examples', 'emtools'
    """
    return {
        "core": _build_tier_core(),
        "agents": _build_tier_agents(),
        "sites": _build_tier_sites(),
        "examples": _build_tier_examples(),
        "emtools": _build_tier_emtools(),
    }


def build_package_context(
    max_agent_doc: int = 240,
    max_core_doc: int = 400,
) -> str:
    r"""
    Build the full pycsamt API reference string.

    Joins all tiers.  Suitable for LLM system
    prompts and offline docstring lookup.

    Parameters
    ----------
    max_agent_doc : int
        Unused — kept for backwards compatibility.
    max_core_doc : int
        Unused — kept for backwards compatibility.
    """
    tiers = build_tiers()
    return "\n".join(tiers.values())


# ── backward-compat re-exports ────────────────────────────────────────────────


def _short_doc(obj: Any, max_chars: int = 300) -> str:
    return _first_para(_full_doc(obj), max_chars)


def _sig_compat(cls: type, max_chars: int = 120) -> str:
    return _sig(cls, max_chars)


# ── module-level singletons ────────────────────────────────────────────────────

# Built once at import — safe in LLM system prompts.
_TIERS: dict[str, str] = build_tiers()

TIER_CORE = _TIERS["core"]
TIER_AGENTS = _TIERS["agents"]
TIER_SITES = _TIERS["sites"]
TIER_EXAMPLES = _TIERS["examples"]
TIER_EMTOOLS = _TIERS["emtools"]

PACKAGE_CONTEXT: str = "\n".join(_TIERS.values())

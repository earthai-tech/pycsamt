# AGENTS.md — working with pyCSAMT v2

Guidance for code assistants (and humans) working in this repository.
This file is also indexed by the RAG assistant, so keep it accurate.

## What this project is

**pyCSAMT v2** is a Python suite for **CSAMT / AMT / MT** electromagnetic
data: from raw EDI files to quality control, distortion correction,
phase-tensor analysis, AI-assisted inversion, and figures — with an
interactive web app, a desktop app, and a conversational *Agent Master*.

- Package: `pycsamt` (import name, **lowercase**). Brand/display name:
  **pyCSAMT**.
- Version: 2.0.0. License: LGPL-3.0. Author: L. Kouadio.

## Repository map (what to touch, what to avoid)

| Area | Path | Notes |
|------|------|-------|
| EM tools (science API) | `pycsamt/emtools/` | `ss.py` (static shift), `qc.py`, `_core.py` |
| Data model | `pycsamt/site/`, `pycsamt/seg/`, `pycsamt/z/` | `Sites`, EDI collection, impedance |
| Forward / inversion | `pycsamt/forward/`, `pycsamt/ai/` | MT1D forward, EMInverter1D |
| Agent layer | `pycsamt/agents/` | routing + per-workflow agents |
| Assistant / RAG | `pycsamt/assistant/` | retrieval, tools, memory, evals |
| Apps (UI) | `pycsamt/app/` | desktop / web / agent_master (Dash) |
| Docs | `docs/source/` | reStructuredText |
| Example data | `data/AMT/WILLY_DATA/` | survey lines L18/L22/L26/L30/L34 PLT |

Do **not** treat `pycsamt/app/` internals as the science API — the
workflows live in `pycsamt/emtools` and `pycsamt/agents`.

## Core data model

```python
from pycsamt.emtools._core import ensure_sites
sites = ensure_sites(edi_dir, recursive=True, verbose=0)   # -> Sites
```

`ensure_sites` is the canonical loader (EDI / directory / Sites). A
station needs valid impedance (Z); files without Z are skipped.

## Common workflows (use these REAL symbols)

| Task | High-level agent | Low-level functions |
|------|------------------|---------------------|
| Load EDI | `MTLoaderAgent` | `ensure_sites` |
| Quality control | `DataQCAgent` | `pycsamt.emtools.qc.build_qc_table`, `qc_flags` |
| Static shift | `StaticShiftAgent` | `pycsamt.emtools.ss.estimate_ss_ama`, `correct_ss_ama`, `apply_ss_factors`, `plot_ss_summary`, `plot_ss_1d_curves` |
| Phase tensor / strike | `PhaseAnalysisAgent` | — |
| AI 1-D inversion | `AIInversionAgent` | `pycsamt.ai.inversion.inv1d.EMInverter1D` |
| Code generation | `CodeGenerationAgent` | — |

Agents return an `AgentResult` (`status`, `summary`, `data`, …) and are
called as `Agent(...).execute({"path": edi_dir, "output_dir": out})`.

## Agent / assistant architecture

```
user request
  → IntentRouter        (question | code | plot | workflow | meta | clarify)
  → ContextInputAgent   (workflow classification; shared keyword registry
                         pycsamt.agents._workflows)
  → WorkflowOrchestratorAgent  (load → qc → … → report chain)
  → AgentResult
```

The **RAG assistant** (`pycsamt/assistant`) grounds answers/code in real
package facts:

```bash
python -m pycsamt.assistant.rag build     # build the index (.pycsamt_rag/)
python -m pycsamt.assistant.rag query "static shift"
python -m pycsamt.assistant.rag eval      # score the assistant
```

- `assistant/rag/` — ingest, retriever, context builder
- `assistant/tools/` — `resolve_line`, `validate_generated_code`,
  `run_static_shift` / `run_qc` / `run_phase_analysis`
- `assistant/memory/` — session / project state, workflow trace
- `projects/willy_project_registry.yml` — survey-line → path + defaults
- `assistant_recipes/*.md` — per-workflow recipes

## Conventions

- Use `pycsamt` (lowercase) in code; `pyCSAMT` only as a display name.
- Keep the workflow keyword table in `pycsamt/agents/_workflows.py` as the
  single source of truth (both `ContextInputAgent` and the orchestrator
  import it).
- Tests live in `*/tests/test_*.py`; run with the `pycsamt-v2` conda env.
- Match the surrounding code style (LGPL header, type hints, docstrings).

## Do / Don't (avoid hallucinations)

- DO call real functions: `estimate_ss_ama`, `correct_ss_ama`,
  `ensure_sites`, `StaticShiftAgent`.
- DON'T invent APIs such as `from pycsamt import static_shift` or
  `pycsamt.run_workflow(...)` — verify a symbol exists first
  (`pycsamt.assistant.tools.validate_generated_code` does this).
- Static-shift factors must be finite and positive; never apply a NaN/0
  factor to impedance.
- AMA static shift assumes ~1-D/2-D structure; on strongly 3-D data it
  may decline (empty factor table) — that is expected, not a bug.
- Resolve survey lines via the project registry, not hard-coded paths.

# Full survey workflow recipe (load → QC → preprocess → invert → report)

## User intent phrases
- load the EDI dataset, run quality control and preprocessing, prepare the
  inversion, run AI-assisted resistivity inversion, evaluate the result,
  and generate a reproducible report
- full pipeline, full workflow, end to end, complete survey processing
- process everything and give me a report

## What the full workflow runs (in order)
1. **Load** — read the EDI files into a site collection (`MTLoaderAgent`).
2. **Quality control** — per-station QC scores and flags (`DataQCAgent`).
3. **Static-shift correction** — preprocessing (`StaticShiftAgent`).
4. **Phase-tensor & strike analysis** (`PhaseAnalysisAgent`).
5. **Denoising** — RPCA noise removal (`DenoisingAgent`).
6. **AI 1-D inversion** — AI-assisted resistivity model (`AIInversionAgent`).
7. **Occam2D preparation** — classical inversion input files (`Occam2DAgent`).
8. **Report** — a reproducible survey report with figures (`ReportAgent`).

## Required inputs
- EDI directory or survey line (e.g. `data/AMT/WILLY_DATA/L22PLT`)
- output directory (every step writes into it; the report cites the outputs)

## Preferred API
```python
from pycsamt.agents.orchestrator import WorkflowOrchestratorAgent

result = WorkflowOrchestratorAgent().execute({
    "request": "run the full pipeline and generate a report",
    "data_path": "data/AMT/WILLY_DATA/L22PLT",
    "output_dir": "outputs/full_run",
})
print(result.summary)
```

In the Agent Master chat, simply describe the compound task — a request
that names several stages (QC + inversion + report) classifies as the
`full` workflow automatically:

> Load this EDI field dataset, perform quality control and preprocessing,
> prepare the inversion, run AI-assisted resistivity inversion, evaluate
> the result, and generate a reproducible report.

## Key symbols
- `pycsamt.agents.orchestrator.WorkflowOrchestratorAgent`
- `pycsamt.agents._workflows.classify_workflow`
- `pycsamt.agents.qc.DataQCAgent`
- `pycsamt.agents.static_shift.StaticShiftAgent`
- `pycsamt.agents.ai_inversion.AIInversionAgent`
- `pycsamt.agents.report.ReportAgent`

## Expected outputs
- per-step results and figures under the output directory
- an AI resistivity model with evaluation metrics
- Occam2D input files ready for a classical run
- a survey report (with the QC, processing, and inversion outcomes)
  that makes the run reproducible

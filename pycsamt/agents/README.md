# pycsamt.agents — AI-assisted MT/AMT workflow automation

A fully modular, LLM-augmented agent system for magnetotelluric data
processing.  Each agent wraps one well-defined processing step; the
`AgentCoordinator` chains them into reproducible workflows.

---

## Table of contents

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [Agent catalogue](#agent-catalogue)
4. [Workflow patterns](#workflow-patterns)
5. [LLM configuration](#llm-configuration)
6. [CLI reference](#cli-reference)
7. [Gradio web interface](#gradio-web-interface)
8. [Model zoo](#model-zoo)
9. [Building custom agents](#building-custom-agents)
10. [Cost tracking](#cost-tracking)

---

## Installation

```bash
pip install pycsamt                        # core (numpy, scipy, matplotlib)
pip install anthropic                      # Claude LLM (default)
pip install openai                         # optional — OpenAI provider
pip install google-generativeai            # optional — Gemini provider
pip install torch                          # optional — AI inversion / denoising
pip install gradio                         # optional — web UI
```

---

## Quick start

### No LLM, no GPU — pure regex parsing

```python
from pycsamt.agents import ContextInputAgent, MTLoaderAgent, AgentCoordinator

ctx    = ContextInputAgent()   # no api_key → regex fallback
loader = MTLoaderAgent()

coord = AgentCoordinator("my_qc_workflow")
coord.add_step("parse", ctx,   description="Parse NL request")
coord.add_step("load",  loader,
               input_fn=lambda r: {"path": r["parse"]["config"]["data_path"]},
               description="Load EDI files")

result = coord.execute(
    {"request": "Load /data/WILLY_EDIs, QC, period 1e-4 to 1 s"},
    dry_run=True,   # preview only — no files touched
)
print(result.summary)
```

### With a Claude API key

```python
from pycsamt.agents import WorkflowOrchestratorAgent

agent = WorkflowOrchestratorAgent(
    api_key="sk-ant-…",
    llm_provider="claude",
)
r = agent.execute({
    "request":    "Load WILLY data, QC, run full phase-tensor analysis",
    "data_path":  "/data/WILLY_EDIs",
    "output_dir": "/out/willy_survey",
})
print(r["workflow_type"])          # "phase_analysis"
print(r.llm_interpretation)       # LLM geological narrative
```

### Pre-trained AI inversion (Phase 5 model zoo)

```python
from pycsamt.agents import AIInversionAgent

agent = AIInversionAgent.from_pretrained("mt1d-resnet-5layer-v1")
r = agent.execute({"path": "/data/WILLY_EDIs"})
print(r["rms_global"])
```

---

## Agent catalogue

### Phase 1 — Foundation

| Agent | Description |
|---|---|
| `ContextInputAgent` | Natural-language → structured workflow JSON (LLM or regex fallback) |
| `MTLoaderAgent` | Load EDI / AVG / J files → `Sites`; per-station completeness + QC table |
| `AgentCoordinator` | Chain agents, checkpoint state, track total cost |

### Phase 2 — Core processing pipeline

| Agent | Description |
|---|---|
| `DataQCAgent` | SNR section, dead-band detection, station QC scores (0–100) |
| `StaticShiftAgent` | AMA / LOESS / spatial-median static-shift detection & correction |
| `PhaseAnalysisAgent` | Full phase-tensor suite: PT section, strike, Bahr skew, Mohr, Argand |
| `ForwardModelAgent` | 1-D / 2-D / 3-D MT forward modelling from a `LayeredModel` |
| `InversionPrepAgent` | Write Occam2D data + mesh + startup files |
| `InversionEvaluationAgent` | Load inversion result, compute RMS & residual PT section |
| `InterpretationAgent` | Map resistivity ranges → lithology; correlate with borehole logs |
| `ReportAgent` | Assemble all figures + tables into Markdown / HTML / PDF |
| `CodeGenerationAgent` | Emit a standalone reproducible Python script from workflow config |

### Phase 3 — Advanced processing & specialised outputs

| Agent | Description |
|---|---|
| `WorkflowOrchestratorAgent` | NL → classify workflow → build & run the correct agent chain |
| `DenoisingAgent` | RPCA / Hampel / EMAP / AI-CAE denoising |
| `AIInversionAgent` | End-to-end 1-D AI inversion (ResNet / CNN / FCN) |
| `Occam2DAgent` | Write Occam2D data + mesh + startup files (higher-level than `InversionPrepAgent`) |
| `ModEmAgent` | Write ModEM3D impedance data file |
| `AnomalyDetectionAgent` | Unsupervised CAE anomaly flagging per (station, frequency) |

### Phase 4 — Deep-learning architectures

| Agent | Description |
|---|---|
| `Inv2DAgent` | U-Net 2-D profile inversion (captures lateral continuity) |
| `EnsembleAgent` | Ensemble 1-D inversion with conformal uncertainty bands |
| `JointInversionAgent` | DRCNN multi-modal joint inversion (MT + TEM / CSAMT / gravity) |

### Phase 5 — Model zoo & pre-trained weights

| Agent | Description |
|---|---|
| `ModelZooAgent` | Browse, download, and deploy pre-trained EM inverter checkpoints |

### Phase 6 — IoT field acquisition

| Agent | Description |
|---|---|
| `IoTFieldAgent` | Monitor a live IoT AMT/MT/CSAMT/CSEM acquisition from edge telemetry: monitoring status, telemetry/station/pipeline tables, sync & power summaries, provenance manifest, and dashboard/edge/power/sync figures |

```python
from pycsamt.agents import IoTFieldAgent

# from live telemetry packets (offline, no LLM key needed)
r = IoTFieldAgent().execute({
    "packets": packets,          # or "session"=FieldSession, or "path"=EDIs
    "survey_id": "WILLY_IOT",
    "manifest": True,            # build + write a signed-capable manifest
    "output_dir": "/out/iot",
})
r["level"]           # 'ok' | 'warn' | 'critical'
r["status"].issues   # concrete reasons the level tripped
r["figures"]["dashboard"]
```

---

## Workflow patterns

### Pattern 1 — Data QC & preprocessing (~5 min)

```
NL → ContextInputAgent → MTLoaderAgent → DataQCAgent
   → StaticShiftAgent → ReportAgent
```

```python
from pycsamt.agents import WorkflowOrchestratorAgent
r = agent.execute({"request": "QC the WILLY EDIs", "data_path": "/data/WILLY"})
```

### Pattern 2 — Phase-tensor analysis (~8 min)

```
MTLoaderAgent → DataQCAgent → StaticShiftAgent
→ PhaseAnalysisAgent → ReportAgent
```

### Pattern 3 — Pre-inversion pipeline (~10 min)

```
MTLoaderAgent → DataQCAgent → StaticShiftAgent → PhaseAnalysisAgent
→ Occam2DAgent → CodeGenerationAgent
```

### Pattern 4 — AI 1-D inversion (~15 min)

```
MTLoaderAgent → DataQCAgent → DenoisingAgent
→ AIInversionAgent → InterpretationAgent → ReportAgent
```

### Pattern 4b — 2-D AI inversion (U-Net)

```python
from pycsamt.agents import Inv2DAgent
r = Inv2DAgent().execute({"path": "/data/WILLY", "output_dir": "/out"})
```

### Pattern 4c — Ensemble with uncertainty

```python
from pycsamt.agents import EnsembleAgent
r = EnsembleAgent(n_estimators=5, epochs=50).execute({"path": "/data/WILLY"})
print(r["coverage"])   # empirical 90 % coverage
```

### Pattern 4d — Multi-modal joint inversion

```python
from pycsamt.agents import JointInversionAgent
r = JointInversionAgent(modalities=["mt", "tem"]).execute({
    "path":           "/data/WILLY",
    "secondary_path": "/data/WILLY_TEM",   # optional
})
```

### Pattern 5 — Full survey (~30 min)

```
MTLoaderAgent → DataQCAgent → StaticShiftAgent → PhaseAnalysisAgent
→ DenoisingAgent → AIInversionAgent → Occam2DAgent → ReportAgent
```

```python
r = agent.execute({"request": "full pipeline for WILLY data", "data_path": "…"})
```

---

## LLM configuration

```python
# Claude (default, recommended)
agent = DataQCAgent(api_key="sk-ant-…", llm_provider="claude",
                    model="claude-sonnet-4-6")

# OpenAI
agent = DataQCAgent(api_key="sk-…", llm_provider="openai",
                    model="gpt-4o")

# Gemini
agent = DataQCAgent(api_key="AIza…", llm_provider="gemini",
                    model="gemini-1.5-pro")

# No LLM (free, regex/rule-based fallback)
agent = DataQCAgent()   # api_key omitted → no LLM calls
```

Every `AgentResult` includes a `cost_estimate_usd` field; call
`result.cost_estimate_usd` to track spend per agent or per workflow.

---

## CLI reference

```bash
# Preview a workflow (no files written, no LLM calls)
python -m pycsamt.agents preview "Load /data/WILLY EDIs, QC, PT analysis"

# List all agents across all phases
python -m pycsamt.agents list

# Show LLM pricing tables
python -m pycsamt.agents pricing

# Browse the model zoo (no network required)
python -m pycsamt.agents zoo

# Download a specific checkpoint
python -m pycsamt.agents zoo mt1d-resnet-5layer-v1

# Launch the Gradio web UI
python -m pycsamt.agents web
python -m pycsamt.agents web --share --port=8080
```

---

## Gradio web interface

```python
from pycsamt.agents.web import launch
launch()    # opens http://localhost:7860
```

Tabs:

| Tab | Agent(s) |
|---|---|
| 💬 Chat | `WorkflowOrchestratorAgent` — plain-English dispatch |
| 📂 Load & QC | `MTLoaderAgent` + `DataQCAgent` |
| ⚡ Static Shift | `StaticShiftAgent` |
| 🔬 Phase Analysis | `PhaseAnalysisAgent` |
| 📡 Forward Model | `ForwardModelAgent` |
| 🤖 AI Inversion (1-D) | `AIInversionAgent` |
| 🗺 2-D AI Inversion | `Inv2DAgent` (U-Net) |
| 📊 Ensemble Inversion | `EnsembleAgent` |
| 🔗 Joint Inversion | `JointInversionAgent` |
| 🏛 Model Zoo | `ModelZooAgent` |
| 📄 Report | `ReportAgent` |

---

## Model zoo

```python
from pycsamt.agents import ModelZooAgent

zoo = ModelZooAgent()

# List all available models
r = zoo.execute({"action": "list"})
for name, desc in r["models"].items():
    print(name, "—", desc)

# Download a checkpoint
r = zoo.execute({"action": "download", "model_name": "mt1d-resnet-5layer-v1"})
print(r["checkpoint_path"])

# Zero-shot prediction on observed EDIs
r = zoo.execute({
    "action":     "predict",
    "model_name": "mt1d-resnet-5layer-v1",
    "path":       "/data/WILLY_EDIs",
    "output_dir": "/out/zoo",
})
print(r["rms_global"])
```

Available checkpoints (weights published when Phase 5 is released):

| Model | Arch | Layers | Solver | Training data |
|---|---|---|---|---|
| `mt1d-resnet-5layer-v1` | ResNet | 5 | MT 1-D | 50 k synthetic profiles |
| `mt1d-cnn-5layer-v1` | CNN | 5 | MT 1-D | 50 k synthetic profiles |
| `mt1d-resnet-7layer-v1` | ResNet | 7 | MT 1-D | 50 k synthetic profiles |
| `csamt1d-resnet-5layer-v1` | ResNet | 5 | CSAMT 1-D | 30 k synthetic profiles |
| `tem1d-fcn-5layer-v1` | FCN | 5 | TEM 1-D | 20 k synthetic profiles |

---

## Building custom agents

```python
from pycsamt.agents import BaseAgent, AgentResult

class MyCustomAgent(BaseAgent):
    SYSTEM_PROMPT = "You are a custom MT processing expert…"

    def __init__(self, *, api_key=None, model=None, llm_provider="claude"):
        super().__init__("MyCustomAgent", api_key=api_key,
                         model=model, llm_provider=llm_provider,
                         section_preset="pseudosection")

    def execute(self, input_data):
        sites = input_data.get("sites")
        # … your processing logic …
        interp = self.query_llm("Interpret this result…", max_tokens=200)
        return AgentResult(
            status="success",
            summary="Custom agent complete.",
            data={"result": 42},
            warnings=[],
            llm_interpretation=interp,
            elapsed_seconds=0.1,
            cost_estimate_usd=self._last_cost,
        )
```

Register in a coordinator:

```python
from pycsamt.agents import MTLoaderAgent, AgentCoordinator

coord = AgentCoordinator("custom_workflow")
coord.add_step("load", MTLoaderAgent(), description="Load EDIs")
coord.add_step("custom", MyCustomAgent(),
               input_fn=lambda r: {"sites": r["load"]["sites"]},
               description="My custom step")
result = coord.execute({"path": "/data/EDIs"})
```

---

## Cost tracking

```python
from pycsamt.agents import estimate_cost, format_cost, PROVIDER_RATES

# estimate cost for a single LLM call
cost = estimate_cost(
    provider="claude",
    model="claude-sonnet-4-6",
    input_tokens=500,
    output_tokens=150,
)
print(format_cost(cost))   # "$0.000525"

# total cost from a coordinator run
result = coord.execute(config)
print(f"Total workflow cost: {format_cost(result.cost_estimate_usd)}")
```

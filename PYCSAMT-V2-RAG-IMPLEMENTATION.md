## 1. What I found in the updated `pycsamt-v2` zip

I inspected the new archive:

```text
pycsamt-archive-updated-new.zip
```

It contains a much cleaner and richer version of the package than before.

### Package size

| Item                       |           Count / observation |
| -------------------------- | ----------------------------: |
| Total files                |                        ~2,440 |
| Python files               |                          ~970 |
| `.pyc` bytecode files      |                          ~648 |
| Documentation `.rst` files |                          ~118 |
| EDI files                  |                          ~131 |
| Fortran files              |                           ~99 |
| Package version            |               `pycsamt 2.0.0` |
| Main package               |                    `pycsamt/` |
| Main docs                  |                `docs/source/` |
| Main data example          | `data/AMT/WILLY_DATA/L18PLT/` |
| LLM/agent layer            |             `pycsamt/agents/` |
| Pipeline layer             |           `pycsamt/pipeline/` |
| Static-shift functions     |       `pycsamt/emtools/ss.py` |

Important: the archive contains many useful files, but also many files that should **not** be indexed in RAG, such as:

```text
__pycache__/
*.pyc
docs/source/_static/
image assets
desktop/web UI assets
generated files
```

Also, `pyproject.toml` points to:

```toml
readme = "README.md"
```

but I did **not** find a root `README.md` in the zip. You should add one, because RAG systems and documentation agents rely strongly on a good root README.

---

# 2. Good news: your new pycsamt-v2 already has an AI-agent foundation

Your new package already contains several important components:

```text
pycsamt/agents/
    router.py
    package_qa.py
    _package_context.py
    code_gen.py
    static_shift.py
    loader.py
    qc.py
    phase_analysis.py
    ai_inversion.py
    inv2d_agent.py
    inv3d_agent.py
    modem_agent.py
    occam2d_agent.py
    pipeline_agent.py
```

This is excellent. You already started building an assistant system.

For example, I found:

| Component               | Role                                                                  |
| ----------------------- | --------------------------------------------------------------------- |
| `IntentRouter`          | Classifies user intent: question, code, plot, workflow, meta, clarify |
| `PackageQAAgent`        | Answers questions about pycsamt using package context                 |
| `_package_context.py`   | Builds runtime package knowledge from docstrings                      |
| `CodeGenerationAgent`   | Generates reproducible Python scripts                                 |
| `StaticShiftAgent`      | Runs static-shift correction                                          |
| `MTLoaderAgent`         | Loads EDI/AVG/J files into `Sites`                                    |
| `_workflows.py`         | Central registry of workflow keywords                                 |
| `pipeline/_registry.py` | Registry of processing steps, including static shift                  |

So the problem is **not that pycsamt-v2 has no assistant layer**.

The problem is that the current assistant is still mostly:

```text
LLM prompt + package context + workflow agents
```

but what you need is:

```text
LLM + RAG + AST code index + workflow tools + data registry + execution validation
```

That is the next step.

---

# 3. What RAG should do in pycsamt-v2

For pycsamt-v2, RAG should be the assistant’s **technical memory**.

It should retrieve:

```text
Which class to use
Which function to call
Which parameters are valid
Which workflow is appropriate
Which example matches the user request
Which data path corresponds to a survey line
Which plotting function should be used
Which output files are expected
```

For example, when the user asks:

```text
write the static shift code for processing my EDI data of line L22PLT
```

The RAG system should retrieve:

```text
pycsamt/agents/static_shift.py
pycsamt/emtools/ss.py
pycsamt/agents/code_gen.py
pycsamt/agents/loader.py
docs/source/agents/processing_agents.rst
docs/source/tutorials/static_shift*.rst
project_registry.yml: L22PLT → /real/path/to/L22PLT
```

Then the assistant should answer using the real package API, not by guessing.

---

# 4. Recommended RAG directory structure

I recommend adding a new assistant/RAG module inside the package:

```text
pycsamt/
│
├── assistant/
│   ├── __init__.py
│   ├── rag/
│   │   ├── __init__.py
│   │   ├── config.py
│   │   ├── ingest.py
│   │   ├── chunkers.py
│   │   ├── ast_indexer.py
│   │   ├── doc_indexer.py
│   │   ├── retriever.py
│   │   ├── reranker.py
│   │   ├── context_builder.py
│   │   ├── prompt_builder.py
│   │   ├── answerer.py
│   │   ├── citations.py
│   │   ├── evals.py
│   │   └── schemas.py
│   │
│   ├── tools/
│   │   ├── __init__.py
│   │   ├── project_registry.py
│   │   ├── data_tools.py
│   │   ├── static_shift_tools.py
│   │   ├── qc_tools.py
│   │   ├── plot_tools.py
│   │   ├── inversion_tools.py
│   │   ├── code_tools.py
│   │   └── validation_tools.py
│   │
│   ├── memory/
│   │   ├── __init__.py
│   │   ├── session_state.py
│   │   ├── project_state.py
│   │   └── workflow_history.py
│   │
│   └── apps/
│       ├── cli.py
│       └── chat_server.py
│
├── assistant_recipes/
│   ├── load_edi.md
│   ├── static_shift.md
│   ├── qc.md
│   ├── phase_tensor.md
│   ├── occam2d.md
│   ├── modem3d.md
│   ├── ai_1d_inversion.md
│   ├── ai_2d_inversion.md
│   ├── ai_3d_inversion.md
│   ├── plotting.md
│   └── report_generation.md
│
├── projects/
│   ├── example_project_registry.yml
│   └── willy_project_registry.yml
│
└── AGENTS.md
```

This structure separates:

| Layer                | Purpose                          |
| -------------------- | -------------------------------- |
| `assistant/rag/`     | Retrieval and indexing           |
| `assistant/tools/`   | Real executable pycsamt tools    |
| `assistant/memory/`  | Session/project state            |
| `assistant_recipes/` | Human-written workflow recipes   |
| `projects/`          | Survey-line registry             |
| `AGENTS.md`          | Instructions for code assistants |

---

# 5. What to index and what to exclude

## A. Exclude from RAG

Do **not** index these files:

```text
**/__pycache__/**
**/*.pyc
docs/source/_static/**
docs/source/images/**
pycsamt/app/desktop/resources/**
pycsamt/app/web/assets/**
*.png
*.svg
*.ico
*.jpg
*.jpeg
*.gif
.pytest_cache/**
.ruff_cache/**
build/**
dist/**
```

These files pollute retrieval.

## B. Index with high priority

Index these files strongly:

```text
pycsamt/agents/*.py
pycsamt/emtools/*.py
pycsamt/api/*.py
pycsamt/site/*.py
pycsamt/seg/*.py
pycsamt/z/*.py
pycsamt/zonge/*.py
pycsamt/pipeline/*.py
pycsamt/forward/*.py
pycsamt/inversion/*.py
pycsamt/models/**/*.py
pycsamt/ai/**/*.py
docs/source/**/*.rst
examples/**/*.py
pycsamt/**/tests/test_*.py
pyproject.toml
LICENSE.md
```

## C. Index EDI data only as metadata

Do **not** embed entire EDI files as normal text. Instead, create metadata summaries:

```json
{
  "line": "L18PLT",
  "path": "data/AMT/WILLY_DATA/L18PLT",
  "n_edi_files": 25,
  "stations": ["18-001A", "18-002U", "..."],
  "data_type": "AMT/EDI",
  "source": "example_data"
}
```

The assistant needs to know that data exists, but it does not need to embed every EDI block unless the user asks for station-level inspection.

---

# 6. The most important addition: project registry

Your assistant cannot understand this:

```text
line L22PLT
```

unless you provide a registry.

The new zip contains example line:

```text
data/AMT/WILLY_DATA/L18PLT/
```

but I did not see `L22PLT` in the inspected top-level sample. So for `L22PLT`, you need a real project registry.

Create:

```text
projects/willy_project_registry.yml
```

Example:

```yaml
project:
  name: WILLY_DATA
  type: AMT
  default_output_root: results/willy

lines:
  L18PLT:
    edi_dir: data/AMT/WILLY_DATA/L18PLT
    station_pattern: "18-*"
    sort_by: lon
    default_workflows:
      - load
      - qc
      - static_shift
      - phase_analysis

  L22PLT:
    edi_dir: /absolute/path/to/WILLY_DATA/L22PLT
    station_pattern: "22-*"
    sort_by: lon
    default_workflows:
      - load
      - qc
      - static_shift
      - phase_analysis
      - occam2d

static_shift_defaults:
  method: ama
  half_window: 3
  weights: tri
  max_skew: 6.0
  pband: null

plot_defaults:
  dpi: 150
  format: png
```

Then create a tool:

```python
def resolve_line(line_name: str) -> dict:
    """
    Resolve a survey line name such as L22PLT to:
    - EDI directory
    - output directory
    - station pattern
    - default processing parameters
    """
```

This is essential.

Without this, the assistant will always guess.

---

# 7. Recommended RAG indexes

Use multiple indexes, not one big mixed index.

```text
rag_store/
│
├── docs_index/
│   └── tutorials, theory, user guide, CLI docs
│
├── api_index/
│   └── public classes, functions, signatures, docstrings
│
├── code_index/
│   └── implementation chunks from .py files
│
├── tests_index/
│   └── real usage from tests
│
├── recipes_index/
│   └── assistant_recipes/*.md
│
├── workflow_index/
│   └── workflows, router, pipeline step registry
│
├── data_index/
│   └── EDI line/station metadata
│
└── project_index/
    └── project registry, user survey aliases
```

Each index has a different job.

| Index            | Used when                                      |
| ---------------- | ---------------------------------------------- |
| `docs_index`     | User asks conceptual or tutorial questions     |
| `api_index`      | User asks how to use a class/function          |
| `code_index`     | User asks implementation details or bug fixes  |
| `tests_index`    | User asks for working code examples            |
| `recipes_index`  | User asks for practical workflows              |
| `workflow_index` | User asks “run/process/apply/prepare”          |
| `data_index`     | User mentions line/station names               |
| `project_index`  | User mentions `L22PLT`, `L18PLT`, survey names |

---

# 8. AST-based indexing is necessary

For Python code, normal text chunking is not enough.

You need to index by:

```text
module
class
function
method
signature
docstring
parameters
return values
examples
related tests
related docs
```

Example metadata for `StaticShiftAgent`:

```json
{
  "symbol": "pycsamt.agents.static_shift.StaticShiftAgent",
  "type": "class",
  "file": "pycsamt/agents/static_shift.py",
  "workflow": "static_shift",
  "summary": "Detect and correct galvanic static shift in MT/AMT data.",
  "input_keys": ["sites", "path", "method", "output_dir"],
  "output_keys": [
    "corrected_sites",
    "shift_factors",
    "rho_before",
    "rho_after",
    "delta_stats",
    "figures",
    "figure_paths"
  ],
  "methods": ["execute"],
  "related_symbols": [
    "pycsamt.emtools.ss.estimate_ss_ama",
    "pycsamt.emtools.ss.correct_ss_ama",
    "pycsamt.emtools.ss.plot_ss_summary",
    "pycsamt.emtools.ss.plot_ss_1d_curves"
  ]
}
```

Example metadata for `estimate_ss_ama`:

```json
{
  "symbol": "pycsamt.emtools.ss.estimate_ss_ama",
  "type": "function",
  "file": "pycsamt/emtools/ss.py",
  "workflow": "static_shift",
  "parameters": [
    "sites",
    "sort_by",
    "half_window",
    "weights",
    "pband",
    "max_skew",
    "robust_freq",
    "robust_overall",
    "recursive",
    "on_dup",
    "strict",
    "verbose",
    "api"
  ],
  "returns": "pandas.DataFrame or APIFrame",
  "return_columns": [
    "station",
    "delta_log10_rho",
    "fac_rho",
    "fac_z",
    "n_used"
  ]
}
```

This is the kind of knowledge your assistant needs.

---

# 9. Proposed RAG ingestion pipeline

Create a command:

```bash
pycsamt-rag build-index \
  --root . \
  --registry projects/willy_project_registry.yml \
  --out .pycsamt_rag
```

Internally, it should do:

```text
1. Scan package files
2. Exclude noise files
3. Parse Python AST
4. Extract functions/classes/docstrings/signatures
5. Parse docs/source/*.rst
6. Parse examples and tests
7. Summarize EDI line metadata
8. Build vector index
9. Build keyword BM25 index
10. Build symbol graph
11. Save manifest
```

The output:

```text
.pycsamt_rag/
│
├── manifest.json
├── symbols.jsonl
├── documents.jsonl
├── chunks.jsonl
├── data_catalog.json
├── workflow_catalog.json
├── project_registry.normalized.json
├── vector_store/
├── bm25_store/
└── graph/
    ├── symbol_edges.jsonl
    └── workflow_edges.jsonl
```

---

# 10. Recommended retrieval logic

When a user asks a question, do **not** retrieve randomly.

Use staged retrieval:

```text
User query
   ↓
IntentRouter
   ↓
Entity extraction
   - line name: L22PLT
   - workflow: static_shift
   - requested output: code
   ↓
Project registry lookup
   ↓
Hybrid retrieval
   - keyword search
   - vector search
   - symbol search
   - recipe search
   ↓
Reranking
   ↓
Context assembly
   ↓
LLM answer
```

Example:

```text
Query:
"write the static shift code for processing my EDI data of line L22PLT"
```

The assistant should classify:

```json
{
  "intent": "code",
  "workflow": "static_shift",
  "line": "L22PLT",
  "requires_data_path": true,
  "requires_execution": false,
  "requires_plotting": true
}
```

Then retrieve:

```text
assistant_recipes/static_shift.md
pycsamt/agents/static_shift.py
pycsamt/emtools/ss.py
pycsamt/agents/code_gen.py
pycsamt/agents/loader.py
projects/willy_project_registry.yml
```

Then generate code.

---

# 11. The assistant should use tools, not only RAG

RAG gives context.

But for a powerful assistant, you need tools.

Add tools like:

```python
resolve_line(line_name)
load_edi_survey(path)
run_static_shift(path, output_dir, method, half_window)
generate_static_shift_script(line_name)
plot_static_shift_summary(path, output_dir)
prepare_occam2d(path, output_dir)
run_phase_tensor_analysis(path, output_dir)
validate_generated_code(script_path)
```

So the assistant workflow becomes:

```text
RAG finds the right API
Tool executes the right function
Validation checks result
LLM explains the result
```

This is much more reliable than:

```text
LLM guesses code
```

---

# 12. Concrete implementation skeleton

## `assistant/rag/schemas.py`

```python
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class RAGChunk:
    id: str
    text: str
    source_path: str
    kind: str
    symbol: str | None = None
    workflow: str | None = None
    module: str | None = None
    priority: int = 1
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class RetrievedContext:
    query: str
    chunks: list[RAGChunk]
    symbols: list[dict[str, Any]]
    project_context: dict[str, Any]
```

---

## `assistant/rag/config.py`

```python
EXCLUDE_PATTERNS = [
    "__pycache__",
    ".pyc",
    "docs/source/_static",
    "docs/source/images",
    "pycsamt/app/desktop/resources",
    "pycsamt/app/web/assets",
    ".pytest_cache",
    ".ruff_cache",
    "build",
    "dist",
]

HIGH_PRIORITY_PATHS = [
    "pycsamt/agents",
    "pycsamt/emtools",
    "pycsamt/api",
    "pycsamt/site",
    "pycsamt/seg",
    "pycsamt/pipeline",
    "pycsamt/ai",
    "docs/source",
    "examples",
]

WORKFLOW_KEYWORDS = {
    "static_shift": [
        "static shift",
        "galvanic",
        "shift correction",
        "estimate_ss_ama",
        "correct_ss_ama",
        "StaticShiftAgent",
    ],
    "qc": [
        "quality control",
        "QC",
        "snr",
        "dead band",
        "DataQCAgent",
    ],
    "phase_analysis": [
        "phase tensor",
        "strike",
        "skew",
        "dimensionality",
    ],
    "inversion": [
        "occam",
        "modem",
        "mare2dem",
        "inversion",
    ],
}
```

---

## `assistant/rag/ast_indexer.py`

```python
from __future__ import annotations

import ast
from pathlib import Path

from .schemas import RAGChunk


def index_python_file(path: Path, root: Path) -> list[RAGChunk]:
    source = path.read_text(encoding="utf-8", errors="ignore")
    tree = ast.parse(source)
    chunks: list[RAGChunk] = []

    module = str(path.relative_to(root)).replace("/", ".").removesuffix(".py")

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            name = node.name
            doc = ast.get_docstring(node) or ""

            start = getattr(node, "lineno", 1)
            end = getattr(node, "end_lineno", start)
            lines = source.splitlines()
            code = "\n".join(lines[start - 1:end])

            symbol = f"{module}.{name}"

            text = f"""
Symbol: {symbol}
Type: {type(node).__name__}
File: {path}
Docstring:
{doc}

Code:
{code}
""".strip()

            chunks.append(
                RAGChunk(
                    id=f"{path}:{start}:{end}",
                    text=text,
                    source_path=str(path),
                    kind="python_symbol",
                    symbol=symbol,
                    module=module,
                    priority=3 if "agents" in str(path) or "emtools" in str(path) else 1,
                    metadata={
                        "start_line": start,
                        "end_line": end,
                        "name": name,
                    },
                )
            )

    return chunks
```

---

## `assistant/tools/project_registry.py`

```python
from __future__ import annotations

from pathlib import Path
import yaml


class ProjectRegistry:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.data = yaml.safe_load(self.path.read_text())

    def resolve_line(self, line_name: str) -> dict:
        lines = self.data.get("lines", {})
        if line_name not in lines:
            known = ", ".join(sorted(lines))
            raise KeyError(
                f"Unknown line {line_name!r}. Known lines: {known}"
            )

        line = dict(lines[line_name])
        project = self.data.get("project", {})
        defaults = {
            "project": project.get("name"),
            "line": line_name,
            "edi_dir": line.get("edi_dir"),
            "sort_by": line.get("sort_by", "lon"),
            "output_root": line.get(
                "output_root",
                f"results/{line_name}"
            ),
            "static_shift": self.data.get("static_shift_defaults", {}),
            "plot": self.data.get("plot_defaults", {}),
        }
        return defaults
```

---

# 13. Static-shift RAG recipe

Create:

```text
assistant_recipes/static_shift.md
```

Content:

````markdown
# Static-shift correction recipe

## User intent phrases

- static shift
- galvanic shift
- correct static shift
- remove static shift
- shift correction
- write static shift code
- process EDI line
- AMA correction

## Required inputs

- EDI directory or survey line name
- output directory
- method: ama, loess, refmedian, bilateral
- sort_by: lon, lat, name
- half_window
- optional period band pband

## Preferred high-level API

Use `StaticShiftAgent` when the user wants to run the workflow.

```python
from pycsamt.agents import StaticShiftAgent

agent = StaticShiftAgent(method="ama", half_window=3)
result = agent.execute({
    "path": edi_dir,
    "output_dir": output_dir,
})
````

## Preferred low-level API

Use `pycsamt.emtools.ss` when the user asks for standalone code.

```python
from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.ss import estimate_ss_ama, correct_ss_ama

sites = ensure_sites(edi_dir)
table = estimate_ss_ama(sites, half_window=3, sort_by="lon")
sites_corr = correct_ss_ama(sites, half_window=3, sort_by="lon")
```

## Expected outputs

* static shift factor table
* corrected Sites object
* summary figure
* per-station 1D curve figure
* optional corrected EDI files

````

This kind of recipe helps the assistant tremendously.

---

# 14. The final assistant flow for pycsamt-v2

The full assistant should work like this:

```text
User:
"write the static shift code for processing my EDI data of line L22PLT"

1. IntentRouter
   intent = code
   workflow = static_shift

2. EntityExtractor
   line = L22PLT
   requested_output = Python script

3. ProjectRegistry
   L22PLT → /real/path/to/L22PLT

4. RAG Retriever
   retrieve:
   - static_shift.md recipe
   - StaticShiftAgent
   - estimate_ss_ama
   - correct_ss_ama
   - plot_ss_summary
   - plot_ss_1d_curves
   - CodeGenerationAgent template

5. ContextBuilder
   build compact prompt with only relevant context

6. LLM
   generate package-correct code

7. Validator
   check imports
   check path exists
   optionally run dry-run

8. Final answer
   return code + explanation + output locations
````

---

# 15. Example final answer the assistant should produce

For the user query:

```text
write the static shift code for processing my EDI data of line L22PLT
```

The assistant should produce something like:

```python
from pathlib import Path

from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.ss import (
    estimate_ss_ama,
    correct_ss_ama,
    plot_ss_summary,
    plot_ss_1d_curves,
)
from pycsamt.emtools._core import _get_z_block, _iter_items, _name

import numpy as np
import matplotlib.pyplot as plt


# Resolved from project registry
line_name = "L22PLT"
edi_dir = Path("/absolute/path/to/WILLY_DATA/L22PLT")
output_dir = Path("results") / line_name / "static_shift"
output_dir.mkdir(parents=True, exist_ok=True)


# 1. Load EDI files
sites = ensure_sites(
    edi_dir,
    recursive=True,
    on_dup="replace",
    strict=False,
    verbose=1,
)

print(f"Loaded {len(list(sites))} stations from {edi_dir}")


# 2. Estimate AMA static-shift factors
ss_table = estimate_ss_ama(
    sites,
    sort_by="lon",
    half_window=3,
    weights="tri",
    pband=None,
    max_skew=6.0,
    robust_freq="median",
    robust_overall="median",
    recursive=False,
    strict=False,
    verbose=1,
)

ss_table.to_csv(
    output_dir / "static_shift_factors_L22PLT.csv",
    index=False,
)

print(ss_table[["station", "delta_log10_rho", "fac_rho", "fac_z", "n_used"]])


# 3. Apply static-shift correction
sites_corr = correct_ss_ama(
    sites,
    sort_by="lon",
    half_window=3,
    weights="tri",
    pband=None,
    max_skew=6.0,
    inplace=False,
    verbose=1,
)


# 4. Collect log10 apparent resistivity before/after

def collect_logrho_det(S):
    rows = []
    freqs = None
    labels = []

    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue

        rxy = 0.2 * np.abs(z[:, 0, 1]) ** 2 / (fr + 1e-24)
        ryx = 0.2 * np.abs(z[:, 1, 0]) ** 2 / (fr + 1e-24)

        logrho = np.log10(np.sqrt(rxy * ryx) + 1e-24)

        rows.append(logrho)
        freqs = fr
        labels.append(_name(ed, i))

    return np.asarray(rows), freqs, labels


logrho_before, freqs, labels = collect_logrho_det(sites)
logrho_after, _, _ = collect_logrho_det(sites_corr)


# 5. Generate static-shift plots
fig_summary = plot_ss_summary(
    logrho_before,
    logrho_after,
    freqs=freqs,
    station_labels=labels,
)

fig_summary.savefig(
    output_dir / "static_shift_summary_L22PLT.png",
    dpi=150,
    bbox_inches="tight",
)

fig_curves = plot_ss_1d_curves(
    logrho_before,
    logrho_after,
    freqs=freqs,
    station_labels=labels,
)

fig_curves.savefig(
    output_dir / "static_shift_1d_curves_L22PLT.png",
    dpi=150,
    bbox_inches="tight",
)

plt.close("all")

print(f"Static-shift correction completed for {line_name}")
print(f"Outputs saved in: {output_dir}")
```

This is the type of response your future assistant should generate automatically.

---

# 16. Recommended model strategy

Use different models for different stages.

| Stage                  | Recommended model type     |
| ---------------------- | -------------------------- |
| Intent routing         | cheap/fast model           |
| Query rewriting        | cheap/fast model           |
| Retrieval reranking    | embedding/reranker model   |
| Code generation        | strongest model            |
| Scientific explanation | strongest model            |
| Validation             | deterministic Python tools |
| Plotting               | pycsamt + Python execution |

You can support:

```text
OpenAI
Gemini
Anthropic
DeepSeek
```

but the important point is to hide them behind one internal interface:

```python
class LLMProvider:
    def complete(...)
    def embed(...)
    def rerank(...)
```

The user should not care which provider is used.

---

# 17. Evaluation set you must create

To know whether the assistant works, create:

```text
assistant/evals/
│
├── rag_questions_static_shift.jsonl
├── rag_questions_loading.jsonl
├── rag_questions_plotting.jsonl
├── rag_questions_inversion.jsonl
├── rag_questions_package_qa.jsonl
└── expected_tools.jsonl
```

Example:

```json
{
  "query": "write the static shift code for processing my EDI data of line L22PLT",
  "expected_intent": "code",
  "expected_workflow": "static_shift",
  "expected_retrieved_symbols": [
    "pycsamt.emtools.ss.estimate_ss_ama",
    "pycsamt.emtools.ss.correct_ss_ama",
    "pycsamt.emtools.ss.plot_ss_summary",
    "pycsamt.emtools.ss.plot_ss_1d_curves"
  ],
  "must_not_contain": [
    "from pycsamt import static_shift",
    "fake_static_shift_function"
  ]
}
```

Then run:

```bash
pycsamt-rag eval --suite assistant/evals/rag_questions_static_shift.jsonl
```

Score:

```text
intent accuracy
workflow accuracy
retrieval recall
symbol recall
code import validity
script execution validity
plot output validity
hallucination rate
```

---

# 18. Priority roadmap

## Phase 1 — Minimum useful RAG

Do this first:

```text
1. Add root README.md
2. Add AGENTS.md
3. Add assistant_recipes/static_shift.md
4. Add project registry for L18PLT and L22PLT
5. Build AST index for pycsamt/agents and pycsamt/emtools
6. Build docs index for docs/source
7. Add hybrid retriever
8. Connect RAG to PackageQAAgent
```

## Phase 2 — Workflow-aware assistant

```text
1. Add resolve_line tool
2. Add generate_static_shift_script tool
3. Add run_static_shift tool
4. Add plot tools
5. Add validation tools
6. Connect IntentRouter → RAG → Tool execution
```

## Phase 3 — Full production assistant

```text
1. Add session memory
2. Add user project memory
3. Add notebook generation
4. Add report generation
5. Add eval suite
6. Add trace logging
7. Add web/desktop integration
```

---

# 19. The key improvement for your current package

Your package already has:

```text
IntentRouter
PackageQAAgent
CodeGenerationAgent
StaticShiftAgent
WorkflowOrchestratorAgent
Pipeline registry
Docs
Tests
Examples
```

So now you should **not rewrite everything**.

Instead, add this missing layer:

```text
pycsamt.assistant.rag
```

and connect it to:

```text
IntentRouter
PackageQAAgent
CodeGenerationAgent
WorkflowOrchestratorAgent
```

The final architecture should be:

```text
User request
   ↓
IntentRouter
   ↓
Entity extractor
   ↓
Project registry
   ↓
RAG retriever
   ↓
Context builder
   ↓
LLM provider
   ↓
pycsamt tool execution
   ↓
validation
   ↓
final answer/code/plots/report
```

That is the structure that will make your pycsamt-v2 assistant powerful, reliable, and useful for real geophysical workflows.

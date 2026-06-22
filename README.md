<div align="center">
  <img src="docs/source/_static/pycsamt-v2-symbol-logo.svg" alt="pyCSAMT v2 logo" width="340"/>
</div>

<p align="center">
  <strong>Scientific Python toolkit for MT, AMT, CSAMT, CSEM, and TDEM survey processing, inversion, interpretation, and reproducible workflow automation.</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/pycsamt/"><img alt="PyPI" src="https://img.shields.io/pypi/v/pycsamt?color=orange"/></a>
  <a href="https://pycsamt.readthedocs.io"><img alt="Docs" src="https://img.shields.io/readthedocs/pycsamt?label=docs"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/actions"><img alt="CI" src="https://img.shields.io/github/actions/workflow/status/earthai-tech/pycsamt/ci.yml?label=CI"/></a>
  <a href="https://opensource.org/licenses/LGPL-3.0"><img alt="License" src="https://img.shields.io/badge/license-LGPL--3.0-blue"/></a>
  <a href="https://www.python.org/downloads/"><img alt="Python" src="https://img.shields.io/pypi/pyversions/pycsamt"/></a>
</p>

---

## Status

pyCSAMT v2 is a rewrite of the original pyCSAMT project. The package is being organized as a scientific Python library with:

- a command-line interface for field and production workflows,
- a declarative pipeline system for reproducible processing,
- physics-based inversion support for electromagnetic methods,
- AI-assisted agents for workflow guidance, QC, inversion preparation, reporting, and automation,

The v2 documentation and API are under active development. Some examples may evolve before the stable v2 release.

---

## What pyCSAMT Does

pyCSAMT supports electromagnetic survey workflows from input data to interpretation-ready deliverables.

| Area | What pyCSAMT v2 provides |
| --- | --- |
| Survey data | EDI, AVG, Jones J, Zonge, TDEM/TEMAVG, station metadata, survey context |
| Data structures | SEG-style EDI objects, site collections, impedance/tipper tensors, metadata models |
| Processing | frequency filtering, tensor rotation, static-shift correction, phase analysis, quality control |
| Forward modelling | 1-D and 2-D EM forward modelling utilities for MT, CSAMT, and TDEM workflows |
| Inversion | preparation, execution wrappers, result inspection, and plotting for supported inversion backends |
| Interpretation | resistivity classification, stratigraphic products, map/profile outputs, export workflows |
| Pipelines | YAML/JSON/Python workflow definitions, presets, step registry, output manifests |
| AI and agents | LLM-assisted workflow orchestration, QC, inversion support, model-zoo access, code/report generation |
| Applications | command-line tools, desktop app entry point, and web/GUI-facing package structure |

---

## Installation

Install the base package:

```bash
pip install pycsamt
```

Install from source for v2 development:

```bash
git clone https://github.com/earthai-tech/pycsamt.git
cd pycsamt
pip install -e ".[dev,docs]"
```

Optional extras are split by workflow:

```bash
# Machine-learning backend
pip install "pycsamt[torch]"
pip install "pycsamt[tensorflow]"

# Geospatial tools
pip install "pycsamt[geo]"

# Documentation build tools
pip install "pycsamt[docs]"

# Desktop and web application dependencies
pip install "pycsamt[app]"

# Everything used by the full v2 development environment
pip install "pycsamt[full]"
```

pyCSAMT v2 currently targets Python 3.9 or newer.

---

## Quick Start

### Command line

The main executable is `pycsamt`:

```bash
pycsamt --help
pycsamt info
```

Common command groups include:

```bash
pycsamt survey set data/edi/
pycsamt edi info data/edi/
pycsamt site info
pycsamt transform avg survey.avg --output-dir outputs/edi/
pycsamt invert build data/edi/ --solver occam2d --workdir runs/occam2d/
pycsamt pipe steps
pycsamt pipe presets
```

The desktop application entry point is:

```bash
pycsamt-gui
```

### Python API

Use the Python API when you want direct control from notebooks, scripts, or applications:

```python
from pycsamt.pipeline import Pipeline, Step

pipe = Pipeline(
    [
        ("notch", Step("NR001", mains_hz=50)),
        ("band", Step("FREQ001")),
        ("align", Step("FREQ004")),
        ("static_shift", Step("SS001")),
        ("rotate", Step("TZ001")),
    ]
)

result = pipe.run(sites, outdir="outputs/profile_l22/")
print(result.summary())
```

For configuration-driven work:

```python
from pycsamt.pipeline import Pipeline

pipe = Pipeline.from_yaml("pipeline.yaml")
result = pipe.run(sites, outdir="outputs/run01/")
```

---

## Pipeline System

The v2 pipeline layer is designed for reproducible survey processing. A pipeline can be created from code, from a preset, or from a configuration file.

```bash
pycsamt pipe init --preset publication_ready --output pipeline.yaml
pycsamt pipe show pipeline.yaml
pycsamt pipe run --config pipeline.yaml --survey data/edi/ --out outputs/run01/
```

Pipeline documentation is organized around:

- concepts: what a pipeline, step, preset, and result mean,
- configuration files: YAML/JSON/Python workflow definitions,
- CLI usage: `pycsamt pipe ...`,
- step catalogue: registered processing operations,
- presets: reusable workflow templates,
- outputs: manifests, figures, tables, logs, and reports.

See `docs/source/pipeline/` for the in-progress v2 pages.

---

## AI-Assisted Agents

The `pycsamt.agents` package provides AI-assisted workflow components. Agents are lazy-loaded so the base package remains usable without LLM client libraries.

Examples of available agent groups:

- foundation agents: request parsing, loading, coordination, and orchestration,
- processing agents: QC, denoising, tensor rotation, static shift, phase/tipper analysis,
- inversion agents: Occam2D, ModEM, MARE2DEM, inversion preparation, evaluation, comparison,
- AI model agents: 1-D/2-D/3-D learned inversion, ensembles, joint inversion, model zoo,
- output agents: interpretation, resistivity maps, report generation, reproducible code generation,
- pipeline agents: pipeline construction and batch survey execution.

Minimal example without an LLM key:

```python
from pycsamt.agents import AgentCoordinator, ContextInputAgent, MTLoaderAgent

coordinator = AgentCoordinator("survey_qc")
coordinator.add_step("context", ContextInputAgent())
coordinator.add_step(
    "load",
    MTLoaderAgent(),
    input_fn=lambda results: {
        "path": results["context"]["config"]["data_path"],
    },
)

result = coordinator.execute(
    {"request": "Load data/edi, run basic QC, and prepare a report."},
)
```

Optional LLM providers include Anthropic Claude, OpenAI, and Google Gemini when their client libraries and API keys are installed/configured.

---

## Documentation

The v2 documentation is being built as a scientific-library documentation set rather than a single long manual.

Local documentation sources live in `docs/source/`:

| Section | Purpose |
| --- | --- |
| `getting_started/` | installation, configuration, data formats, first survey |
| `tutorials/` | step-by-step practical recipes |
| `user_guide/` | processing, inversion, AI inversion, interpretation |
| `pipeline/` | reproducible workflow system |
| `agents/` | AI-assisted agents and orchestration |
| `cli/` | command-line reference |
| `api/` | generated and curated API reference |
| `theory/` | scientific background |
| `development/` | API policy, docstring style, docs build process |
| `release_notes/` | version-specific release notes |

Build the docs locally:

```bash
cd docs
sphinx-build -b html source _build/html
```

The public documentation is planned at <https://pycsamt.readthedocs.io>.

---

## Development

Install development and documentation dependencies:

```bash
pip install -e ".[dev,docs]"
```

Run tests:

```bash
pytest
```

Build documentation:

```bash
cd docs
sphinx-build -b html source _build/html
```

Developer-facing documentation is being written in `docs/source/development/`, especially:

- API policy and package structure,
- NumPy-style docstring conventions,
- Sphinx documentation build workflow.

---

## Citation

If you use pyCSAMT in research or professional studies, please cite the relevant software and method papers:

```bibtex
@article{Kouadio2022,
  author  = {Kouadio, K. Laurent and Liu, Rong and Liu, Chum-ning and
             Mi, Binbin and Malory, Albert O.},
  title   = {pyCSAMT: An open-source Python library for controlled source
             audio-frequency magnetotelluric data processing and
             pseudostratigraphic log generation},
  journal = {Journal of Applied Geophysics},
  year    = {2022},
  volume  = {201},
  pages   = {104647},
  doi     = {10.1016/j.jappgeo.2022.104647}
}

@article{Kouadio2023,
  author  = {Kouadio, K. Laurent and others},
  title   = {Recovering the electrical resistivity from MT data using a
             deep resistivity convolutional neural network},
  journal = {Journal of Geophysical Research: Solid Earth},
  year    = {2023},
  doi     = {10.1029/2023JB027538}
}
```

---

## Contributing

Contributions are welcome while v2 is being stabilized. Start with:

- `docs/source/development/api_policy.rst` for public API rules,
- `docs/source/development/docstring_style.rst` for documentation style,
- `docs/source/development/documentation_build.rst` for docs build guidance,
- `docs/source/contributing.rst` for project contribution notes.

Bug reports and feature requests can be opened at <https://github.com/earthai-tech/pycsamt/issues>.

---

## License

pyCSAMT is distributed under the GNU Lesser General Public License v3.0 or later. See `LICENSE.md` for the full license text.

---

<p align="center">
  Developed by <a href="https://github.com/earthai-tech">earthai-tech</a> &nbsp;|&nbsp;
  Lead: <a href="mailto:etanoyau@gmail.com">Kouadio K. Laurent</a>
</p>

<div align="center">
  <img src="docs/source/_static/pycsamt-v2-symbol-logo.svg" alt="pyCSAMT v2" width="380"/>
  <br/>
  <em>Scientific Python for electromagnetic geophysics &mdash; processing, inversion, AI agents, and apps.</em>
  <br/><br/>
  <a href="https://pypi.org/project/pycsamt/"><img alt="PyPI" src="https://img.shields.io/pypi/v/pycsamt?color=orange&label=PyPI"/></a>
  <a href="https://pycsamt.readthedocs.io"><img alt="Docs" src="https://img.shields.io/readthedocs/pycsamt?label=docs"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/actions"><img alt="CI" src="https://img.shields.io/github/actions/workflow/status/earthai-tech/pycsamt/ci.yml?label=CI"/></a>
  <a href="https://codecov.io/gh/earthai-tech/pycsamt"><img alt="Coverage" src="https://img.shields.io/codecov/c/github/earthai-tech/pycsamt?label=coverage"/></a>
  <a href="https://opensource.org/licenses/LGPL-3.0"><img alt="License" src="https://img.shields.io/badge/license-LGPL--3.0-blue"/></a>
  <img alt="Semver" src="https://img.shields.io/badge/semver-2.0.0-informational"/>
  <br/>
  <img alt="Python" src="https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue?logo=python&logoColor=white"/>
  <img alt="Backends" src="https://img.shields.io/badge/backend-PyTorch%20%7C%20TensorFlow-EE4C2C?logo=pytorch&logoColor=white"/>
  <img alt="Platform" src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey"/>
  <br/>
  <a href="https://pypi.org/project/pycsamt/"><img alt="Downloads" src="https://img.shields.io/pypi/dm/pycsamt?label=downloads&color=brightgreen"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/commits/v2"><img alt="Last commit" src="https://img.shields.io/github/last-commit/earthai-tech/pycsamt/v2"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/issues"><img alt="Issues" src="https://img.shields.io/github/issues/earthai-tech/pycsamt"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/pulls"><img alt="PRs welcome" src="https://img.shields.io/badge/PRs-welcome-brightgreen"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/stargazers"><img alt="Stars" src="https://img.shields.io/github/stars/earthai-tech/pycsamt?style=flat&color=yellow"/></a>
  <img alt="Code size" src="https://img.shields.io/github/languages/code-size/earthai-tech/pycsamt"/>
</div>

---

## Overview

pyCSAMT v2 is a ground-up rewrite of the original pyCSAMT library, designed as a
full-lifecycle scientific toolkit for controlled-source and natural-source
electromagnetic (EM) geophysics. It covers every stage of an EM survey project:
reading multi-format field data, applying industry-standard processing corrections,
driving 1-D through 3-D inversions, and producing publication-ready outputs — all
from a coherent, scikit-learn-inspired Python API.

What sets v2 apart is its integration of **physics-informed neural networks (PINNs)**
and **hybrid deep-learning inverters** alongside classical inversion backends
(Occam2D, ModEM, MARE2DEM), a declarative **pipeline system** for reproducible
workflows, and a suite of **AI-assisted agents** that can orchestrate entire survey
projects through natural-language instructions. An interactive **web dashboard**
(Dash) and a **desktop GUI** bring these capabilities to users who prefer a
point-and-click experience.

---

## What pyCSAMT v2 Does

| Layer | Capabilities |
|---|---|
| **Data I/O** | EDI, AVG, Jones J, Zonge, TDEM/TEMAVG, MARE2DEM `.emdata`, `.resistivity`, `.poly` |
| **Processing** | Notch filtering, static-shift correction, tensor rotation, phase-tensor analysis, quality control |
| **Forward modelling** | 1-D and 2-D EM forward utilities for MT, CSAMT, and TDEM survey design |
| **Classical inversion** | Occam2D, ModEM, MARE2DEM wrappers — input builders, runners, result loaders, section plots |
| **AI inversion** | PINN 1-D / 2-D / 3-D, hybrid deep-learning inverters, model-zoo access |
| **Interpretation** | Resistivity classification, pseudostratigraphic logs, profile and map outputs |
| **Pipeline system** | YAML / JSON / Python workflow definitions, step registry, presets, reproducible manifests |
| **AI agents** | LLM-driven orchestration for QC, inversion prep, report and code generation (Anthropic, OpenAI, Gemini) |
| **Applications** | CLI (`pycsamt`), interactive web dashboard (Dash), desktop GUI |

---

## Installation

The base package installs with a single command and has minimal dependencies:

```bash
pip install pycsamt
```

For the complete v2 environment — ML backends, apps, geospatial tools, and docs:

```bash
pip install "pycsamt[full]"
```

Install from source to track the active development branch:

```bash
git clone https://github.com/earthai-tech/pycsamt.git
cd pycsamt
git checkout v2
pip install -e ".[full]"
```

Individual extras let you pull in only what your workflow needs:

| Extra | Installs |
|---|---|
| `torch` | PyTorch backend for PINN / hybrid inverters |
| `tensorflow` | TensorFlow / Keras backend |
| `geo` | GeoPandas, Fiona, pyproj for spatial workflows |
| `app` | Dash, Flask, and desktop GUI dependencies |
| `docs` | Sphinx, PyData theme, numpydoc |
| `full` | All of the above |

Requires **Python 3.9 or newer**.

---

## Quick Start

### Processing pipeline

Build a repeatable processing chain from named steps and run it against a site
collection in a single call. The result object carries per-step outputs, a
processing log, and an exportable manifest.

```python
from pycsamt.pipeline import Pipeline, Step

pipe = Pipeline([
    ("notch",        Step("NR001", mains_hz=50)),
    ("band",         Step("FREQ001")),
    ("static_shift", Step("SS001")),
    ("rotate",       Step("TZ001")),
])
result = pipe.run(sites, outdir="outputs/run01/")
print(result.summary())
```

Pipelines can also be defined in YAML and executed from the CLI:

```bash
pycsamt pipe run --config pipeline.yaml --survey data/edi/ --out outputs/run01/
```

### MARE2DEM inversion results

Load an inversion run directory and plot observed versus predicted MT responses
or the resistivity section — the scanner finds all output files automatically:

```python
from pycsamt.models.mare2dem import InversionResult, PlotResponse, PlotModel

result  = InversionResult("runs/demo_mt/")
PlotResponse(result).plot(max_rx=6, savefig="response.pdf")
PlotModel(result).plot(cmap="turbo_r", savefig="section.pdf")
```

### AI-assisted agent

Give a plain-language task to the agent master and receive a structured result —
the orchestrator routes subtasks to specialist agents for loading, QC, inversion
setup, and reporting:

```python
from pycsamt.agents import AgentMaster

master = AgentMaster(provider="anthropic")
report = master.run(
    "Load data/edi/, flag stations with RMS > 2, build an Occam2D input "
    "for profile L22, launch inversion, and produce a PDF report."
)
```

### PINN / hybrid inverter

Physics-informed inversion requires only the site data and a forward operator:

```python
from pycsamt.ai.inversion import PINN2D

inv = PINN2D(n_layers=64, epochs=3000, backend="torch")
model = inv.fit(sites, frequencies=freqs)
model.plot_section()
```

---

## Applications

pyCSAMT v2 ships with three ready-to-use application entry points.

**Command-line interface** — covers the full workflow from the terminal:

```bash
pycsamt survey set data/edi/
pycsamt edi info data/edi/
pycsamt invert build data/edi/ --solver occam2d --workdir runs/occam2d/
pycsamt pipe presets
```

**Web dashboard** — an interactive Dash application for data exploration,
pipeline configuration, agent chat, and inversion monitoring:

```bash
pycsamt-web
```

**Desktop GUI** — a native desktop application for point-and-click access
to the same processing and inversion capabilities:

```bash
pycsamt-gui
```

---

## Citation

If pyCSAMT contributes to published research, please cite the relevant works:

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

Bug reports and feature requests are welcome at the
[issue tracker](https://github.com/earthai-tech/pycsamt/issues).
For code contributions, please read the developer guide in
`docs/source/development/` before opening a pull request.

---

## License

pyCSAMT is distributed under the
[GNU Lesser General Public License v3.0](https://opensource.org/licenses/LGPL-3.0)
or later. See [`LICENSE.md`](LICENSE.md) for the full text.

---

<p align="center">
  Developed by <a href="https://github.com/earthai-tech">earthai-tech</a>
  &nbsp;&mdash;&nbsp;
  Lead developer: <a href="mailto:etanoyau@gmail.com">Kouadio K. Laurent</a>
</p>

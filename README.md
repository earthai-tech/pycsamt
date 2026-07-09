<div align="center">
  <img src="docs/source/_static/logo/pycsamt-v2-symbol-logo.svg" alt="pyCSAMT v2" width="380"/>
  <br/>
  <em>Scientific Python for electromagnetic geophysics &mdash; processing, inversion, AI agents, and apps.</em>
  <br/><br/>
  <a href="https://pypi.org/project/pycsamt/"><img alt="PyPI" src="https://img.shields.io/pypi/v/pycsamt?color=orange&label=PyPI"/></a>
  <a href="https://pycsamt.readthedocs.io/en/latest/?badge=latest"><img alt="Docs" src="https://readthedocs.org/projects/pycsamt/badge/?version=latest"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/actions/workflows/ci.yml?query=branch%3Av2"><img alt="CI" src="https://github.com/earthai-tech/pycsamt/actions/workflows/ci.yml/badge.svg?branch=v2"/></a>
  <a href="https://codecov.io/github/earthai-tech/pycsamt?branch=v2"><img alt="Coverage" src="https://codecov.io/github/earthai-tech/pycsamt/graph/badge.svg?branch=v2"/></a>
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

**pyCSAMT v2** is a full-lifecycle toolkit for controlled-source and
natural-source electromagnetic (EM) geophysics — CSAMT, AMT, MT, and TDEM.
One coherent, scikit-learn-inspired API takes a survey from raw field files
through quality control, corrections, and inversion to interpreted,
publication-ready results.

<div align="center">

**[Documentation](https://pycsamt.readthedocs.io)** ·
[Getting started](https://pycsamt.readthedocs.io/en/latest/getting_started/index.html) ·
[User guide](https://pycsamt.readthedocs.io/en/latest/user_guide/index.html) ·
[Examples](https://pycsamt.readthedocs.io/en/latest/examples/index.html) ·
[API reference](https://pycsamt.readthedocs.io/en/latest/api/index.html)

</div>

## ✨ Highlights

- 📥 **Data I/O & QC** — EDI, Zonge AVG, Jones J, TDEM, and MARE2DEM files in
  one site model; frequency audits and noisy-station flagging.
- 📡 **IoT-enabled field acquisition** — station telemetry with edge QC
  (powerline harmonics, SNR, contact resistance, frequency coverage),
  clock-synchronisation audit, power monitoring, pluggable transports
  (file/HTTP built in; MQTT/serial/WebSocket optional), acquisition
  provenance, a field-network simulator, and dashboards — feeding directly
  into the processing pipeline.
- 🎚️ **Processing & corrections** — a catalogue of 25 methods in six
  categories: notch filtering, static-shift removal, tensor rotation,
  phase-tensor analysis, and more.
- 🧱 **Forward modelling** — synthetic layered-earth and 2-D models, forward
  responses, realistic noise, and datasets for survey design or training.
- 🧠 **Inversion, classical & AI** — Occam2D, ModEM, and MARE2DEM end to end,
  plus physics-informed neural networks (PINN 1-D/2-D/3-D) and hybrid
  deep-learning inverters.
- 🗺️ **Interpretation & mapping** — resistivity classification,
  pseudostratigraphic logs, station maps, pseudosections, and 3-D quick looks.
- 🤖 **Pipelines, agents & apps** — reproducible YAML/JSON/Python workflows,
  LLM-driven agents (Anthropic, OpenAI, Gemini), a web dashboard, and a
  desktop GUI.

## 🚀 Installation

```bash
pip install pycsamt           # core, minimal dependencies
pip install "pycsamt[full]"   # + ML backends, apps, geospatial, docs
```

Requires **Python 3.9+**.

<details>
<summary><b>Optional extras &amp; source install</b></summary>

| Extra | Installs |
|---|---|
| `torch` | PyTorch backend for PINN / hybrid inverters |
| `tensorflow` | TensorFlow / Keras backend |
| `geo` | pyproj, xarray, contextily for maps and reprojection |
| `agents` | Anthropic, OpenAI/DeepSeek, and Gemini SDKs for AI agents |
| `app` | Desktop GUI (PySide6) + Dash web dashboard |
| `docs` | Sphinx, PyData theme, numpydoc |
| `full` | All of the above (prefers the PyTorch backend) |

Track the active development branch:

```bash
git clone https://github.com/earthai-tech/pycsamt.git
cd pycsamt && git checkout v2
pip install -e ".[full]"
```

</details>

## ⚡ Quick start

Chain named processing steps into a repeatable pipeline — the result carries
per-step outputs, a log, and an exportable manifest:

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

Or hand the whole job to an AI agent in plain language:

```python
from pycsamt.agents import AgentMaster

master = AgentMaster(provider="anthropic")
report = master.run(
    "Load data/edi/, flag stations with RMS > 2, build an Occam2D input "
    "for profile L22, launch inversion, and produce a PDF report."
)
```

<details>
<summary><b>More examples — inversion results, PINN, IoT, CLI</b></summary>

Load a MARE2DEM run directory and plot responses or the resistivity section:

```python
from pycsamt.models.mare2dem import InversionResult, PlotResponse, PlotModel

result = InversionResult("runs/demo_mt/")
PlotResponse(result).plot(max_rx=6, savefig="response.pdf")
PlotModel(result).plot(cmap="turbo_r", savefig="section.pdf")
```

Physics-informed inversion from site data and a forward operator:

```python
from pycsamt.ai.inversion import PINN2D

inv = PINN2D(n_layers=64, epochs=3000, backend="torch")
model = inv.fit(sites, frequencies=freqs)
model.plot_section()
```

Ingest IoT field telemetry, quality-control it at the edge, and hand it
straight to processing — with a reproducible provenance manifest:

```python
from pycsamt.iot import FieldSession, simulate_iot_network, plot_field_dashboard

# real edge telemetry, or a simulated network for demos/tests
packets = simulate_iot_network(n_stations=24, profiles=["L1", "L3"])
session = FieldSession("SSL2026")
session.add_packets(packets)

session.assess()                                 # stream QC -> MonitoringStatus
session.to_pipeline_input()                      # hand-off for processing
session.export_manifest("field_manifest.json")   # reproducible provenance
plot_field_dashboard(session, output_path="dashboard.png")
```

The full workflow is also scriptable from the terminal:

```bash
pycsamt survey set data/edi/
pycsamt invert build data/edi/ --solver occam2d --workdir runs/occam2d/
pycsamt pipe run --config pipeline.yaml --survey data/edi/ --out outputs/run01/
```

</details>

## 🖥️ Interfaces

| Python API | CLI | Web dashboard | Desktop GUI |
|:---:|:---:|:---:|:---:|
| `import pycsamt` | `pycsamt --help` | `pycsamt-web` | `pycsamt-desktop` |

The same engine drives all four — script it, automate it, or point and click.

## 📖 Citation

If pyCSAMT contributes to published research, please cite
[Kouadio et al. (2022), *J. Applied Geophysics*](https://doi.org/10.1016/j.jappgeo.2022.104647).

<details>
<summary><b>BibTeX</b></summary>

```bibtex
@article{Kouadio2022,
  author  = {Kouadio, K. L. and Liu, R. and Mi, B. and Liu, C.},
  title   = {pyCSAMT: An alternative Python toolbox for groundwater exploration
             using controlled-source audio-frequency magnetotelluric},
  journal = {Journal of Applied Geophysics},
  year    = {2022},
  doi     = {10.1016/j.jappgeo.2022.104647}
}

@article{Kouadio2023,
  author  = {Kouadio, K. L. and Liu, R. and Malory, A. O. and Liu, W. and Liu, C.},
  title   = {A novel approach for water reservoir mapping using
             controlled-source audio-frequency magnetotelluric in Xingning area,
             Hunan Province, China},
  journal = {Geophysical Prospecting},
  year    = {2023},
  doi     = {10.1111/1365-2478.13385}
}
```

</details>

## 🤝 Contributing & license

Bug reports and feature requests are welcome on the
[issue tracker](https://github.com/earthai-tech/pycsamt/issues); see the
[developer guide](https://pycsamt.readthedocs.io/en/latest/development/index.html)
before opening a pull request. Distributed under the
[LGPL-3.0](https://opensource.org/licenses/LGPL-3.0) or later — see
[`LICENSE.md`](LICENSE.md).

---

<p align="center">
  Developed by <a href="https://github.com/earthai-tech">earthai-tech</a>
  &nbsp;&mdash;&nbsp;
  Lead developer: <a href="mailto:etanoyau@gmail.com">Kouadio K. Laurent</a>
</p>

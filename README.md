<div align="center">
  <img src="docs/source/_static/pycsamt-v2-symbol-logo.svg" alt="pyCSAMT v2" width="320"/>
</div>

<p align="center">
  <strong>Python toolkit for MT, AMT, CSAMT, CSEM &amp; TDEM — from raw survey data to AI-assisted inversion and interpretation.</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/pycsamt/"><img alt="PyPI" src="https://img.shields.io/pypi/v/pycsamt?color=orange"/></a>
  <a href="https://pycsamt.readthedocs.io"><img alt="Docs" src="https://img.shields.io/readthedocs/pycsamt?label=docs"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/actions"><img alt="CI" src="https://img.shields.io/github/actions/workflow/status/earthai-tech/pycsamt/ci.yml?label=CI"/></a>
  <a href="https://opensource.org/licenses/LGPL-3.0"><img alt="License" src="https://img.shields.io/badge/license-LGPL--3.0-blue"/></a>
  <a href="https://www.python.org/downloads/"><img alt="Python" src="https://img.shields.io/pypi/pyversions/pycsamt"/></a>
</p>

---

pyCSAMT v2 is a complete rewrite built around a clean scientific-Python API.
It covers the full electromagnetic survey lifecycle — data loading, processing,
forward modelling, inversion, and interpretation — and adds a declarative
pipeline layer and AI-assisted agents powered by large language models.

**Key capabilities**

| | |
|---|---|
| **Data** | EDI, AVG, Jones J, Zonge, TDEM/TEMAVG, MARE2DEM `.emdata` |
| **Processing** | notch filter, static-shift correction, tensor rotation, phase analysis, QC |
| **Inversion** | Occam2D, ModEM, MARE2DEM wrappers; PINN and hybrid deep-learning inverters |
| **Agents** | LLM-driven orchestration for QC, inversion prep, code and report generation |
| **Pipeline** | YAML/JSON/Python workflow definitions with presets and reproducible outputs |
| **Apps** | CLI (`pycsamt`), web dashboard (Dash), desktop GUI |

---

## Install

```bash
pip install pycsamt          # base
pip install "pycsamt[full]"  # everything
```

From source (v2 development):

```bash
git clone https://github.com/earthai-tech/pycsamt.git
cd pycsamt
pip install -e ".[full]"
```

Optional extras: `torch`, `tensorflow`, `geo`, `app`, `docs`.
Requires Python 3.9+.

---

## Quick look

```python
# Reproducible processing pipeline
from pycsamt.pipeline import Pipeline, Step

pipe = Pipeline([
    ("notch",        Step("NR001", mains_hz=50)),
    ("band",         Step("FREQ001")),
    ("static_shift", Step("SS001")),
    ("rotate",       Step("TZ001")),
])
result = pipe.run(sites, outdir="outputs/run01/")

# AI-assisted inversion agent
from pycsamt.agents import AgentMaster

master = AgentMaster(provider="anthropic")
report = master.run(
    "Load data/edi, flag noisy stations, build an Occam2D input, "
    "run inversion, and return a resistivity section plot."
)

# MARE2DEM results in three lines
from pycsamt.models.mare2dem import InversionResult, PlotResponse

result = InversionResult("runs/demo_mt/")
PlotResponse(result).plot(max_rx=4, savefig="response.pdf")
```

---

## Citation

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
```

---

## License

LGPL-3.0 — see [`LICENSE.md`](LICENSE.md).

---

<p align="center">
  Developed by <a href="https://github.com/earthai-tech">earthai-tech</a> &nbsp;|&nbsp;
  Lead: <a href="mailto:etanoyau@gmail.com">Kouadio K. Laurent</a>
</p>

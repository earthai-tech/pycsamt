<div align="center">
  <img src="docs/source/_static/pycsamt_logo.svg" alt="pyCSAMT logo" width="320"/>
</div>

<h1 align="center">pyCSAMT v2</h1>

<p align="center">
  <strong>Python toolkit for audio-frequency magnetotelluric (MT / AMT / CSAMT / CSEM) data processing, inversion, and geological interpretation</strong>
</p>

<p align="center">
  <a href="https://pypi.org/project/pycsamt/"><img alt="PyPI" src="https://img.shields.io/pypi/v/pycsamt?color=orange"/></a>
  <a href="https://pycsamt.readthedocs.io"><img alt="Docs" src="https://img.shields.io/readthedocs/pycsamt?label=docs"/></a>
  <a href="https://github.com/earthai-tech/pycsamt/actions"><img alt="CI" src="https://img.shields.io/github/actions/workflow/status/earthai-tech/pycsamt/ci.yml?label=CI"/></a>
  <a href="https://opensource.org/licenses/LGPL-3.0"><img alt="License" src="https://img.shields.io/badge/license-LGPL--3.0-blue"/></a>
  <a href="https://www.python.org/downloads/"><img alt="Python" src="https://img.shields.io/pypi/pyversions/pycsamt"/></a>
</p>

---

## Overview

**pyCSAMT v2** is a complete rewrite of the original pyCSAMT library, designed to support the full lifecycle of electromagnetic (EM) survey data from field acquisition to geological interpretation.

| Scope | Supported methods |
|-------|------------------|
| Controlled-source AMT | CSAMT, CSEM |
| Natural-source MT | AMT, broadband MT |
| Inversion targets | 1-D, 2-D resistivity models |
| AI-assisted inversion | Physics-informed CNN (DRCNN), joint multi-survey |

### Key modules

| Package | Purpose |
|---------|---------|
| `pycsamt.io` | EDI / AVG / J read-write; `EDICollection` batch loader |
| `pycsamt.processing` | Static-shift correction, phase tensor, Z decomposition |
| `pycsamt.models` | OCCAM2D and ModEM I/O; forward modelling helpers |
| `pycsamt.inversion` | Inversion wrappers (`Occam2DInv`, `ModEMInv`) |
| `pycsamt.forward` | 1-D / 2-D analytic and finite-difference forward solvers |
| `pycsamt.ai` | Dual-backend (PyTorch / TensorFlow) deep-learning inversion |
| `pycsamt.interp` | Post-inversion interpretation: stratigraphic logs, calibration, export |

---

## Quick start

```python
import pycsamt

# --- Load a set of EDI files ---
from pycsamt.io import EDICollection
coll = EDICollection.from_dir("data/edi/")
print(coll)          # EDICollection(n_sites=47, freq_range=[0.001, 10000] Hz)

# --- Run OCCAM2D inversion ---
from pycsamt.models import OccamModel
model = OccamModel(n_layers=60, target_rms=1.0)
result = model.fit(coll)

# --- Interpret the 2-D resistivity model ---
from pycsamt.interp import ResistivityModel, ModelCalibrator
rm = ResistivityModel.from_occam2d(result)
cal = ModelCalibrator(ptol=0.30, verbose=True).fit(rm)

# Export pseudostratigraphic logs for Oasis Montaj
from pycsamt.interp import export
export.to_oasis_montaj_xyz(cal.stratigraphic_logs(), "output/profile.xyz")
```

---

## Installation

### Stable release (PyPI)

```bash
pip install pycsamt
```

### With deep-learning support

```bash
# PyTorch backend
pip install "pycsamt[torch]"

# TensorFlow backend
pip install "pycsamt[tensorflow]"

# Both backends + geospatial extras
pip install "pycsamt[torch,tensorflow,geo]"
```

### Development install

```bash
git clone https://github.com/earthai-tech/pycsamt.git
cd pycsamt
pip install -e ".[dev]"
```

---

## AI-assisted inversion

pyCSAMT v2 ships a dual-backend deep learning module that can be switched between
PyTorch and TensorFlow at runtime:

```python
from pycsamt.backends import set_backend
set_backend("torch")          # or "tensorflow"

from pycsamt.ai.inversion import Inv2DNet
net = Inv2DNet(n_stations=30, n_depth=60, n_freq=54)
net.fit(X_train, y_train, epochs=200)
rho_pred = net.predict(X_test)
net.save("inv2d_checkpoint.npz")
```

---

## Geological interpretation

The `pycsamt.interp` package converts raw inversion resistivity grids into
actionable lithological products:

```python
from pycsamt.interp import ResistivityModel, ModelCalibrator, RockDatabase
from pycsamt.interp.plot import PlotStratigraphicLog

rm  = ResistivityModel.from_occam2d(result)
cal = ModelCalibrator(ptol=0.25).fit(rm)

logs = cal.stratigraphic_logs(db=RockDatabase.default())
PlotStratigraphicLog(logs[0]).plot()
```

Supported export formats: **Oasis Montaj XYZ**, **CSV**, **LAS 2.0**, **VTK**.

---

## Documentation

Full documentation is hosted at **<https://pycsamt.readthedocs.io>** (under construction for v2).

---

## Citation

If you use pyCSAMT in your research, please cite:

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
  doi     = {10.1016/j.jappgeo.2022.104647},
}

@article{Kouadio2023,
  author  = {Kouadio, K. Laurent and others},
  title   = {Recovering the electrical resistivity from MT data using a
             deep resistivity convolutional neural network},
  journal = {Journal of Geophysical Research: Solid Earth},
  year    = {2023},
  doi     = {10.1029/2023JB027538},
}
```

---

## Contributing

Contributions are welcome — see [CONTRIBUTING](docs/source/contributing.rst) for guidelines.
Bug reports and feature requests: <https://github.com/earthai-tech/pycsamt/issues>

---

## License

pyCSAMT is distributed under the **GNU Lesser General Public License v3.0 or later**.
See [LICENSE.md](LICENSE.md) for the full text.

---

<p align="center">
  Developed by <a href="https://github.com/earthai-tech">earthai-tech</a> &nbsp;|&nbsp;
  Lead: <a href="mailto:etanoyau@gmail.com">Kouadio K. Laurent</a>
</p>

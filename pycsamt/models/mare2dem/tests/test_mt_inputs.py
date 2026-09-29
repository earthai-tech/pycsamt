"""MARE2DEM MT inversion inputs from EDI data (build_mt_inputs) and the
file-convention fixes found by running the real solver (2026-09-25).

The real run (CSAMT line, gfortran/OpenMPI/MKL build in WSL) showed:

* bounds were written as log10 (``-1, 5``) where MARE2DEM expects linear
  ohm-m (``0.1, 100000``, as in the vendored demo) -- the bandpass
  transform stopped bounding, trial models reached log10 rho ~ 74 and
  every worker segfaulted;
* ``Mare2DEMRunner`` stripped the iteration from ``mare2dem.0`` (it
  treated ``.0`` as an extension), and MARE2DEM refused to start;
* ``InversionResult`` loaded the first ``.resistivity``/``.resp`` found --
  the starting model -- instead of the last iteration.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.models.mare2dem import Mare2DEMConfig, build_mt_inputs, mt_core_grid
from pycsamt.models.mare2dem.iotools.resistivity import (
    ResistivityFile,
    read_resistivity,
)

ROOT = Path(__file__).resolve().parents[4]  # repository root
CSAMT = ROOT / Path("data/CSAMT")


def _sites():
    if not CSAMT.is_dir():
        pytest.skip("bundled CSAMT EDIs not found")
    from pycsamt.io import read_transfer_function
    from pycsamt.site.base import Sites

    return Sites([read_transfer_function(p)
                  for p in sorted(CSAMT.glob("*.edi"))])


def test_core_grid_spans_receivers_and_skin_depths():
    y, z = mt_core_grid([0.0, 100.0, 200.0, 300.0], [1000.0, 1.0],
                        initial_rho=100.0)
    assert y[0] < 0 and y[-1] > 300 and np.allclose(np.diff(y), 50.0)
    delta_hi = 503.0 * np.sqrt(100.0 / 1.0)
    assert z[0] == 0 and z[-1] >= 1.5 * delta_hi
    assert np.all(np.diff(np.diff(z)) >= -1e-9)  # thicknesses grow
    assert np.diff(z)[0] == pytest.approx(503.0 * np.sqrt(0.1) / 4)


def test_core_grid_rejects_empty():
    with pytest.raises(ValueError):
        mt_core_grid([], [1.0])


def test_build_mt_inputs_writes_a_runnable_set(tmp_path):
    out = build_mt_inputs(_sites(), tmp_path,
                          Mare2DEMConfig(initial_rho=50.0, max_iterations=4))
    files = out["files"]
    for key in ("data", "poly", "model", "settings"):
        assert files[key].is_file(), key
    # The runner argument keeps the iteration number (mare2dem.0)
    assert out["run_stem"] == "mare2dem.0"
    rf = read_resistivity(files["model"])
    # linear bounds, like MARE2DEM's own demo files
    assert list(rf.global_bounds) == [0.1, 100000.0]
    assert np.all(rf.bounds == 0)  # per-region: defer to global
    n_core = (len(out["y_edges"]) - 1) * (len(out["z_edges"]) - 1)
    assert out["n_parameters"] == n_core
    assert rf.num_regions == n_core + 2  # + fixed air and padding
    rho = np.asarray(rf.resistivity, float).ravel()
    assert np.allclose(rho[:n_core], 50.0)  # linear ohm-m, not log10
    assert out["padding"] >= 50000.0
    header = files["model"].read_text().splitlines()[:8]
    assert any("mare2dem.poly" in h for h in header)
    assert any("mare2dem.emdata" in h for h in header)


def test_write_resistivity_default_bounds_are_linear(tmp_path):
    from pycsamt.models.mare2dem.iotools.resistivity import write_resistivity

    rf = ResistivityFile(resistivity_file=str(tmp_path / "m.0.resistivity"))
    rf.resistivity = np.array([[10.0]])
    rf.free_parameter = np.array([[1]])
    rf.bounds = np.zeros((1, 2))
    rf.prejudice = np.zeros((1, 2))
    text = write_resistivity(rf, tmp_path / "m.0.resistivity").read_text()
    assert "0.1, 100000" in text


def test_input_builder_stub_writes_linear_rho(tmp_path):
    from pycsamt.models.mare2dem import InputBuilder

    p = InputBuilder(config=Mare2DEMConfig(initial_rho=100.0)) \
        .write_resistivity(tmp_path, filename="x.resistivity")
    rf = read_resistivity(p)
    assert float(np.asarray(rf.resistivity).ravel()[0]) == pytest.approx(100.0)


@pytest.mark.parametrize("arg, stem", [
    ("mare2dem.0", "mare2dem.0"),
    ("mare2dem.0.resistivity", "mare2dem.0"),
    ("run/demo.12.resistivity", "demo.12"),
    ("mare2dem", "mare2dem"),
])
def test_runner_keeps_iteration_in_stem(arg, stem):
    from pycsamt.models.mare2dem.runner import _run_stem

    assert _run_stem(arg) == stem


def test_result_loads_last_iteration(tmp_path):
    demo = ROOT / Path("data/mare2dem/demo_mt_inversion")
    if not demo.is_dir():
        pytest.skip("bundled MARE2DEM demo not found")
    from pycsamt.models.mare2dem import InversionResult

    res = InversionResult(demo)
    assert res.model_files[-1].name == "demo.6.resistivity"
    assert res.response_files[-1].name == "demo.6.resp"


def test_iteration_sort_key():
    from pycsamt.models.mare2dem.results import _iteration_of

    names = ["a.10.resistivity", "a.2.resistivity", "a.resistivity"]
    ordered = sorted((Path(n) for n in names), key=_iteration_of)
    assert [p.name for p in ordered] == ["a.resistivity", "a.2.resistivity",
                                         "a.10.resistivity"]


def test_plot_colour_limits_use_log10_of_linear_bounds(tmp_path):
    demo = ROOT / Path("data/mare2dem/demo_mt_inversion")
    if not demo.is_dir():
        pytest.skip("bundled MARE2DEM demo not found")
    import matplotlib

    matplotlib.use("Agg")
    from pycsamt.models.mare2dem import InversionResult
    from pycsamt.models.mare2dem.plot import PlotModel

    fig = PlotModel(InversionResult(demo)).plot()
    coll = [c for ax in fig.axes for c in ax.collections]
    assert coll
    lo, hi = coll[0].get_clim()
    assert lo == pytest.approx(-1.0) and hi == pytest.approx(5.0)

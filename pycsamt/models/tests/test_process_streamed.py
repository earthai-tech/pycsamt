"""Tests for pycsamt.models._process.run_streamed (live solver console).

Real child processes (the running Python interpreter) stand in for the
Fortran solvers.
"""

from __future__ import annotations

import subprocess
import sys
import time

import pytest

from pycsamt.models._process import ProcessCancelled, run_streamed, solver_env

PY = sys.executable


def test_lines_are_streamed_and_teed(tmp_path):
    got = []
    tee = tmp_path / "out.log"
    code = run_streamed(
        [PY, "-u", "-c", "import sys\nprint('iter 1 rms 3.2')\n"
                   "print('boom', file=sys.stderr)\nprint('iter 2')"],
        cwd=tmp_path, on_output=got.append, tee=tee)
    assert code == 0
    assert got[0] == "iter 1 rms 3.2" and "boom" in got and got[-1] == "iter 2"
    assert tee.read_text().splitlines() == got


def test_exit_code_is_returned(tmp_path):
    assert run_streamed([PY, "-c", "raise SystemExit(3)"], cwd=tmp_path) == 3


def test_cancel_stops_a_silent_process_quickly(tmp_path):
    """Cancel must not wait for output: a solver can be silent for
    minutes between iterations."""
    t0 = time.monotonic()
    calls = {"n": 0}

    def cancel():
        calls["n"] += 1
        return time.monotonic() - t0 > 0.5

    with pytest.raises(ProcessCancelled):
        run_streamed([PY, "-c", "import time; time.sleep(60)"],
                     cwd=tmp_path, cancel=cancel)
    assert time.monotonic() - t0 < 20
    assert calls["n"] >= 2


def test_timeout_kills(tmp_path):
    with pytest.raises(subprocess.TimeoutExpired):
        run_streamed([PY, "-c", "import time; time.sleep(60)"],
                     cwd=tmp_path, timeout=0.5)


def test_solver_env_unbuffers_gfortran():
    env = solver_env({"PATH": "x"})
    assert env["GFORTRAN_UNBUFFERED_ALL"] == "y" and env["PATH"] == "x"
    assert solver_env({"GFORTRAN_UNBUFFERED_ALL": "n"})[
        "GFORTRAN_UNBUFFERED_ALL"] == "n"


def test_modem_runner_streams_through_hook(tmp_path, monkeypatch):
    import pycsamt.models.modem.runner as rmod

    seen = {}

    def fake(cmd, **kw):
        seen.update(kw, cmd=cmd)
        kw["on_output"]("NLCG iteration 1")
        return 0

    monkeypatch.setattr(rmod, "run_streamed", fake)
    monkeypatch.setattr(rmod, "_resolve_binary", lambda name, wd: name)
    lines = []
    rmod.ModEmRunner(tmp_path).run(
        "m.rho", "d.dat", "c.inv", mode="2d", load_result=False,
        on_output=lines.append)
    assert lines == ["NLCG iteration 1"]
    assert seen["cancel"] is None and "-I" in seen["cmd"]


def test_mare2dem_runner_streams_through_hook(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem.runner as rmod
    from pycsamt.models.mare2dem.config import Mare2DEMConfig

    seen = {}

    def fake(cmd, **kw):
        seen.update(kw, cmd=cmd)
        return 0

    monkeypatch.setattr(rmod, "run_streamed", fake)
    cfg = Mare2DEMConfig(binary="wsl:/opt/MARE2DEM", n_procs=2)
    stop = lambda: False  # noqa: E731
    rmod.Mare2DEMRunner(tmp_path, config=cfg).run(
        "mare2dem", load_result=False, cancel=stop)
    assert seen["cmd"][:2] == ["wsl", "-e"] and seen["cancel"] is stop
    assert "GFORTRAN_UNBUFFERED_ALL=y" in seen["cmd"][-1]


def test_mare2dem_nonzero_exit_raises(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem.runner as rmod
    from pycsamt.models.mare2dem.config import Mare2DEMConfig

    monkeypatch.setattr(rmod, "run_streamed", lambda cmd, **kw: 2)
    cfg = Mare2DEMConfig(binary="wsl:/opt/MARE2DEM")
    with pytest.raises(subprocess.CalledProcessError):
        rmod.Mare2DEMRunner(tmp_path, config=cfg).run(
            "mare2dem", load_result=False, on_output=print)


def test_modem_args_are_relative_inside_workdir(tmp_path, monkeypatch):
    """Mod2DMT/Mod3DMT read each argument into an 80-character buffer; a
    long absolute run-folder path was cut and the solver stopped with
    "Please specify a valid inverse control file" (found running the real
    Mod2DMT from the desktop Inversion Studio, 2026-09-25)."""
    import pycsamt.models.modem.runner as rmod

    wd = tmp_path / ("x" * 90)
    wd.mkdir()
    seen = {}
    monkeypatch.setattr(rmod, "run_streamed",
                        lambda cmd, **kw: seen.setdefault("cmd", cmd) and 0)
    monkeypatch.setattr(rmod, "_resolve_binary", lambda name, w: name)
    runner = rmod.ModEmRunner(wd)
    runner.run(wd / "m.rho", wd / "d.dat", wd / "c.inv", mode="2d",
               load_result=False, on_output=print)
    assert seen["cmd"][-3:] == ["m.rho", "d.dat", "c.inv"]
    outside = tmp_path / "elsewhere.rho"
    assert runner._arg(outside) == str(outside)  # left untouched


# ── ModEM 2-D inputs from scalar CSAMT data (found with the real Mod2DMT) ──


def _csamt_sites():
    from pathlib import Path

    root = Path(__file__).resolve().parents[3] / "data" / "CSAMT"
    if not root.is_dir():
        pytest.skip("bundled CSAMT EDIs not found")
    from pycsamt.io import read_transfer_function
    from pycsamt.site.base import Sites

    return Sites([read_transfer_function(p) for p in sorted(root.glob("*.edi"))])


def test_modem_errors_finite_for_edis_without_errors(tmp_path):
    """Scalar CSAMT EDIs have Zxy only and no error blocks; a NaN-blind
    max() wrote every error as NAN and Mod2DMT's objective was NaN."""
    import numpy as np

    from pycsamt.models.modem import InputBuilder, ModEmConfig, ModEmData

    cfg = ModEmConfig(mode="2d", component_type="TE_Impedance")
    files = InputBuilder(config=cfg).build(_csamt_sites(), tmp_path)
    assert "NAN" not in files["data"].read_text().upper()
    data = ModEmData.read(files["data"])
    errs = np.array([r[-1] for b in data.blocks for r in b["rows"]]
                    if data.blocks and "rows" in data.blocks[0] else [1.0])
    assert np.isfinite(errs).all() and (errs > 0).all()


def test_modem_2d_mesh_has_no_zero_width_cells(tmp_path):
    """csa350 and csa400 share a position; the zero gap between them became
    a zero-width column."""
    from pycsamt.models.modem import InputBuilder, ModEmConfig, ModEmModel2D

    cfg = ModEmConfig(mode="2d", component_type="TE_Impedance")
    files = InputBuilder(config=cfg).build(_csamt_sites(), tmp_path)
    model = ModEmModel2D.read(files["model"])
    assert (model.x_widths > 0).all()

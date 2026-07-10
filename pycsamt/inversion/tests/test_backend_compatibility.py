# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Compatibility tests for inversion backends.

The optional SimPEG/pyGIMLi smoke tests are intentionally gated by an
environment variable so normal CI does not require heavy external engines.
Run them with ``PYCSAMT_RUN_OPTIONAL_BACKENDS=1`` in an environment where the
backend packages are installed.
"""

from __future__ import annotations

import importlib.util
import os

import numpy as np
import pytest

from pycsamt.inversion import (
    InversionConfig,
    InversionMesh,
    InversionResult,
    InversionWorkflow,
    StartingModel,
    available_backends,
    get_backend,
)

_RUN_OPTIONAL = os.environ.get("PYCSAMT_RUN_OPTIONAL_BACKENDS") == "1"


def test_backend_registry_instantiates_all_adapters():
    expected = {"builtin", "occam2d", "modem", "simpeg", "pygimli"}
    assert expected.issubset(set(available_backends()))
    for name in expected:
        backend_cls = get_backend(name)
        cfg = InversionConfig(backend=name, dimension="3d" if name == "modem" else "1d")
        if name == "occam2d":
            cfg = InversionConfig(backend=name, dimension="2d")
        assert backend_cls(cfg).name == name


def test_backend_result_conversion_contracts_for_1d_models():
    model = StartingModel([100.0, 250.0, 500.0], [50.0, 150.0])
    for backend in ("builtin", "simpeg", "pygimli"):
        result = InversionResult(
            method="mt",
            dimension="1d",
            backend=backend,
            model=model,
            mesh=InversionMesh.for_1d([25.0, 125.0, 275.0]),
            rms=1.0,
        )
        converted = result.to_resistivity_model()
        assert converted.rho_2d.shape == (3, 1)
        assert converted.method == f"{backend}:mt:1d"


def test_backend_result_conversion_contracts_for_2d_sections():
    rho_2d = np.log10(np.array([[100.0, 120.0], [50.0, 80.0]]))
    model = {
        "rho_2d": rho_2d,
        "x_centers": [0.0, 500.0],
        "z_centers": [25.0, 125.0],
        "station_x": [0.0, 500.0],
        "station_names": ["S1", "S2"],
    }
    for backend in ("builtin", "simpeg", "pygimli", "occam2d"):
        result = InversionResult(
            method="mt",
            dimension="2d",
            backend=backend,
            model=model,
            mesh=InversionMesh(dimension="2d", x_centers=[0.0, 500.0], z_centers=[25.0, 125.0]),
            rms=1.0,
        )
        converted = result.to_resistivity_model()
        assert converted.rho_2d.shape == (2, 2)
        assert converted.station_names == ["S1", "S2"]
        assert converted.method == f"{backend}:mt"


def test_external_ready_results_fail_conversion_without_model(tmp_path):
    for backend, dimension in (("occam2d", "2d"), ("modem", "3d")):
        result = InversionResult(
            method="mt",
            dimension=dimension,
            backend=backend,
            status="ready",
            files={"workdir": str(tmp_path)},
        )
        with pytest.raises(ValueError, match="no model"):
            result.to_resistivity_model()


@pytest.mark.skipif(
    not _RUN_OPTIONAL,
    reason="set PYCSAMT_RUN_OPTIONAL_BACKENDS=1 to run installed SimPEG smoke tests",
)
@pytest.mark.skipif(
    importlib.util.find_spec("simpeg") is None,
    reason="SimPEG is not installed",
)
def test_installed_simpeg_backend_smoke(tmp_path):
    data = {
        "freqs": [1.0, 10.0],
        "rho_a": [100.0, 120.0],
        "phase": [45.0, 46.0],
    }
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data=data,
        n_layers=2,
        max_iter=1,
        starting_model=StartingModel([100.0, 200.0], [100.0]),
        backend_options={"n_cells": 4, "depth_max": 500.0, "estimate_beta": False},
        workdir=str(tmp_path / "simpeg"),
    )
    result = InversionWorkflow(cfg).run()
    assert result.backend == "simpeg"
    assert result.dimension == "1d"
    assert np.isfinite(result.rms)


@pytest.mark.skipif(
    not _RUN_OPTIONAL,
    reason="set PYCSAMT_RUN_OPTIONAL_BACKENDS=1 to run installed pyGIMLi smoke tests",
)
@pytest.mark.skipif(
    importlib.util.find_spec("pygimli") is None,
    reason="pyGIMLi is not installed",
)
def test_installed_pygimli_backend_smoke(tmp_path):
    data = {
        "freqs": [1.0, 10.0],
        "rho_a": [100.0, 120.0],
        "phase": [45.0, 46.0],
    }
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="pygimli",
        data=data,
        n_layers=2,
        max_iter=1,
        starting_model=StartingModel([100.0, 200.0], [100.0]),
        backend_options={"lam": 20.0},
        workdir=str(tmp_path / "pygimli"),
    )
    result = InversionWorkflow(cfg).run()
    assert result.backend == "pygimli"
    assert result.dimension == "1d"
    assert np.isfinite(result.rms)

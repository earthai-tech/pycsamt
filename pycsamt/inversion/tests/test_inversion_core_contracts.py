from __future__ import annotations

import json
import zipfile
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.inversion.base import BaseInversionBackend
from pycsamt.inversion.config import InversionConfig
from pycsamt.inversion.export import (
    to_archive,
    to_csv,
    to_geojson,
    to_npz,
    to_vtk,
)
from pycsamt.inversion.mesh import (
    InversionMesh,
    build_1d_tensor_mesh,
    core_rho_from_start,
    depth_widths,
)
from pycsamt.inversion.model import StartingModel
from pycsamt.inversion.objective import (
    ErrorModel,
    component_errors,
    weighted_rms,
)
from pycsamt.inversion.plot import plot_model, plot_rms
from pycsamt.inversion.regularization import (
    Regularization,
    regularization_residual,
)
from pycsamt.inversion.results import (
    InversionHistory,
    InversionResult,
    InversionUncertainty,
)
from pycsamt.inversion.workflow import InversionWorkflow


def _section_result() -> InversionResult:
    rho = np.log10(np.array([
        [100.0, 150.0, 200.0],
        [80.0, 120.0, 180.0],
    ]))
    return InversionResult(
        method="mt",
        dimension="2d",
        backend="builtin",
        model={
            "rho_2d": rho,
            "x_centers": [0.0, 500.0, 1000.0],
            "z_centers": [50.0, 150.0],
            "station_x": [0.0, 500.0, 1000.0],
            "station_names": ["S1", "S2", "S3"],
        },
        rms=1.2,
        objective=4.5,
        uncertainty={
            "confidence": np.full_like(rho, 0.8),
            "model_std": np.full_like(rho, 0.1),
        },
        history=[
            {"iteration": 0, "objective": 9.0, "rms": 2.1},
            {"iteration": 1, "objective": 4.5, "rms": 1.2},
        ],
        metadata={"station_rms": [1.1, 1.2, 1.3]},
    )


def test_config_validation_and_template_roundtrip(tmp_path):
    cfg = InversionConfig(
        method="AMT",
        dimension="2D",
        backend="builtin",
        output_dir=tmp_path / "out",
        backend_options={"alpha_x": 2.0},
    )
    assert cfg.method == "amt"
    assert cfg.dimension == "2d"
    assert cfg.output_path == tmp_path / "out"

    path = cfg.write_template(tmp_path / "inversion.json", fmt="json")
    loaded = InversionConfig.from_file(path)
    assert loaded.method == "amt"
    assert loaded.dimension == "2d"
    assert loaded.backend_options["alpha_x"] == 2.0

    with pytest.raises(ValueError):
        InversionConfig(method="gravity")
    with pytest.raises(ValueError):
        InversionConfig(n_layers=1)
    with pytest.raises(ValueError):
        InversionConfig(phase_error=0.0)


def test_starting_model_coercion_validation_and_layered_adapter():
    obj = SimpleNamespace(
        resistivity=[50.0, 100.0, 300.0],
        thickness=[25.0, 75.0],
        name="object_model",
    )
    model = StartingModel.coerce(obj)
    assert model.name == "object_model"
    np.testing.assert_allclose(model.depths, [0.0, 25.0, 100.0])
    layered = model.to_layered_model()
    assert layered.n_layers == 3

    same = StartingModel.coerce(model)
    assert same is model
    with pytest.raises(ValueError):
        StartingModel([100.0], [])
    with pytest.raises(ValueError):
        StartingModel([100.0, 200.0], [0.0])


def test_mesh_builders_and_core_rho_padding():
    widths = depth_widths(120.0, 3, {"growth_factor": 1.0})
    np.testing.assert_allclose(widths, [40.0, 40.0, 40.0])

    class TensorMesh:
        def __init__(self, widths, origin="0"):
            self.widths = widths
            self.origin = origin

    start = StartingModel([100.0, 300.0], [50.0])
    mesh, z = build_1d_tensor_mesh(
        start,
        {"n_cells": 3, "depth_max": 120.0, "growth_factor": 1.0},
        TensorMesh,
    )
    assert mesh.origin == "0"
    np.testing.assert_allclose(z, [20.0, 60.0, 100.0])

    rho = core_rho_from_start(start, nz=2, nx=3)
    assert rho.shape == (2, 3)
    np.testing.assert_allclose(rho[:, 0], [100.0, 300.0])

    with pytest.raises(ValueError):
        InversionMesh(dimension="4d")
    with pytest.raises(ValueError):
        depth_widths(0.0, 3)


def test_objective_errors_masks_and_weighted_rms_edges():
    model = ErrorModel(
        rho_relative=0.1,
        rho_absolute=2.0,
        phase_absolute=4.0,
        masks={"station": [True, False]},
    )
    rho = np.array([[10.0, 100.0], [20.0, 200.0]])
    err = model.errors(rho, component="rho")
    np.testing.assert_allclose(err, [[2.0, 10.0], [2.0, 20.0]])
    np.testing.assert_array_equal(
        model.mask(rho, component="rho"),
        [[True, True], [False, False]],
    )
    np.testing.assert_allclose(
        component_errors(rho, SimpleNamespace(error_floor=0.2), component="rho"),
        np.maximum(np.abs(rho) * 0.2, 1e-12),
    )
    assert weighted_rms(
        [1.0, np.nan, 3.0],
        [1.0, 2.0, 1.0],
        [1.0, 1.0, 2.0],
    ) == pytest.approx(np.sqrt(0.5))
    assert np.isnan(weighted_rms([np.nan], [0.0], [1.0]))


def test_regularization_residual_shapes_and_kinds():
    values = np.array([[1.0, 2.0, 4.0], [1.5, 2.5, 5.0]])
    smooth = regularization_residual(
        values,
        regularization=Regularization(kind="smooth", alpha_x=4.0, alpha_z=9.0),
    )
    assert smooth.size == 7
    damped = regularization_residual(
        values,
        reference=np.ones_like(values),
        regularization=Regularization(kind="damped", alpha_s=1.0),
    )
    assert damped.size > smooth.size
    blocky = regularization_residual(
        values,
        regularization=Regularization(kind="blocky", alpha_x=1.0, alpha_z=1.0),
    )
    assert blocky.size == smooth.size
    none = regularization_residual(values, regularization=Regularization(kind="none"))
    assert none.size == 0


def test_result_containers_convert_dicts_and_arrays():
    result = _section_result()
    assert result.converged is True
    assert isinstance(result.uncertainty, InversionUncertainty)
    assert isinstance(result.history, InversionHistory)
    arrays = result.history.arrays()
    np.testing.assert_allclose(arrays["objective"], [9.0, 4.5])

    model = result.to_resistivity_model()
    assert model.rho_2d.shape == (2, 3)
    assert model.station_names == ["S1", "S2", "S3"]

    layered = InversionResult(
        method="mt",
        dimension="1d",
        backend="builtin",
        model=StartingModel([100.0, 250.0, 500.0], [20.0, 80.0]),
    )
    converted = layered.to_resistivity_model()
    assert converted.rho_2d.shape == (3, 1)
    np.testing.assert_allclose(converted.rho_2d[:, 0], np.log10([100.0, 250.0, 500.0]))

    with pytest.raises(ValueError):
        InversionResult("mt", "1d", "builtin").to_resistivity_model()


def test_exports_write_expected_payloads(tmp_path):
    result = _section_result()
    native = tmp_path / "native.dat"
    native.write_text("native\n", encoding="utf-8")
    result.files["mesh"] = str(native)

    csv_path = to_csv(result, tmp_path / "model.csv", log_rho=False)
    assert csv_path.read_text(encoding="utf-8").splitlines()[0].startswith("x_m,z_m,rho_ohm_m")

    npz_path = to_npz(result, tmp_path / "model.npz")
    with np.load(npz_path) as data:
        assert data["rho_2d"].shape == (2, 3)
        assert "uncertainty_confidence" in data
        assert "history_objective" in data

    geojson_path = to_geojson(result, tmp_path / "model.geojson")
    geo = json.loads(geojson_path.read_text(encoding="utf-8"))
    assert geo["type"] == "FeatureCollection"
    assert len(geo["features"]) == 6
    assert "uncertainty_confidence" in geo["features"][0]["properties"]

    vtk_path = to_vtk(result, tmp_path / "model.vtk")
    vtk_text = vtk_path.read_text(encoding="utf-8")
    assert "RECTILINEAR_GRID" in vtk_text
    assert "rho_log10_ohm_m" in vtk_text

    archive_path = to_archive(result, tmp_path / "snapshot.zip")
    with zipfile.ZipFile(archive_path) as zf:
        names = set(zf.namelist())
    assert {"metadata.json", "result.npz", "model.csv"}.issubset(names)
    assert any(name.startswith("native_files/mesh_") for name in names)


def test_plot_helpers_use_existing_axes_and_save(tmp_path):
    result = _section_result()
    fig, (ax_model, ax_rms) = plt.subplots(1, 2, figsize=(8, 3))
    out_model = tmp_path / "model.png"
    out_rms = tmp_path / "rms.png"

    returned_model = plot_model(
        result,
        ax=ax_model,
        colorbar=False,
        savepath=str(out_model),
        savefig_kw={"dpi": 80},
    )
    returned_rms = plot_rms(
        result,
        ax=ax_rms,
        savepath=str(out_rms),
        savefig_kw={"dpi": 80},
    )

    assert returned_model is ax_model
    assert returned_rms is ax_rms
    assert out_model.exists()
    assert out_rms.exists()
    assert ax_rms.get_ylabel() == "Weighted RMS"
    plt.close(fig)


def test_backend_support_checks_and_workflow_merges_config(monkeypatch):
    class DummyBackend(BaseInversionBackend):
        name = "dummy"
        supports = (("mt", "1d"),)

        def run(self, data=None):
            self.check_supported()
            return InversionResult(
                self.config.method,
                self.config.dimension,
                self.name,
                status="prepared",
            )

    cfg = InversionConfig(method="mt", dimension="1d", backend="builtin")
    backend = DummyBackend(cfg)
    assert backend.run().status == "prepared"

    backend.config = InversionConfig(method="tdem", dimension="1d", backend="builtin")
    with pytest.raises(NotImplementedError):
        backend.check_supported()

    monkeypatch.setattr(
        InversionConfig,
        "to_backend",
        lambda self: DummyBackend(self),
    )
    wf = InversionWorkflow({"method": "mt"}, dimension="1d", backend="builtin")
    assert wf.config.method == "mt"
    assert wf.config.dimension == "1d"
    assert wf.run().backend == "dummy"

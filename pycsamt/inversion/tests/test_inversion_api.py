# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from pycsamt.forward import (
    Grid2D,
    LayeredModel,
    MT1DForward,
    MT2DForward,
    TEM1DForward,
)
from pycsamt.inversion import (
    EMData,
    ErrorModel,
    InversionConfig,
    InversionHistory,
    InversionResult,
    InversionUncertainty,
    InversionWorkflow,
    Regularization,
    StartingModel,
    available_backends,
    component_errors,
    component_mask,
    error_model_from_config,
    export,
    plot,
    pygimli_lambda,
    regularization_from_config,
    regularization_residual,
    regularization_weight,
    run_inversion,
)
from pycsamt.inversion.amt import (
    AMT1DInversion,
    AMT2DInversion,
)
from pycsamt.inversion.backends import modem as modem_backend
from pycsamt.inversion.backends import (
    occam2d as occam2d_backend,
)
from pycsamt.inversion.backends import (
    pygimli as pygimli_backend,
)
from pycsamt.inversion.backends import (
    simpeg as simpeg_backend,
)
from pycsamt.inversion.backends.builtin import (
    Builtin1DBackend,
)
from pycsamt.inversion.backends.modem import ModEMBackend
from pycsamt.inversion.backends.occam2d import Occam2DBackend
from pycsamt.inversion.backends.pygimli import PyGIMLiBackend
from pycsamt.inversion.backends.simpeg import SimPEGBackend
from pycsamt.inversion.emap import EMAP2DInversion
from pycsamt.inversion.mesh import (
    InversionMesh,
    build_1d_tensor_mesh,
    build_3d_tensor_mesh,
    build_fd2d_grid,
    core_rho_from_start,
    depth_widths,
)
from pycsamt.inversion.mt import (
    MT1DInversion,
    MT2DInversion,
    MT3DInversion,
)
from pycsamt.inversion.objective import (
    relative_errors,
    weighted_rms,
)
from pycsamt.inversion.tdem import (
    TDEM1DInversion,
    TDEM2DInversion,
)


def test_config_template_roundtrip(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="builtin",
        n_layers=3,
        workdir=str(tmp_path / "run"),
    )
    path = cfg.write_template(tmp_path / "inv_config.py")
    loaded = InversionConfig.from_file(path)
    assert loaded.method == "mt"
    assert loaded.dimension == "1d"
    assert loaded.n_layers == 3


def test_emdata_coerce_aliases():
    data = EMData.coerce(
        {
            "method": "AMT",
            "freqs": [1.0, 10.0],
            "rho_a": [100.0, 120.0],
            "phase": [45.0, 47.0],
            "stations": ["S1"],
        }
    )
    assert data.method == "amt"
    assert data.has_mt_response
    assert data.n_samples == 2
    assert data.station_names == ["S1"]


def test_emdata_public_docstrings_describe_coercion_contracts():
    assert "Shape convention" in EMData.__doc__
    assert "response-like objects" in EMData.__doc__
    assert "from pycsamt.inversion.data import EMData" in EMData.__doc__
    assert "References" in EMData.__doc__
    assert "periods" in EMData.from_dict.__doc__
    assert (
        "from pycsamt.inversion.data import EMData"
        in EMData.from_dict.__doc__
    )
    assert "to_soundings()" in EMData.coerce.__doc__
    assert "class Response" in EMData.coerce.__doc__


def test_export_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.export import to_csv" in export.to_csv.__doc__
    )
    assert "RFC 4180" in export.to_csv.__doc__
    assert (
        "from pycsamt.inversion.export import to_npz" in export.to_npz.__doc__
    )
    assert "numpy.savez_compressed" in export.to_npz.__doc__
    assert (
        "from pycsamt.inversion.export import to_geojson"
        in export.to_geojson.__doc__
    )
    assert "RFC 7946" in export.to_geojson.__doc__
    assert (
        "from pycsamt.inversion.export import to_vtk" in export.to_vtk.__doc__
    )
    assert "Visualization Toolkit" in export.to_vtk.__doc__
    assert (
        "from pycsamt.inversion.export import to_geotiff"
        in export.to_geotiff.__doc__
    )
    assert "GeoTIFF" in export.to_geotiff.__doc__
    assert (
        "from pycsamt.inversion.export import to_archive"
        in export.to_archive.__doc__
    )
    assert "metadata.json" in export.to_archive.__doc__


def test_mesh_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.mesh import InversionMesh"
        in InversionMesh.__doc__
    )
    assert "Oldenburg" in InversionMesh.__doc__
    assert (
        "from pycsamt.inversion.mesh import depth_widths"
        in depth_widths.__doc__
    )
    assert "Ward" in depth_widths.__doc__
    assert (
        "from pycsamt.inversion.mesh import build_1d_tensor_mesh"
        in build_1d_tensor_mesh.__doc__
    )
    assert (
        "from pycsamt.inversion.mesh import build_3d_tensor_mesh"
        in build_3d_tensor_mesh.__doc__
    )
    assert (
        "from pycsamt.inversion.mesh import build_fd2d_grid"
        in build_fd2d_grid.__doc__
    )
    assert (
        "from pycsamt.inversion.mesh import core_rho_from_start"
        in core_rho_from_start.__doc__
    )


def test_model_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.model import StartingModel"
        in StartingModel.__doc__
    )
    assert "ReferenceModel" in StartingModel.__doc__
    assert "Constable" in StartingModel.__doc__
    assert (
        "from pycsamt.inversion.model import StartingModel"
        in StartingModel.default.__doc__
    )
    assert (
        "from pycsamt.inversion.model import StartingModel"
        in StartingModel.from_dict.__doc__
    )
    assert (
        "from pycsamt.inversion.model import StartingModel"
        in StartingModel.coerce.__doc__
    )
    assert (
        "from pycsamt.inversion.model import StartingModel"
        in StartingModel.to_layered_model.__doc__
    )


def test_objective_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.objective import ErrorModel"
        in ErrorModel.__doc__
    )
    assert "Tarantola" in ErrorModel.__doc__
    assert (
        "from pycsamt.inversion.objective import ErrorModel"
        in ErrorModel.errors.__doc__
    )
    assert (
        "from pycsamt.inversion.objective import ErrorModel"
        in ErrorModel.mask.__doc__
    )
    assert (
        "from pycsamt.inversion.objective import relative_errors"
        in relative_errors.__doc__
    )
    assert (
        "from pycsamt.inversion.objective import error_model_from_config"
        in error_model_from_config.__doc__
    )
    assert (
        "from pycsamt.inversion.objective import component_errors"
        in component_errors.__doc__
    )
    assert (
        "from pycsamt.inversion.objective import component_mask"
        in component_mask.__doc__
    )
    assert (
        "from pycsamt.inversion.objective import weighted_rms"
        in weighted_rms.__doc__
    )


def test_regularization_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.regularization import Regularization"
        in Regularization.__doc__
    )
    assert "Tikhonov" in Regularization.__doc__
    assert (
        "from pycsamt.inversion.regularization import Regularization"
        in Regularization.validate.__doc__
    )
    assert (
        "from pycsamt.inversion.regularization import regularization_from_config"
        in regularization_from_config.__doc__
    )
    assert (
        "from pycsamt.inversion.regularization import regularization_weight"
        in regularization_weight.__doc__
    )
    assert (
        "from pycsamt.inversion.regularization import pygimli_lambda"
        in pygimli_lambda.__doc__
    )
    assert (
        "from pycsamt.inversion.regularization import regularization_residual"
        in regularization_residual.__doc__
    )
    assert "Farquharson" in regularization_residual.__doc__


def test_results_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.results import InversionHistory"
        in InversionHistory.__doc__
    )
    assert "Constable" in InversionHistory.__doc__
    assert (
        "from pycsamt.inversion.results import InversionHistory"
        in InversionHistory.arrays.__doc__
    )
    assert (
        "from pycsamt.inversion.results import InversionUncertainty"
        in InversionUncertainty.__doc__
    )
    assert (
        "from pycsamt.inversion.results import InversionResult"
        in InversionResult.__doc__
    )
    assert "to_resistivity_model" in InversionResult.__doc__
    assert (
        "from pycsamt.inversion.results import InversionResult"
        in InversionResult.converged.__doc__
    )
    assert (
        "from pycsamt.inversion.results import InversionResult"
        in InversionResult.to_resistivity_model.__doc__
    )
    assert (
        "from pycsamt.inversion.results import InversionResult"
        in InversionResult.summary.__doc__
    )


def test_plot_public_docstrings_include_examples_and_references():
    assert "Available plots" in plot.__doc__
    assert "from pycsamt.inversion.plot import plot_model" in plot.__doc__
    assert (
        "from pycsamt.inversion.plot import plot_model"
        in plot.plot_model.__doc__
    )
    assert "to_resistivity_model" in plot.plot_model.__doc__
    assert "Tufte" in plot.plot_model.__doc__
    assert (
        "from pycsamt.inversion.plot import plot_rms" in plot.plot_rms.__doc__
    )
    assert "station_rms" in plot.plot_rms.__doc__
    assert "Aster" in plot.plot_rms.__doc__


def test_workflow_public_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.workflow import InversionWorkflow"
        in InversionWorkflow.__doc__
    )
    assert (
        "from pycsamt.inversion.workflow import run_inversion"
        in InversionWorkflow.__doc__
    )
    assert "Aster" in InversionWorkflow.__doc__
    assert (
        "from pycsamt.inversion.workflow import InversionWorkflow"
        in InversionWorkflow.__init__.__doc__
    )
    assert (
        "from pycsamt.inversion.workflow import InversionWorkflow"
        in InversionWorkflow.run.__doc__
    )
    assert (
        "from pycsamt.inversion.workflow import run_inversion"
        in run_inversion.__doc__
    )


def test_builtin_backend_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.backends.builtin import Builtin1DBackend"
        in Builtin1DBackend.__doc__
    )
    assert "finite-difference 2-D MT profile mode" in Builtin1DBackend.__doc__
    assert "Constable" in Builtin1DBackend.__doc__
    assert (
        "from pycsamt.inversion.backends.builtin import Builtin1DBackend"
        in Builtin1DBackend.run.__doc__
    )
    assert "scipy" in Builtin1DBackend.run.__doc__.lower()


def test_modem_backend_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.backends.modem import ModEMBackend"
        in ModEMBackend.__doc__
    )
    assert "ModEM ecosystem" in ModEMBackend.__doc__
    assert "Kelbert" in ModEMBackend.__doc__
    assert (
        "from pycsamt.inversion.backends.modem import ModEMBackend"
        in ModEMBackend.run.__doc__
    )
    assert "run_external" in ModEMBackend.run.__doc__


def test_occam2d_backend_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.backends.occam2d import Occam2DBackend"
        in Occam2DBackend.__doc__
    )
    assert "Occam2D-style workflows" in Occam2DBackend.__doc__
    assert "deGroot-Hedlin" in Occam2DBackend.__doc__
    assert (
        "from pycsamt.inversion.backends.occam2d import Occam2DBackend"
        in Occam2DBackend.run.__doc__
    )
    assert "run_external" in Occam2DBackend.run.__doc__


def test_pygimli_backend_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.backends.pygimli import PyGIMLiBackend"
        in PyGIMLiBackend.__doc__
    )
    assert "stitched station-by-station 1-D" in PyGIMLiBackend.__doc__
    assert "2-D finite-element pyGIMLi EM inversion" in PyGIMLiBackend.__doc__
    assert "mt_operator" in PyGIMLiBackend.__doc__
    assert "Rücker" in PyGIMLiBackend.__doc__
    assert (
        "from pycsamt.inversion.backends.pygimli import PyGIMLiBackend"
        in PyGIMLiBackend.run.__doc__
    )
    assert "pygimli" in PyGIMLiBackend.run.__doc__.lower()


def test_simpeg_backend_docstrings_include_examples_and_references():
    assert (
        "from pycsamt.inversion.backends.simpeg import SimPEGBackend"
        in SimPEGBackend.__doc__
    )
    assert "Simulation3DPrimarySecondary" in SimPEGBackend.__doc__
    assert "log conductivity" in SimPEGBackend.__doc__
    assert "target_chifact" in SimPEGBackend.__doc__
    assert "Cockett" in SimPEGBackend.__doc__
    assert (
        "from pycsamt.inversion.backends.simpeg import SimPEGBackend"
        in SimPEGBackend.run.__doc__
    )
    assert "rho_3d" in SimPEGBackend.run.__doc__


def test_emdata_coerce_forward_response_object():
    class Response:
        freqs = np.array([1.0, 10.0])
        rho_a = np.array([100.0, 120.0])
        phase = np.array([45.0, 47.0])
        station_name = "S01"

    data = EMData.coerce(Response(), method="mt")
    assert data.has_mt_response
    np.testing.assert_allclose(data.frequencies, [1.0, 10.0])
    np.testing.assert_allclose(data.rho_a, [100.0, 120.0])
    assert data.station_names == ["S01"]


def test_emdata_coerce_2d_response_transposes_frequency_station_arrays():
    class Response2D:
        freqs = np.array([1.0, 10.0])
        rho_a_te = np.array([[100.0, 110.0, 120.0], [200.0, 210.0, 220.0]])
        phase_te = np.array([[45.0, 46.0, 47.0], [50.0, 51.0, 52.0]])
        station_names = ["S1", "S2", "S3"]
        x_stations = [0.0, 500.0, 1000.0]

    data = EMData.coerce(Response2D(), method="mt")
    assert data.rho_a.shape == (3, 2)
    np.testing.assert_allclose(data.rho_a[1], [110.0, 210.0])
    assert data.station_names == ["S1", "S2", "S3"]


def test_emdata_coerce_station_response_collection():
    class Station:
        def __init__(self, name, x, rho):
            self.name = name
            self.x = x
            self.freqs = np.array([1.0, 10.0])
            self.rho_a = np.asarray(rho, dtype=float)
            self.phase = np.array([45.0, 47.0])

    data = EMData.coerce(
        [
            Station("S1", 0.0, [100.0, 120.0]),
            Station("S2", 500.0, [110.0, 130.0]),
        ]
    )
    assert data.rho_a.shape == (2, 2)
    assert data.station_names == ["S1", "S2"]
    np.testing.assert_allclose(data.station_x, [0.0, 500.0])


def test_emdata_coerce_tdem_soundings_collection():
    class Sounding:
        def __init__(self, name, x, values):
            self.station_name = name
            self.x = x
            self.time_gates = np.array([1e-4, 1e-3])
            self.data = np.asarray(values, dtype=float)
            self.error = np.array([1e-10, 1e-11])

    data = EMData.coerce(
        [
            Sounding("T1", 0.0, [1e-8, 1e-9]),
            Sounding("T2", 100.0, [2e-8, 2e-9]),
        ]
    )
    assert data.method == "tdem"
    assert data.has_tdem_response
    assert data.values.shape == (2, 2)
    assert data.station_names == ["T1", "T2"]


def test_emdata_coerce_tdem_survey_object():
    class Survey:
        def to_soundings(self):
            class Sounding:
                station_name = "T1"
                x = 10.0
                time_gates = np.array([1e-4, 1e-3])
                data = np.array([1e-8, 1e-9])

            return [Sounding()]

    data = EMData.coerce(Survey())
    assert data.method == "tdem"
    np.testing.assert_allclose(data.times, [1e-4, 1e-3])
    assert data.station_names == ["T1"]


def test_backend_registry():
    names = available_backends()
    assert "builtin" in names
    assert "occam2d" in names
    assert "simpeg" in names


def test_simpeg_backend_missing_dependency_is_clear():
    if importlib.util.find_spec("simpeg") is not None:
        pytest.skip(
            "SimPEG is installed; missing-dependency path not applicable."
        )
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0], "phase": [45.0]},
    )
    with pytest.raises(ImportError, match="SimPEG"):
        InversionWorkflow(cfg).run()


def test_simpeg_observation_packing_order_and_errors():
    data = EMData.coerce(
        {
            "freqs": [1.0, 10.0],
            "rho_a": [100.0, 200.0],
            "phase": [45.0, 50.0],
        }
    )
    cfg = InversionConfig(error_floor=0.05, phase_error=2.0)
    values, errors = simpeg_backend._pack_nsem_observations(data, cfg)
    np.testing.assert_allclose(values, [100.0, 45.0, 200.0, 50.0])
    np.testing.assert_allclose(errors, [5.0, 2.0, 10.0, 2.0])


def test_simpeg_3d_observation_packing_and_station_locations():
    data = EMData.coerce(
        {
            "freqs": [1.0, 10.0],
            "rho_a": [[100.0, 200.0], [110.0, 210.0]],
            "phase": [[45.0, 50.0], [46.0, 51.0]],
            "station_x": [0.0, 500.0],
            "metadata": {"station_y": [10.0, 20.0]},
        }
    )
    cfg = InversionConfig(error_floor=0.05, phase_error=2.0)
    values, errors = simpeg_backend._pack_nsem_observations(data, cfg)
    np.testing.assert_allclose(
        values,
        [100.0, 110.0, 45.0, 46.0, 200.0, 210.0, 50.0, 51.0],
    )
    np.testing.assert_allclose(
        errors, [5.0, 5.5, 2.0, 2.0, 10.0, 10.5, 2.0, 2.0]
    )
    locations = simpeg_backend._station_locations(data)
    np.testing.assert_allclose(
        locations, [[0.0, 10.0, 0.0], [500.0, 20.0, 0.0]]
    )


def test_simpeg_3d_model_reshape_helper():
    class Mesh:
        shape_cells = (2, 2, 2)

    log_sigma = np.log(np.arange(1.0, 9.0))
    rho = simpeg_backend._rho_3d_from_log_sigma(log_sigma, Mesh())
    assert rho.shape == (2, 2, 2)
    np.testing.assert_allclose(
        np.sort(rho.reshape(-1)),
        np.sort(1.0 / np.arange(1.0, 9.0)),
    )


def test_mesh_builders_centralize_tensor_mesh_setup():
    class TensorMesh:
        def __init__(self, widths, origin="0"):
            self.widths = [np.asarray(w, dtype=float) for w in widths]
            self.origin = origin
            self.nC = int(np.prod([w.size for w in self.widths]))
            self.shape_cells = tuple(w.size for w in self.widths)

    start = StartingModel([100.0, 200.0], [50.0])
    widths = depth_widths(
        100.0, 4, {"min_cell_size": 5.0, "growth_factor": 1.0}
    )
    np.testing.assert_allclose(widths, [25.0, 25.0, 25.0, 25.0])

    mesh_1d, z = build_1d_tensor_mesh(
        start,
        {"n_cells": 4, "depth_max": 100.0, "growth_factor": 1.0},
        TensorMesh,
    )
    assert mesh_1d.shape_cells == (4,)
    np.testing.assert_allclose(z, [12.5, 37.5, 62.5, 87.5])

    mesh_3d, centers = build_3d_tensor_mesh(
        [0.0, 100.0],
        [10.0, 20.0],
        {
            "nx": 2,
            "ny": 2,
            "nz": 2,
            "depth_max": 100.0,
            "x_pad": 0.0,
            "y_pad": 0.0,
        },
        TensorMesh,
    )
    assert mesh_3d.shape_cells == (2, 2, 2)
    assert set(centers) == {"x", "y", "z", "z_depth"}


def test_fd2d_grid_builder_uses_padding_and_station_offset():
    start = StartingModel([100.0, 200.0], [50.0])
    grid, core_shape = build_fd2d_grid(
        start,
        [10.0, 1010.0],
        {"nx": 2, "n_pad": 1, "x_margin": 0.0, "x_max": 1000.0},
        Grid2D,
    )
    assert core_shape == (2, 2)
    assert grid.resistivity.shape == (3, 4)
    np.testing.assert_allclose(
        core_rho_from_start(start, 2, 2), [[100.0, 100.0], [200.0, 200.0]]
    )
    assert hasattr(grid, "_pycsamt_x_offset")


def test_regularization_helpers_from_config_and_residual():
    cfg = InversionConfig(
        regularization="blocky",
        backend_options={
            "alpha_x": 2.0,
            "alpha_z": 3.0,
            "regularization_weight": 4.0,
            "reference_weight": 0.5,
            "lam": 7.0,
        },
    )
    reg = regularization_from_config(cfg)
    assert reg.kind == "blocky"
    assert reg.alpha_x == 2.0
    assert reg.alpha_z == 3.0
    assert regularization_weight(cfg) == 4.0
    assert pygimli_lambda(cfg) == 7.0

    values = np.array([[1.0, 2.0], [1.5, 3.0]])
    residual = regularization_residual(
        values,
        reference=np.ones_like(values),
        regularization=reg,
        axes=("x", "z"),
    )
    assert residual.size == 8
    assert np.all(np.isfinite(residual))


def test_error_model_component_floors_and_masks():
    cfg = InversionConfig(
        error_floor=0.05,
        phase_error=2.0,
        backend_options={
            "error_model": {
                "rho_absolute": 3.0,
                "phase_relative": 0.01,
                "tdem_absolute": 1e-10,
                "masks": {"station": [True, False]},
            }
        },
    )
    model = error_model_from_config(cfg)
    assert isinstance(model, ErrorModel)
    np.testing.assert_allclose(
        component_errors([10.0, 100.0], cfg, component="rho"),
        [3.0, 5.0],
    )
    np.testing.assert_allclose(
        component_errors([100.0, 400.0], cfg, component="phase"),
        [2.0, 4.0],
    )
    np.testing.assert_allclose(
        component_errors([1e-12, 1e-8], cfg, component="tdem"),
        [1e-10, 5e-10],
    )
    mask = component_mask(np.ones((2, 3)), cfg, component="rho")
    assert mask.shape == (2, 3)
    assert mask[0].all()
    assert not mask[1].any()


def test_simpeg_error_model_honors_masks_and_absolute_phase():
    data = EMData.coerce(
        {
            "freqs": [1.0],
            "rho_a": [[100.0], [200.0]],
            "phase": [[45.0], [50.0]],
        }
    )
    cfg = InversionConfig(
        error_floor=0.05,
        phase_error=2.0,
        backend_options={"masks": {"station": [True, False]}},
    )
    values, errors = simpeg_backend._pack_nsem_observations(data, cfg)
    np.testing.assert_allclose(values, [100.0, 200.0, 45.0, 50.0])
    assert errors[0] == 5.0
    assert errors[1] > 1e20
    assert errors[2] == 2.0
    assert errors[3] > 1e20


def test_pygimli_tdem_error_model_relative_floor():
    data = EMData.coerce(
        {
            "method": "tdem",
            "times": [1e-4, 1e-3],
            "values": [1e-8, 1e-9],
        }
    )
    cfg = InversionConfig(
        method="tdem",
        backend_options={"tdem_relative": 0.1, "tdem_absolute": 1e-10},
    )
    errors = pygimli_backend._tdem_errors(data, cfg)
    np.testing.assert_allclose(errors, [0.1, 0.1])


def test_simpeg_regularization_helper_uses_shared_alphas():
    class FakeReg:
        def __init__(self, mesh, mapping=None):
            self.mesh = mesh
            self.mapping = mapping
            self.alpha_s = None
            self.alpha_x = None
            self.alpha_y = None
            self.alpha_z = None
            self.reference_model = None

    class Modules:
        class regularization:
            WeightedLeastSquares = FakeReg

    class Mesh:
        nC = 3

    cfg = InversionConfig(
        regularization="smooth",
        backend_options={
            "alpha_s": 0.1,
            "alpha_x": 2.0,
            "alpha_y": 2.5,
            "alpha_z": 3.0,
            "reference_model": [100.0, 200.0, 400.0],
        },
    )
    reg = simpeg_backend._build_regularization(Mesh(), object(), Modules, cfg)
    assert reg.alpha_s == 0.1
    assert reg.alpha_x == 2.0
    assert reg.alpha_y == 2.5
    assert reg.alpha_z == 3.0
    np.testing.assert_allclose(
        reg.reference_model, np.log(1.0 / np.array([100.0, 200.0, 400.0]))
    )


def test_builtin_1d_regularization_metadata_smoke(tmp_path):
    freqs = np.logspace(-1, 1, 5)
    response = MT1DForward(freqs=freqs).run(
        LayeredModel([100.0, 300.0], [500.0])
    )
    cfg = InversionConfig(
        data={
            "freqs": freqs,
            "rho_a": response.rho_a,
            "phase": response.phase,
        },
        n_layers=2,
        max_iter=2,
        regularization="damped",
        backend_options={"regularization_weight": 0.1, "alpha_s": 0.5},
        reference_model=StartingModel([120.0, 250.0], [600.0]),
        workdir=str(tmp_path / "reg"),
    )
    result = InversionWorkflow(cfg).run()
    assert result.metadata["regularization"]["kind"] == "damped"
    assert result.metadata["regularization_weight"] == 0.1


def test_modem_backend_dry_run_command_ready(tmp_path):
    for name in ("m0.ws", "data.dat", "control.inv", "covariance.cov"):
        (tmp_path / name).write_text("", encoding="utf-8")
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        backend_options={
            "binary_3d": "Mod3DMT",
            "files": {
                "model": "m0.ws",
                "data": "data.dat",
                "control": "control.inv",
                "covariance": "covariance.cov",
            },
        },
    )
    result = InversionWorkflow(cfg).run()
    assert result.status in {"ready", "loaded"}
    assert result.metadata["command"] is not None
    assert "Mod3DMT" in result.metadata["command"]
    assert result.metadata["run_external"] is False


def test_modem_backend_reports_missing_runner_files(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="modem",
        workdir=str(tmp_path),
        backend_options={
            "files": {
                "model": "m0.ws",
                "data": "data.dat",
                "control": "control.inv",
                "covariance": "covariance.cov",
            },
        },
    )
    result = InversionWorkflow(cfg).run()
    assert result.status == "prepared"
    assert any("missing files" in warning for warning in result.warnings)
    assert result.metadata["command"] is None


def test_modem_file_validation_helpers(tmp_path):
    model = tmp_path / "m0.ws"
    model.write_text("", encoding="utf-8")
    files = {
        "model": str(model),
        "data": str(tmp_path / "data.dat"),
        "control": str(tmp_path / "control.inv"),
    }
    missing = modem_backend._missing_runner_files(
        files,
        required=("model", "data", "control"),
    )
    assert missing == [
        f"data={tmp_path / 'data.dat'}",
        f"control={tmp_path / 'control.inv'}",
    ]
    assert modem_backend._relative_to_workdir(str(model), tmp_path) == "m0.ws"


def test_occam2d_backend_dry_run_command_ready(tmp_path):
    for name in (
        "OccamDataFile.dat",
        "Occam2DMesh",
        "Occam2DModel",
        "Startup",
    ):
        (tmp_path / name).write_text("", encoding="utf-8")
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        backend_options={"binary_name": "Occam2D"},
    )
    result = InversionWorkflow(cfg).run()
    assert result.status in {"ready", "loaded"}
    assert result.metadata["command"] == "Occam2D Startup"
    assert result.metadata["run_external"] is False
    assert result.files["startup"].endswith("Startup")


def test_occam2d_backend_reports_missing_runner_files(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="occam2d",
        workdir=str(tmp_path),
        backend_options={
            "files": {
                "data": "OccamDataFile.dat",
                "mesh": "Occam2DMesh",
                "model": "Occam2DModel",
                "startup": "Startup",
            },
        },
    )
    result = InversionWorkflow(cfg).run()
    assert result.status == "prepared"
    assert any("missing files" in warning for warning in result.warnings)
    assert result.metadata["command"] is None


def test_occam2d_file_validation_helpers(tmp_path):
    startup = tmp_path / "Startup"
    startup.write_text("", encoding="utf-8")
    files = {
        "data": str(tmp_path / "OccamDataFile.dat"),
        "mesh": str(tmp_path / "Occam2DMesh"),
        "model": str(tmp_path / "Occam2DModel"),
        "startup": str(startup),
    }
    missing = occam2d_backend._missing_runner_files(files)
    assert missing == [
        f"data={tmp_path / 'OccamDataFile.dat'}",
        f"mesh={tmp_path / 'Occam2DMesh'}",
        f"model={tmp_path / 'Occam2DModel'}",
    ]
    assert (
        occam2d_backend._command_string("Occam2D", "Startup")
        == "Occam2D Startup"
    )


def test_pygimli_backend_missing_dependency_is_clear():
    if importlib.util.find_spec("pygimli") is not None:
        pytest.skip(
            "pyGIMLi is installed; missing-dependency path not applicable."
        )
    cfg = InversionConfig(
        method="tdem",
        dimension="1d",
        backend="pygimli",
        data={"times": [1e-4], "values": [1e-9]},
    )
    with pytest.raises(ImportError, match="pyGIMLi"):
        InversionWorkflow(cfg).run()


def test_pygimli_observation_packing_helpers():
    data = EMData.coerce(
        {
            "freqs": [1.0, 10.0],
            "rho_a": [100.0, 200.0],
            "phase": [45.0, 50.0],
        }
    )
    values = pygimli_backend._pack_mt_observations(data)
    errors = pygimli_backend._pack_mt_errors(
        data,
        InversionConfig(error_floor=0.05, phase_error=2.0),
    )
    np.testing.assert_allclose(values, [100.0, 200.0, 45.0, 50.0])
    np.testing.assert_allclose(errors[:2], [0.05, 0.05])
    assert np.all(errors[2:] > 0)


def test_pygimli_mt_operator_fallbacks_select_block_model():
    class FakeEM:
        class MT1dBlockModelling:
            def __init__(self, **kwargs):
                if "nLayers" not in kwargs:
                    raise TypeError("nLayers required")
                self.kwargs = kwargs

    start = StartingModel([100.0, 200.0, 300.0], [50.0, 150.0])
    fop, start_model, operator_name, parameterization = (
        pygimli_backend._make_mt_forward_operator(
            FakeEM,
            np.array([1.0, 0.1]),
            start,
            {"mt_operator": ["Missing", "MT1dBlockModelling"]},
        )
    )
    assert operator_name == "MT1dBlockModelling"
    assert parameterization == "thickness_resistivity"
    assert fop.kwargs["nLayers"] == 3
    np.testing.assert_allclose(
        start_model, [50.0, 150.0, 100.0, 200.0, 300.0]
    )


def test_pygimli_tdem_operator_constructor_fallback():
    class FakeEM:
        class TDEMSmoothModelling:
            def __init__(self, **kwargs):
                if "t" not in kwargs:
                    raise TypeError("old pyGIMLi uses t")
                self.kwargs = kwargs

    fop, operator_name = pygimli_backend._make_tdem_forward_operator(
        FakeEM,
        np.array([1e-4, 1e-3]),
        np.array([100.0]),
        {},
        tx_area=10.0,
        rx_area=None,
    )
    assert operator_name == "TDEMSmoothModelling"
    np.testing.assert_allclose(fop.kwargs["t"], [1e-4, 1e-3])
    assert fop.kwargs["txArea"] == 10.0


def test_pygimli_inversion_run_fallback_uses_relative_error():
    class FakeInversion:
        def run(self, observed, **kwargs):
            if "errorVals" in kwargs:
                raise TypeError("old pyGIMLi uses relativeError")
            assert "relativeError" in kwargs
            return np.asarray(kwargs["startModel"], dtype=float) + 1.0

    recovered = pygimli_backend._run_inversion(
        FakeInversion(),
        np.array([1.0, 2.0]),
        startModel=np.array([10.0, 20.0]),
        errorVals=np.array([0.1, 0.1]),
        lam=20.0,
        maxIter=3,
        verbose=False,
    )
    np.testing.assert_allclose(recovered, [11.0, 21.0])


def test_expected_phase_tree_modules_exist():
    import importlib

    modules = [
        "pycsamt.inversion.mt.inv1d",
        "pycsamt.inversion.mt.inv2d",
        "pycsamt.inversion.mt.inv3d",
        "pycsamt.inversion.amt.inv1d",
        "pycsamt.inversion.amt.inv2d",
        "pycsamt.inversion.emap.inv2d",
        "pycsamt.inversion.tdem.inv1d",
        "pycsamt.inversion.tdem.inv2d",
    ]
    for module in modules:
        assert importlib.import_module(module)


def test_method_wrapper_defaults_do_not_duplicate_constructor_logic():
    wrappers = [
        (MT1DInversion, "mt", "1d", "builtin"),
        (MT2DInversion, "mt", "2d", "builtin"),
        (MT3DInversion, "mt", "3d", "modem"),
        (AMT1DInversion, "amt", "1d", "builtin"),
        (AMT2DInversion, "amt", "2d", "builtin"),
        (EMAP2DInversion, "emap", "2d", "builtin"),
        (TDEM1DInversion, "tdem", "1d", "builtin"),
        (TDEM2DInversion, "tdem", "2d", "builtin"),
    ]
    for cls, method, dimension, backend in wrappers:
        inv = cls(max_iter=1)
        assert inv.config.method == method
        assert inv.config.dimension == dimension
        assert inv.config.backend == backend


def test_builtin_1d_inversion_smoke(tmp_path):
    freqs = np.logspace(-2, 2, 12)
    true_model = LayeredModel([80.0, 20.0, 500.0], [250.0, 900.0])
    response = MT1DForward(freqs=freqs).run(true_model)
    data = {
        "freqs": freqs,
        "rho_a": response.rho_a,
        "phase": response.phase,
    }
    start = StartingModel([100.0, 50.0, 300.0], [200.0, 1000.0])
    cfg = InversionConfig(
        data=data,
        n_layers=3,
        starting_model=start,
        max_iter=25,
        workdir=str(tmp_path / "builtin"),
    )
    result = InversionWorkflow(cfg).run()
    assert result.backend == "builtin"
    assert result.model.n_layers == 3
    assert np.isfinite(result.rms)
    assert result.to_resistivity_model().rho_2d.shape == (3, 1)
    assert result.uncertainty is not None
    assert result.uncertainty.model_std.shape == (3, 1)
    assert result.uncertainty.confidence.shape == (3, 1)
    assert np.nanmin(result.uncertainty.confidence) >= 0.0
    assert np.nanmax(result.uncertainty.confidence) <= 1.0
    assert result.history is not None
    assert result.history.n_records >= 1
    arrays = result.history.arrays()
    assert "objective" in arrays
    assert "rms" in arrays


def test_mt_convenience_wrapper(tmp_path):
    freqs = np.logspace(-1, 1, 5)
    response = MT1DForward(freqs=freqs).run(
        LayeredModel([100.0, 300.0], [500.0])
    )
    inv = MT1DInversion(
        {"freqs": freqs, "rho_a": response.rho_a, "phase": response.phase},
        n_layers=2,
        max_iter=5,
        workdir=str(tmp_path / "mt"),
    )
    result = inv.run()
    assert result.method == "mt"
    assert result.dimension == "1d"


def test_builtin_tdem_1d_inversion_smoke(tmp_path):
    times = np.logspace(-5, -3, 4)
    true_model = LayeredModel([80.0, 300.0], [120.0])
    fwd_options = {"loop_radius": 25.0, "n_freqs": 8, "n_lam": 12}
    response = TEM1DForward(times=times, **fwd_options).run(true_model)
    inv = TDEM1DInversion(
        {
            "times": times,
            "values": response.dBz_dt,
        },
        n_layers=2,
        starting_model=StartingModel([100.0, 250.0], [150.0]),
        max_iter=2,
        backend_options=fwd_options,
        workdir=str(tmp_path / "tdem"),
    )
    result = inv.run()
    assert result.method == "tdem"
    assert result.dimension == "1d"
    assert result.backend == "builtin"
    assert result.model.n_layers == 2
    assert np.isfinite(result.rms)


def test_builtin_2d_stitched_profile_and_exports(tmp_path):
    import json
    import zipfile

    freqs = np.logspace(-2, 2, 10)
    models = [
        LayeredModel([80.0, 20.0, 500.0], [250.0, 900.0]),
        LayeredModel([100.0, 35.0, 600.0], [300.0, 850.0]),
        LayeredModel([120.0, 60.0, 700.0], [350.0, 800.0]),
    ]
    responses = [MT1DForward(freqs=freqs).run(model) for model in models]
    data = {
        "method": "mt",
        "freqs": freqs,
        "rho_a": np.vstack([resp.rho_a for resp in responses]),
        "phase": np.vstack([resp.phase for resp in responses]),
        "station_names": ["S1", "S2", "S3"],
        "station_x": [0.0, 500.0, 1000.0],
    }
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="builtin",
        data=data,
        n_layers=3,
        starting_model=StartingModel([100.0, 50.0, 500.0], [300.0, 900.0]),
        max_iter=15,
        workdir=str(tmp_path / "profile"),
    )
    result = InversionWorkflow(cfg).run()
    model = result.to_resistivity_model()
    assert result.dimension == "2d"
    assert model.rho_2d.shape == (3, 3)
    assert model.station_names == ["S1", "S2", "S3"]
    assert np.isfinite(result.rms)

    csv_path = export.to_csv(result, tmp_path / "profile.csv")
    npz_path = export.to_npz(result, tmp_path / "profile.npz")
    geojson_path = export.to_geojson(result, tmp_path / "profile.geojson")
    vtk_path = export.to_vtk(result, tmp_path / "profile.vtk")
    archive_path = export.to_archive(result, tmp_path / "profile.zip")
    assert csv_path.exists()
    assert npz_path.exists()
    assert geojson_path.exists()
    assert vtk_path.exists()
    assert archive_path.exists()
    with np.load(npz_path) as loaded:
        assert "uncertainty_confidence" in loaded
        assert (
            loaded["uncertainty_confidence"].shape
            == result.uncertainty.confidence.shape
        )
        assert "history_objective" in loaded
    with geojson_path.open(encoding="utf-8") as fh:
        geojson = json.load(fh)
    assert geojson["type"] == "FeatureCollection"
    assert len(geojson["features"]) == model.rho_2d.size
    assert "uncertainty_confidence" in geojson["features"][0]["properties"]
    vtk_text = vtk_path.read_text(encoding="utf-8")
    assert "RECTILINEAR_GRID" in vtk_text
    assert "SCALARS rho_log10_ohm_m" in vtk_text
    with zipfile.ZipFile(archive_path) as zf:
        assert {"metadata.json", "result.npz", "model.csv"}.issubset(
            set(zf.namelist())
        )

    ax_model = plot.plot_model(result, section="compact", colorbar=False)
    ax_rms = plot.plot_rms(result)
    assert ax_model.get_ylabel() == "Depth (m)"
    assert ax_rms.get_ylabel() == "Weighted RMS"


def test_geotiff_export_requires_rasterio_or_writes(tmp_path):
    result = InversionResult(
        method="mt",
        dimension="2d",
        backend="builtin",
        model={
            "rho_2d": np.array([[2.0, 2.1], [2.2, 2.3]]),
            "x_centers": np.array([0.0, 100.0]),
            "z_centers": np.array([25.0, 75.0]),
        },
    )
    out = tmp_path / "model.tif"
    if importlib.util.find_spec("rasterio") is None:
        with pytest.raises(ImportError, match="rasterio"):
            export.to_geotiff(result, out)
    else:
        assert export.to_geotiff(result, out).exists()


def test_builtin_2d_fd_physics_smoke(tmp_path):
    freqs = np.array([1.0, 10.0])
    grid = Grid2D.halfspace(
        rho=100.0,
        nx=2,
        nz=2,
        x_max=1000.0,
        z_max=1000.0,
        n_pad=0,
        n_stations=2,
    )
    response = MT2DForward(freqs, grid, verbose=False).run()
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="builtin",
        data={
            "freqs": freqs,
            "rho_a": response.rho_a_te.T,
            "phase": response.phase_te.T,
            "station_x": [0.0, 1000.0],
            "station_names": ["S1", "S2"],
        },
        n_layers=2,
        starting_model=StartingModel([100.0, 100.0], [500.0]),
        max_iter=1,
        backend_options={
            "profile_mode": "fd2d",
            "nx": 2,
            "n_pad": 0,
            "x_margin": 0.0,
            "x_max": 1000.0,
            "components": ("te",),
            "regularization_weight": 0.0,
            "forward_verbose": False,
        },
        workdir=str(tmp_path / "fd2d"),
    )
    result = InversionWorkflow(cfg).run()
    assert result.metadata["profile_mode"] == "fd2d"
    assert result.mesh.metadata["engine"] == "builtin_fd2d"
    assert result.model["rho_2d"].shape == (2, 2)
    assert result.uncertainty is not None
    assert result.uncertainty.sensitivity.shape == (2, 2)
    assert result.history is not None
    assert result.history.metadata["mode"] == "fd2d"
    assert np.isfinite(result.rms)


def test_inversion_uncertainty_result_container_coerces_dict():
    result = InversionResult(
        method="mt",
        dimension="1d",
        backend="builtin",
        uncertainty={
            "model_std": [[0.1], [0.2]],
            "confidence": [[0.9], [0.8]],
            "station_confidence": [0.85],
            "depth_confidence": [0.9, 0.8],
        },
    )
    assert isinstance(result.uncertainty, InversionUncertainty)
    assert result.uncertainty.model_std.shape == (2, 1)


def test_inversion_history_container_arrays_and_result_coercion():
    result = InversionResult(
        method="mt",
        dimension="1d",
        backend="builtin",
        history={
            "records": [
                {
                    "iteration": 0,
                    "phi_d": 4.0,
                    "phi_m": 1.0,
                    "objective": 5.0,
                },
                {
                    "iteration": 1,
                    "phi_d": 2.0,
                    "phi_m": 0.5,
                    "objective": 2.5,
                },
            ],
            "metadata": {"type": "test"},
        },
    )
    assert isinstance(result.history, InversionHistory)
    arrays = result.history.arrays()
    np.testing.assert_allclose(arrays["phi_d"], [4.0, 2.0])
    np.testing.assert_allclose(arrays["objective"], [5.0, 2.5])

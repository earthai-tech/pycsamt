"""Round-trip tests for the first Occam1D implementation layer."""

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DConfig,
    Occam1DData,
    Occam1DModel,
    Occam1DStartup,
    detect_file_type,
)


def test_data_roundtrip_preserves_missing_values(tmp_path):
    data = Occam1DData(
        [100.0, 10.0, 1.0],
        [10.0, np.nan, 100.0],
        [40.0, 45.0, np.nan],
        [0.05, np.nan, 0.1],
        [2.0, 2.0, np.nan],
        station="S01",
    )
    path = data.write(tmp_path / "Data")
    restored = Occam1DData.read(path)
    np.testing.assert_allclose(restored.frequency, data.frequency)
    np.testing.assert_allclose(
        restored.resistivity, data.resistivity, equal_nan=True
    )
    np.testing.assert_allclose(restored.phase, data.phase, equal_nan=True)
    assert restored.station == "S01"
    assert detect_file_type(path) == "data"


def test_model_and_startup_roundtrip(tmp_path):
    model = Occam1DModel.build(8, 5.0, 2000.0, resistivity=80.0)
    assert model.thickness[0] == pytest.approx(5.0)
    assert model.depth[-1] == pytest.approx(2000.0)
    assert np.all(np.diff(model.thickness[:-1]) >= 0)
    model_path = model.write(tmp_path / "Model")
    restored_model = Occam1DModel.read(model_path)
    np.testing.assert_allclose(restored_model.depth, model.depth)
    np.testing.assert_allclose(restored_model.resistivity, model.resistivity)
    assert detect_file_type(model_path) == "model"

    cfg = Occam1DConfig(n_layers=8, starting_resistivity=80.0)
    startup = Occam1DStartup.from_model(model, cfg)
    startup_path = startup.write(tmp_path / "Startup")
    restored = Occam1DStartup.read(startup_path)
    np.testing.assert_allclose(restored.parameters, startup.parameters)
    assert restored.iteration == 0
    assert detect_file_type(startup_path) == "startup"


def test_invalid_frequency_is_rejected():
    with pytest.raises(ValueError, match="frequency"):
        Occam1DData([0.0], [10.0], [45.0], [0.05], [2.0])


def test_config_rejects_inconsistent_depths():
    with pytest.raises(ValueError, match="depth_max"):
        Occam1DConfig(first_thickness=10.0, depth_max=5.0).validate()

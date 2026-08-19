"""Tests for Occam1D sounding extraction and input generation."""

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DBatch,
    Occam1DConfig,
    Occam1DData,
    Occam1DInputBuilder,
    Occam1DModel,
    Occam1DStartup,
)


class SiteStub:
    def __init__(self, name="S01"):
        self.name = name
        self.freq = np.array([1.0, 100.0, 10.0])
        self.z = np.zeros((3, 2, 2), dtype=complex)
        self.z[:, 0, 1] = [1 + 1j, 2 + 2j, 3 + 3j]
        self.z[:, 1, 0] = -self.z[:, 0, 1]
        self.z_err = np.full((3, 2, 2), 0.01)
        self.rho = None
        self.phase = None


def test_single_builder_writes_complete_input_set(tmp_path):
    cfg = Occam1DConfig(
        mode="xy",
        n_layers=6,
        first_thickness=5.0,
        depth_max=1000.0,
        freq_min=5.0,
    )
    builder = Occam1DInputBuilder(
        SiteStub(), tmp_path / "run", cfg
    ).build()
    assert builder.is_ready
    assert builder.data.frequency.tolist() == [100.0, 10.0]
    assert isinstance(
        Occam1DData.read(tmp_path / "run" / cfg.data_file), Occam1DData
    )
    assert isinstance(
        Occam1DModel.read(tmp_path / "run" / cfg.model_file), Occam1DModel
    )
    assert isinstance(
        Occam1DStartup.read(tmp_path / "run" / cfg.startup_file),
        Occam1DStartup,
    )


def test_determinant_builder_uses_impedance_tensor(tmp_path):
    builder = Occam1DInputBuilder(
        SiteStub(),
        tmp_path,
        Occam1DConfig(n_layers=5, depth_max=500.0),
    ).build()
    assert np.all(np.isfinite(builder.data.resistivity))
    assert np.all(builder.data.resistivity_error >= 0.05)


def test_multiple_sites_require_selection(tmp_path):
    with pytest.raises(ValueError, match="multiple sites"):
        Occam1DInputBuilder([SiteStub("A"), SiteStub("B")], tmp_path).build()


def test_batch_creates_station_directories(tmp_path):
    batch = Occam1DBatch(
        [SiteStub("A 01"), SiteStub("B/02")],
        tmp_path,
        Occam1DConfig(mode="yx", n_layers=5, depth_max=500.0),
    ).build_all()
    assert batch.is_ready
    assert (tmp_path / "A_01" / "Startup").is_file()
    assert (tmp_path / "B_02" / "Startup").is_file()


def _small_batch(tmp_path):
    return Occam1DBatch(
        [SiteStub("A01"), SiteStub("B02")],
        tmp_path,
        Occam1DConfig(
            mode="yx", n_layers=5, depth_max=500.0, max_iterations=3
        ),
    ).build_all()


def test_invert_all_sequential_returns_a_summary_per_station(tmp_path):
    batch = _small_batch(tmp_path)
    outcome = batch.invert_all(n_jobs=1, export_text=False)
    assert set(outcome["results"]) == {"A01", "B02"}
    assert outcome["failures"] == {}
    for summary in outcome["results"].values():
        assert summary["iterations"] >= 1
        assert (tmp_path / summary["station"] / "native-restart.json").is_file()


def test_invert_all_parallel_matches_sequential(tmp_path):
    pytest.importorskip("joblib")
    sequential = _small_batch(tmp_path / "seq").invert_all(
        n_jobs=1, export_text=False
    )
    parallel = _small_batch(tmp_path / "par").invert_all(
        n_jobs=2, export_text=False
    )
    for station in sequential["results"]:
        assert (
            sequential["results"][station]["final_rms"]
            == pytest.approx(
                parallel["results"][station]["final_rms"], rel=1e-8
            )
        )


def test_invert_all_requires_built_stations(tmp_path):
    batch = Occam1DBatch(
        SiteStub("A01"),
        tmp_path,
        Occam1DConfig(mode="yx", n_layers=5, depth_max=500.0),
    )
    with pytest.raises(RuntimeError, match="no completed station builders"):
        batch.invert_all()

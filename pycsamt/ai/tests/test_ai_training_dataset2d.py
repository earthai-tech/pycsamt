"""Contracts for :mod:`pycsamt.ai.training.dataset2d`.

Both ``zxy`` (TE mode) and ``zyx`` (TM mode) are validated by this
module; see its module docstring for the convergence evidence behind
that claim, and :mod:`pycsamt.forward.tests.test_em2d` for the
underlying ``_assemble_tm`` fix this depends on.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.geology import GeologyGrid
from pycsamt.ai.training.dataset2d import (
    Maxwell2DDataset,
    Maxwell2DDatasetConfig,
    Maxwell2DSample,
    build_2d_maxwell_problem,
    generate_2d_maxwell_dataset,
)
from pycsamt.forward.maxwell.benchmarks import half_space_impedance


def _small_grid():
    return GeologyGrid.regular_2d(nx=6, nz=4, dx_m=300, dz_m=150)


def _config(**overrides):
    values = dict(
        dataset_id="demo-2d-v1",
        grid=_small_grid(),
        correlation_length_x_m=(600.0, 900.0),
        correlation_length_z_m=(150.0, 300.0),
        frequencies_hz=[10.0, 3.0],
        station_x_m=[600.0, 900.0, 1200.0],
        n_realizations=2,
        seed=0,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    values.update(overrides)
    return Maxwell2DDatasetConfig(**values)


def test_config_defaults_to_both_modes():
    config = _config()
    assert config.components == ("zxy", "zyx")


def test_build_problem_accepts_an_explicit_geological_model():
    grid = _small_grid()
    model = np.full(grid.shape, 100.0)
    problem = build_2d_maxwell_problem(
        grid,
        model,
        [10.0, 3.0],
        [600.0, 900.0],
        station_names=["A", "B"],
        metadata={"family": "layered"},
    )
    assert problem.components == ("zxy", "zyx")
    assert problem.receivers.names == ("A", "B")
    assert problem.metadata["family"] == "layered"
    assert problem.conductivity_s_m.shape == problem.mesh.shape


def test_build_problem_rejects_invalid_model_and_station_names():
    grid = _small_grid()
    with pytest.raises(ValueError, match="resistivity_ohm_m"):
        build_2d_maxwell_problem(
            grid,
            np.ones((2, 2)),
            [10.0],
            [600.0],
        )
    with pytest.raises(ValueError, match="station_names"):
        build_2d_maxwell_problem(
            grid,
            np.full(grid.shape, 100.0),
            [10.0],
            [600.0, 900.0],
            station_names=["duplicate", "duplicate"],
        )


def test_config_rejects_non_2d_grid_and_bad_surface():
    grid_3d = GeologyGrid.regular_3d(
        nx=2, ny=2, nz=2, dx_m=1, dy_m=1, dz_m=1
    )
    with pytest.raises(TypeError, match="2-D GeologyGrid"):
        _config(grid=grid_3d)
    deep_grid = GeologyGrid.regular_2d(
        nx=6, nz=4, dx_m=300, dz_m=150, z_origin_m=10.0
    )
    with pytest.raises(ValueError, match="depth 0"):
        _config(grid=deep_grid)


def test_config_rejects_bad_correlation_ranges():
    with pytest.raises(ValueError, match="0 < low <= high"):
        _config(correlation_length_x_m=(900.0, 600.0))
    with pytest.raises(ValueError, match="0 < low <= high"):
        _config(correlation_length_z_m=(0.0, 100.0))


def test_config_rejects_bad_frequencies_and_stations():
    with pytest.raises(ValueError, match="frequencies_hz"):
        _config(frequencies_hz=[1.0, 1.0])
    with pytest.raises(ValueError, match="station_x_m"):
        _config(station_x_m=[-500.0])
    with pytest.raises(ValueError, match="station_x_m"):
        _config(station_x_m=[600.0, 600.0])


def test_config_rejects_bad_components_realizations_seed():
    with pytest.raises(ValueError, match="components"):
        _config(components=("zzz",))
    with pytest.raises(ValueError, match="n_realizations"):
        _config(n_realizations=0)
    with pytest.raises(ValueError, match="seed"):
        _config(seed=-1)


def test_config_rejects_bad_std_and_max_mesh_cells():
    with pytest.raises(ValueError, match="log_resistivity_std"):
        _config(log_resistivity_std=0.0)
    with pytest.raises(ValueError, match="mesh_safety_factor"):
        _config(mesh_safety_factor=0.0)
    with pytest.raises(ValueError, match="max_mesh_cells"):
        _config(max_mesh_cells=1)


def test_to_dict_round_trips_key_fields():
    config = _config()
    record = config.to_dict()
    assert record["dataset_id"] == "demo-2d-v1"
    assert record["components"] == ["zxy", "zyx"]
    assert record["max_mesh_cells"] == config.max_mesh_cells


def test_generate_dataset_end_to_end_both_modes():
    config = _config()
    dataset = generate_2d_maxwell_dataset(config)
    assert isinstance(dataset, Maxwell2DDataset)
    assert dataset.rejected == ()
    assert len(dataset.samples) == 2
    sample = dataset.samples[0]
    assert isinstance(sample, Maxwell2DSample)
    assert sample.survey.shape == (3, 2, 2)
    assert sample.survey.components == ("zxy", "zyx")
    assert np.all(np.isfinite(sample.resistivity_ohm_m))
    assert np.all(sample.resistivity_ohm_m > 0)
    assert sample.relative_residual < 1e-6
    assert sample.mesh_cells > 0


def test_generate_dataset_te_mode_only_still_available():
    config = _config(components=("zxy",))
    dataset = generate_2d_maxwell_dataset(config)
    sample = dataset.samples[0]
    assert sample.survey.shape == (3, 2, 1)
    assert sample.survey.components == ("zxy",)


def test_generate_dataset_is_seed_deterministic():
    first = generate_2d_maxwell_dataset(_config(seed=7))
    second = generate_2d_maxwell_dataset(_config(seed=7))
    np.testing.assert_array_equal(
        first.samples[0].resistivity_ohm_m,
        second.samples[0].resistivity_ohm_m,
    )
    np.testing.assert_array_equal(
        first.samples[0].survey.impedance,
        second.samples[0].survey.impedance,
    )
    assert first.samples[0].seed == second.samples[0].seed


def test_generate_dataset_different_seeds_differ():
    first = generate_2d_maxwell_dataset(_config(seed=1))
    second = generate_2d_maxwell_dataset(_config(seed=2))
    assert not np.array_equal(
        first.samples[0].resistivity_ohm_m,
        second.samples[0].resistivity_ohm_m,
    )


def test_generate_dataset_split_covers_every_sample():
    config = _config(
        n_realizations=6, validation_fraction=0.2, test_fraction=0.2
    )
    dataset = generate_2d_maxwell_dataset(config)
    ids = {s.realization_id for s in dataset.samples}
    split_ids = (
        set(dataset.split.train)
        | set(dataset.split.validation)
        | set(dataset.split.test)
    )
    assert ids == split_ids
    assert dataset.manifest.sample_count == len(dataset.samples)


def test_dataset_select_partitions_are_disjoint():
    config = _config(
        n_realizations=6, validation_fraction=0.2, test_fraction=0.2
    )
    dataset = generate_2d_maxwell_dataset(config)
    train_ids = {s.realization_id for s in dataset.select("train")}
    test_ids = {s.realization_id for s in dataset.select("test")}
    assert train_ids.isdisjoint(test_ids)
    assert (
        len(dataset.select("train"))
        + len(dataset.select("validation"))
        + len(dataset.select("test"))
        == len(dataset.samples)
    )


def test_max_mesh_cells_guard_raises_a_clear_error():
    config = _config(
        grid=GeologyGrid.regular_2d(nx=6, nz=4, dx_m=1.0, dz_m=150),
        station_x_m=[3.0],
        mesh_safety_factor=8.0,
        max_mesh_cells=10,
    )
    with pytest.raises(ValueError, match="max_mesh_cells"):
        generate_2d_maxwell_dataset(config)


def test_both_modes_match_analytic_half_space_response():
    grid = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=200, dz_m=100)
    config = Maxwell2DDatasetConfig(
        dataset_id="halfspace-check",
        grid=grid,
        correlation_length_x_m=(1.0, 2.0),
        correlation_length_z_m=(1.0, 2.0),
        frequencies_hz=[10.0, 1.0],
        station_x_m=[400.0],
        n_realizations=1,
        seed=0,
        log_resistivity_mean=2.0,
        log_resistivity_std=1e-6,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    dataset = generate_2d_maxwell_dataset(config)
    sample = dataset.samples[0]
    analytic = half_space_impedance(100.0, config.frequencies_hz)
    zxy = sample.survey.impedance[0, :, 0]
    zyx = sample.survey.impedance[0, :, 1]
    te_error = np.abs(zxy - analytic) / np.abs(analytic)
    # Zyx = -Zxy for an isotropic half-space (standard MT sign convention).
    tm_error = np.abs(zyx - (-analytic)) / np.abs(analytic)
    assert np.all(te_error < 0.05)
    assert np.all(tm_error < 0.05)

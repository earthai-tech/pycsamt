"""Contracts for :mod:`pycsamt.ai.training.dataset3d`.

Mirrors :mod:`pycsamt.ai.tests.test_ai_training_dataset2d`, but for the
genuinely 3-D (:class:`~pycsamt.forward.maxwell.mt3d.MT3DAdapter`)
path — see that module's docstring for why a padded, non-uniform
solver mesh is used instead of a uniform one within the small research
cell budget.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.geology import GeologyGrid
from pycsamt.ai.training.dataset3d import (
    Maxwell3DDataset,
    Maxwell3DDatasetConfig,
    Maxwell3DSample,
    generate_3d_maxwell_dataset,
)
from pycsamt.forward.maxwell.benchmarks import half_space_impedance


def _small_grid():
    return GeologyGrid.regular_3d(nx=4, ny=4, nz=6, dx_m=200, dy_m=200, dz_m=100)


def _config(**overrides):
    values = dict(
        dataset_id="demo-3d-v1",
        grid=_small_grid(),
        correlation_length_x_m=(400.0, 800.0),
        correlation_length_y_m=(400.0, 800.0),
        correlation_length_z_m=(100.0, 200.0),
        frequencies_hz=[50.0, 20.0],
        station_xy_m=[[300.0, 300.0], [500.0, 500.0], [700.0, 700.0]],
        n_realizations=2,
        seed=0,
        mesh_safety_factor=2.0,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    values.update(overrides)
    return Maxwell3DDatasetConfig(**values)


def test_config_defaults_to_te_like_mode():
    config = _config()
    assert config.components == ("zxy", "zyx")


def test_config_rejects_non_3d_grid_and_bad_surface():
    grid_2d = GeologyGrid.regular_2d(nx=4, nz=4, dx_m=1, dz_m=1)
    with pytest.raises(TypeError, match="3-D GeologyGrid"):
        _config(grid=grid_2d)
    deep_grid = GeologyGrid.regular_3d(
        nx=4, ny=4, nz=6, dx_m=200, dy_m=200, dz_m=100, z_origin_m=10.0
    )
    with pytest.raises(ValueError, match="depth 0"):
        _config(grid=deep_grid)


def test_config_rejects_bad_correlation_ranges():
    with pytest.raises(ValueError, match="0 < low <= high"):
        _config(correlation_length_x_m=(900.0, 600.0))
    with pytest.raises(ValueError, match="0 < low <= high"):
        _config(correlation_length_y_m=(0.0, 100.0))
    with pytest.raises(ValueError, match="0 < low <= high"):
        _config(correlation_length_z_m=(0.0, 100.0))


def test_config_rejects_bad_frequencies_and_stations():
    with pytest.raises(ValueError, match="frequencies_hz"):
        _config(frequencies_hz=[1.0, 1.0])
    with pytest.raises(ValueError, match="station_xy_m"):
        _config(station_xy_m=[[-500.0, 300.0]])
    with pytest.raises(ValueError, match="station_xy_m"):
        _config(station_xy_m=[[300.0, 300.0], [300.0, 300.0]])


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
    with pytest.raises(ValueError, match="cells_per_skin_depth"):
        _config(cells_per_skin_depth=0.0)


def test_to_dict_round_trips_key_fields():
    config = _config()
    record = config.to_dict()
    assert record["dataset_id"] == "demo-3d-v1"
    assert record["components"] == ["zxy", "zyx"]
    assert record["max_mesh_cells"] == config.max_mesh_cells
    assert record["cells_per_skin_depth"] == config.cells_per_skin_depth


def test_generate_dataset_end_to_end():
    config = _config()
    dataset = generate_3d_maxwell_dataset(config)
    assert isinstance(dataset, Maxwell3DDataset)
    assert dataset.rejected == ()
    assert len(dataset.samples) == 2
    sample = dataset.samples[0]
    assert isinstance(sample, Maxwell3DSample)
    assert sample.survey.shape == (3, 2, 2)
    assert sample.survey.components == ("zxy", "zyx")
    assert np.all(np.isfinite(sample.resistivity_ohm_m))
    assert np.all(sample.resistivity_ohm_m > 0)
    assert sample.mesh_cells > 0


def test_generate_dataset_te_like_mode_only_still_available():
    config = _config(components=("zxy",))
    dataset = generate_3d_maxwell_dataset(config)
    sample = dataset.samples[0]
    assert sample.survey.shape == (3, 2, 1)
    assert sample.survey.components == ("zxy",)


def test_generate_dataset_is_seed_deterministic():
    first = generate_3d_maxwell_dataset(_config(seed=7))
    second = generate_3d_maxwell_dataset(_config(seed=7))
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
    first = generate_3d_maxwell_dataset(_config(seed=1))
    second = generate_3d_maxwell_dataset(_config(seed=2))
    assert not np.array_equal(
        first.samples[0].resistivity_ohm_m,
        second.samples[0].resistivity_ohm_m,
    )


def test_generate_dataset_split_covers_every_sample():
    config = _config(
        n_realizations=6, validation_fraction=0.2, test_fraction=0.2
    )
    dataset = generate_3d_maxwell_dataset(config)
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
    dataset = generate_3d_maxwell_dataset(config)
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
    config = _config(mesh_safety_factor=8.0, max_mesh_cells=50)
    with pytest.raises(ValueError, match="max_mesh_cells"):
        generate_3d_maxwell_dataset(config)


def test_matches_analytic_half_space_response():
    grid = GeologyGrid.regular_3d(nx=4, ny=4, nz=4, dx_m=200, dy_m=200, dz_m=100)
    config = Maxwell3DDatasetConfig(
        dataset_id="halfspace-check",
        grid=grid,
        correlation_length_x_m=(1.0, 2.0),
        correlation_length_y_m=(1.0, 2.0),
        correlation_length_z_m=(1.0, 2.0),
        # Frequencies matched to MT3DAdapter's own validated range (see
        # forward/tests/test_maxwell_mt3d.py's calibrated mesh, <=2 Hz).
        # cells_per_skin_depth (see test_higher_cells_per_skin_depth_
        # improves_high_frequency_accuracy below) narrows but does not
        # fully close the higher-frequency gap within a small cell
        # budget -- see the module docstring's caveat.
        frequencies_hz=[2.0, 1.0],
        station_xy_m=[[400.0, 400.0]],
        n_realizations=1,
        seed=0,
        log_resistivity_mean=2.0,
        log_resistivity_std=1e-6,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    dataset = generate_3d_maxwell_dataset(config)
    sample = dataset.samples[0]
    analytic = half_space_impedance(100.0, config.frequencies_hz)
    zxy = sample.survey.impedance[0, :, 0]
    zyx = sample.survey.impedance[0, :, 1]
    te_error = np.abs(zxy - analytic) / np.abs(analytic)
    # Zyx = -Zxy for an isotropic half-space (standard MT sign convention).
    tm_error = np.abs(zyx - (-analytic)) / np.abs(analytic)
    assert np.all(te_error < 0.05)
    assert np.all(tm_error < 0.05)


def _high_frequency_half_space_error(cells_per_skin_depth: float) -> np.ndarray:
    grid = GeologyGrid.regular_3d(nx=4, ny=4, nz=4, dx_m=200, dy_m=200, dz_m=100)
    config = Maxwell3DDatasetConfig(
        dataset_id="hf-check",
        grid=grid,
        correlation_length_x_m=(1.0, 2.0),
        correlation_length_y_m=(1.0, 2.0),
        correlation_length_z_m=(1.0, 2.0),
        frequencies_hz=[50.0, 20.0],
        station_xy_m=[[400.0, 400.0]],
        n_realizations=1,
        seed=0,
        log_resistivity_mean=2.0,
        log_resistivity_std=1e-6,
        cells_per_skin_depth=cells_per_skin_depth,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    dataset = generate_3d_maxwell_dataset(config)
    sample = dataset.samples[0]
    analytic = half_space_impedance(100.0, config.frequencies_hz)
    zxy = sample.survey.impedance[0, :, 0]
    return np.abs(zxy - analytic) / np.abs(analytic)


def test_higher_cells_per_skin_depth_improves_high_frequency_accuracy():
    """Regression guard for the frequency-aware core-resolution fix:
    a higher cells_per_skin_depth (finer core, same domain-extent
    safety factor) must reduce -- not just change -- the 20-50 Hz
    half-space error, even though it does not fully close the gap
    to the <=2 Hz range's accuracy within a small cell budget (see
    the module docstring's "Frequency-dependent mesh accuracy"
    section for the measured numbers this guards).
    """
    coarse_error = _high_frequency_half_space_error(4.0)
    default_error = _high_frequency_half_space_error(8.0)
    assert np.all(default_error < coarse_error)

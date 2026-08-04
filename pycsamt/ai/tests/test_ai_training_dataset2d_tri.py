"""Contracts for :mod:`pycsamt.ai.training.dataset2d_tri`.

``Mare2DEMAdapter`` needs a real compiled MARE2DEM binary to produce
genuine physics (see its own module docstring), which does not exist
in this environment. Mesh generation and geological-field sampling
need no external binary at all and are tested directly; the full
``generate_2d_tri_maxwell_dataset()`` pipeline is exercised end-to-end
against a scripted stand-in executable (mirrors
``test_maxwell_mare2dem.py``'s technique), proving the data plumbing
is correct without claiming real Maxwell physics.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from pycsamt.ai.training.dataset2d_tri import (
    MaxwellTri2DDataset,
    MaxwellTri2DDatasetConfig,
    MaxwellTri2DSample,
    generate_2d_tri_maxwell_dataset,
)
from pycsamt.forward.maxwell import ExternalRunPolicy
from pycsamt.forward.maxwell.mare2dem import Mare2DEMAdapter
from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh
from pycsamt.models.mare2dem.config import Mare2DEMConfig


def _config(**overrides):
    values = dict(
        dataset_id="demo-2dtri-v1",
        x_range_m=(0.0, 1000.0),
        z_range_m=(0.0, 500.0),
        correlation_length_x_m=(200.0, 400.0),
        correlation_length_z_m=(100.0, 200.0),
        frequencies_hz=[10.0, 1.0],
        station_x_m=[200.0, 500.0, 800.0],
        n_realizations=2,
        seed=0,
        mesh_target_cell_m=250.0,
        field_grid_cell_m=100.0,
        validation_fraction=0.0,
        test_fraction=0.0,
    )
    values.update(overrides)
    return MaxwellTri2DDatasetConfig(**values)


# ---------------------------------------------------------------------------
# build_graded_tri_mesh (re-exercised here against this module's own config
# defaults; the exhaustive unit tests live in
# pycsamt/forward/tests/test_maxwell_tri_mesh_gen.py)
# ---------------------------------------------------------------------------


def test_build_graded_tri_mesh_places_stations_exactly_on_nodes():
    mesh = build_graded_tri_mesh(
        (0.0, 1000.0), (0.0, 500.0), [200.0, 500.0, 800.0], surface_cell_m=100.0
    )
    for station_x in (200.0, 500.0, 800.0):
        assert np.any(
            np.isclose(mesh.nodes_m[:, 0], station_x)
            & np.isclose(mesh.nodes_m[:, 1], 0.0)
        )


def test_build_graded_tri_mesh_assigns_one_region_per_triangle():
    mesh = build_graded_tri_mesh(
        (0.0, 1000.0), (0.0, 500.0), [500.0], surface_cell_m=100.0
    )
    assert sorted(mesh.region_ids.tolist()) == list(
        range(1, mesh.n_triangles + 1)
    )


def test_build_graded_tri_mesh_rejects_bad_ranges():
    with pytest.raises(ValueError, match="increasing"):
        build_graded_tri_mesh(
            (100.0, 0.0), (0.0, 500.0), [50.0], surface_cell_m=100.0
        )


def test_build_graded_tri_mesh_grows_with_depth():
    mesh = build_graded_tri_mesh(
        (0.0, 1000.0), (0.0, 500.0), [200.0, 500.0, 800.0], surface_cell_m=30.0
    )
    centroids = mesh.triangle_centroids_m
    areas = mesh.triangle_areas_m2
    shallow = areas[centroids[:, 1] < 40.0]
    deep = areas[centroids[:, 1] > 400.0]
    assert shallow.size > 0 and deep.size > 0
    assert deep.mean() > shallow.mean()


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------


def test_config_defaults_to_both_modes():
    assert _config().components == ("zxy", "zyx")


def test_config_rejects_nonzero_surface():
    with pytest.raises(ValueError, match="depth 0"):
        _config(z_range_m=(10.0, 500.0))


def test_config_rejects_bad_correlation_ranges():
    with pytest.raises(ValueError, match="0 < low"):
        _config(correlation_length_x_m=(400.0, 100.0))


def test_config_rejects_stations_outside_domain():
    with pytest.raises(ValueError, match="station_x_m"):
        _config(station_x_m=[200.0, 5000.0])


def test_config_rejects_bad_components():
    with pytest.raises(ValueError, match="components"):
        _config(components=("zxx",))


def test_to_dict_round_trips_key_fields():
    config = _config()
    data = config.to_dict()
    assert data["dataset_id"] == "demo-2dtri-v1"
    assert data["x_range_m"] == [0.0, 1000.0]
    assert data["topo_x_m"] is None and data["topo_z_m"] is None


def test_config_accepts_topo_without_requiring_flat_zero_surface():
    # z_range_m[0]=10 (not 0) would normally be rejected -- but the
    # "must be depth 0" rule only applies to the flat-surface default.
    config = _config(
        z_range_m=(0.0, 500.0),
        topo_x_m=[0.0, 500.0, 1000.0],
        topo_z_m=[-20.0, -50.0, -10.0],
    )
    assert config.topo_x_m.tolist() == [0.0, 500.0, 1000.0]
    data = config.to_dict()
    assert data["topo_z_m"] == [-20.0, -50.0, -10.0]


def test_config_rejects_topo_given_only_partially():
    with pytest.raises(ValueError, match="topo_x_m and topo_z_m"):
        _config(topo_x_m=[0.0, 1000.0])


def test_config_rejects_mismatched_topo_arrays():
    with pytest.raises(ValueError, match="topo_x_m/topo_z_m"):
        _config(topo_x_m=[0.0, 500.0, 1000.0], topo_z_m=[0.0, 0.0])


# ---------------------------------------------------------------------------
# generate_2d_tri_maxwell_dataset: scripted stand-in end-to-end
# ---------------------------------------------------------------------------

_FAKE_MARE2DEM_SCRIPT = textwrap.dedent(
    """
    import sys
    from pathlib import Path

    stem = sys.argv[-1]
    workdir = Path.cwd()
    lines = (workdir / f"{stem}.emdata").read_text().splitlines()

    data_lines = []
    in_data = False
    for line in lines:
        s = line.strip()
        if not s or s.startswith("!"):
            continue
        if s.lower().startswith("# data"):
            in_data = True
            continue
        if in_data:
            data_lines.append(s)

    out = ["Format:  EMResp_2.3\\n", f"# Data:       {len(data_lines)}\\n"]
    for line in data_lines:
        parts = line.split()
        code, freq_idx, tx_idx, rx_idx = (int(v) for v in parts[:4])
        data_val, std_err = float(parts[4]), float(parts[5])
        predicted = 1.0 if code in (113, 115) else 1.0
        out.append(
            f"{code:7d} {freq_idx:7d} {tx_idx:7d} {rx_idx:7d} "
            f"{data_val:22.15g} {std_err:22.15g} {predicted:20.15g} 0.0\\n"
        )

    (workdir / f"{stem}.0.resp").write_text("".join(out))
    """
)


def _fake_adapter(tmp_path: Path) -> Mare2DEMAdapter:
    script = tmp_path / "fake_mare2dem.py"
    script.write_text(_FAKE_MARE2DEM_SCRIPT)
    adapter = Mare2DEMAdapter(
        config=Mare2DEMConfig(use_mpi=False),
        run_policy=ExternalRunPolicy(sys.executable, workdir=str(tmp_path)),
    )
    adapter._build_command = lambda problem, workdir, executable, context: [
        str(executable),
        str(script),
        context["stem"],
    ]
    return adapter


def test_generate_dataset_end_to_end_against_scripted_stand_in(tmp_path):
    dataset = generate_2d_tri_maxwell_dataset(
        _config(), adapter=_fake_adapter(tmp_path)
    )
    assert isinstance(dataset, MaxwellTri2DDataset)
    assert len(dataset.samples) >= 1
    assert dataset.rejected == ()
    for sample in dataset.samples:
        assert isinstance(sample, MaxwellTri2DSample)
        assert sample.resistivity_ohm_m.shape == (dataset.mesh.n_triangles,)
        assert sample.survey.shape == (3, 2, 2)


def test_generate_dataset_with_topo_places_receivers_on_the_terrain(tmp_path):
    topo_x = [0.0, 250.0, 500.0, 750.0, 1000.0]
    topo_z = [0.0, -30.0, -60.0, -30.0, 0.0]
    dataset = generate_2d_tri_maxwell_dataset(
        _config(topo_x_m=topo_x, topo_z_m=topo_z),
        adapter=_fake_adapter(tmp_path),
    )
    assert len(dataset.samples) >= 1
    expected_z = np.interp([200.0, 500.0, 800.0], topo_x, topo_z)
    assert not np.any(np.isclose(expected_z, 0.0))
    coords = dataset.samples[0].survey.coordinates_m
    np.testing.assert_allclose(coords[:, 2], expected_z)


def test_dataset_split_covers_every_sample(tmp_path):
    dataset = generate_2d_tri_maxwell_dataset(
        _config(n_realizations=4, validation_fraction=0.25),
        adapter=_fake_adapter(tmp_path),
    )
    sample_ids = {s.realization_id for s in dataset.samples}
    split_ids = (
        set(dataset.split.train)
        | set(dataset.split.validation)
        | set(dataset.split.test)
    )
    assert sample_ids == split_ids


def test_generate_dataset_is_seed_deterministic(tmp_path):
    config = _config()
    d1 = generate_2d_tri_maxwell_dataset(config, adapter=_fake_adapter(tmp_path))
    d2 = generate_2d_tri_maxwell_dataset(config, adapter=_fake_adapter(tmp_path))
    r1 = sorted(d1.samples, key=lambda s: s.realization_id)
    r2 = sorted(d2.samples, key=lambda s: s.realization_id)
    for s1, s2 in zip(r1, r2):
        np.testing.assert_array_equal(
            s1.resistivity_ohm_m, s2.resistivity_ohm_m
        )


def test_dataset_select_partitions_are_disjoint(tmp_path):
    dataset = generate_2d_tri_maxwell_dataset(
        _config(
            n_realizations=6, validation_fraction=0.2, test_fraction=0.2
        ),
        adapter=_fake_adapter(tmp_path),
    )
    train_ids = {s.realization_id for s in dataset.select("train")}
    val_ids = {s.realization_id for s in dataset.select("validation")}
    test_ids = {s.realization_id for s in dataset.select("test")}
    assert not (train_ids & val_ids)
    assert not (train_ids & test_ids)
    assert not (val_ids & test_ids)

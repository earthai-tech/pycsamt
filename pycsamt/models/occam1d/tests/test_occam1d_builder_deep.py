"""Transactional and batch-isolation tests for Occam1D input builders."""

import hashlib
import json

import numpy as np
import pytest

from pycsamt.models.occam1d import (
    Occam1DBatch,
    Occam1DBuildState,
    Occam1DConfig,
    Occam1DInputBuilder,
)


class Site:
    def __init__(self, name="S01", *, valid=True):
        self.name = name
        self.freq = np.array([100.0, 10.0]) if valid else None
        self.rho = np.full((2, 2, 2), 100.0) if valid else None
        self.phase = np.full((2, 2, 2), 45.0) if valid else None


def _config(**overrides):
    values = {
        "mode": "xy",
        "n_layers": 4,
        "first_thickness": 5.0,
        "depth_max": 100.0,
    }
    values.update(overrides)
    return Occam1DConfig(**values)


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_prepare_has_no_filesystem_side_effect_and_copies_config(tmp_path):
    original = _config(mode="xy")
    workdir = tmp_path / "run"
    builder = Occam1DInputBuilder(Site(), workdir, original)

    returned = builder.prepare(mode="yx")

    assert returned is builder
    assert builder.build_state is Occam1DBuildState.PREPARED
    assert not workdir.exists()
    assert original.mode == "xy"
    assert builder.config.mode == "yx"
    assert builder.data is not None
    assert builder.model is not None
    assert builder.startup is not None


def test_commit_writes_verified_manifest_and_binds_paths(tmp_path):
    builder = Occam1DInputBuilder(
        Site(), tmp_path / "run", _config()
    ).build()

    assert builder.is_ready
    assert builder.is_bound
    assert builder.build_state is Occam1DBuildState.READY
    manifest = json.loads(builder.paths.manifest.read_text(encoding="utf8"))
    assert manifest["schema"].endswith("/v1")
    assert manifest["station"] == "S01"
    assert manifest["data"]["n_observations"] == 4
    assert manifest["model"]["n_layers"] == 4
    assert manifest["files"][builder.config.data_file] == _digest(
        builder.paths.data
    )
    assert manifest["files"][builder.config.model_file] == _digest(
        builder.paths.model
    )
    assert manifest["files"][builder.config.startup_file] == _digest(
        builder.paths.startup
    )
    assert builder.data.path == builder.paths.data
    assert builder.model.path == builder.paths.model
    assert builder.startup.path == builder.paths.startup
    assert builder.diagnostics()["is_ready"] is True
    assert "state      : ready" in builder.summary()


def test_repeated_identical_build_has_deterministic_manifest(tmp_path):
    workdir = tmp_path / "run"
    first = Occam1DInputBuilder(Site(), workdir, _config()).build()
    text = first.paths.manifest.read_text(encoding="utf8")

    second = Occam1DInputBuilder(Site(), workdir, _config()).build()

    assert second.paths.manifest.read_text(encoding="utf8") == text


def test_overwrite_false_refuses_existing_managed_files(tmp_path):
    workdir = tmp_path / "run"
    Occam1DInputBuilder(Site(), workdir, _config()).build()
    builder = Occam1DInputBuilder(Site(), workdir, _config()).prepare()

    with pytest.raises(FileExistsError, match="Managed output"):
        builder.commit(overwrite=False)
    assert builder.build_state is Occam1DBuildState.PREPARED


def test_commit_requires_preparation(tmp_path):
    builder = Occam1DInputBuilder(Site(), tmp_path / "run", _config())

    with pytest.raises(RuntimeError, match="prepare"):
        builder.commit()


def test_station_selection_by_name_and_index(tmp_path):
    sites = [Site("A"), Site("B")]
    by_name = Occam1DInputBuilder(
        sites, tmp_path / "name", _config(), station="b"
    ).build()
    by_index = Occam1DInputBuilder(
        sites, tmp_path / "index", _config(), station=0
    ).build()

    assert by_name.data.station == "B"
    assert by_index.data.station == "A"
    with pytest.raises(IndexError, match="outside"):
        Occam1DInputBuilder(
            sites, tmp_path / "bad", _config(), station=2
        ).prepare()
    with pytest.raises(ValueError, match="non-negative"):
        Occam1DInputBuilder(
            sites, tmp_path / "bad", _config(), station=-1
        )


def test_duplicate_station_names_are_ambiguous(tmp_path):
    builder = Occam1DInputBuilder(
        [Site("A"), Site("A")],
        tmp_path / "run",
        _config(),
        station="A",
    )

    with pytest.raises(ValueError, match="not unique"):
        builder.prepare()


def test_reserved_manifest_filename_is_rejected_before_write(tmp_path):
    builder = Occam1DInputBuilder(
        Site(),
        tmp_path / "run",
        _config(data_file="occam1d_manifest.json"),
    )

    with pytest.raises(ValueError, match="reserved manifest"):
        builder.prepare()
    assert not (tmp_path / "run").exists()


def test_unknown_override_does_not_mutate_builder_config(tmp_path):
    builder = Occam1DInputBuilder(Site(), tmp_path / "run", _config())
    before = builder.config.to_dict()

    with pytest.raises(TypeError, match="Unknown"):
        builder.prepare(does_not_exist=1)
    assert builder.config.to_dict() == before
    assert builder.build_state is Occam1DBuildState.FAILED


def test_failed_commit_restores_every_prior_managed_file(
    tmp_path, monkeypatch
):
    workdir = tmp_path / "run"
    existing = Occam1DInputBuilder(Site(), workdir, _config()).build()
    paths = (
        existing.paths.data,
        existing.paths.model,
        existing.paths.startup,
        existing.paths.manifest,
    )
    before = {path: path.read_bytes() for path in paths}
    replacement = Occam1DInputBuilder(
        Site(), workdir, _config(starting_resistivity=250.0)
    ).prepare()

    real_replace = __import__("os").replace
    calls = {"count": 0}

    def fail_second(source, destination):
        calls["count"] += 1
        if calls["count"] == 2:
            raise OSError("simulated commit failure")
        return real_replace(source, destination)

    monkeypatch.setattr("os.replace", fail_second)
    with pytest.raises(OSError, match="simulated"):
        replacement.commit()

    assert replacement.build_state is Occam1DBuildState.FAILED
    assert {path: path.read_bytes() for path in paths} == before


def test_batch_configuration_is_isolated_per_child(tmp_path):
    original = _config(mode="xy")
    batch = Occam1DBatch(
        [Site("A"), Site("B")], tmp_path, original
    ).build_all(mode="yx")

    assert batch.is_ready
    assert original.mode == "xy"
    assert batch.config.mode == "yx"
    assert all(builder.config is not batch.config for builder in batch.builders)
    assert all(builder.data.mode == "yx" for builder in batch.builders)
    assert batch.diagnostics()["completed"] == 2
    assert "completed=2" in batch.summary()


def test_batch_rejects_sanitized_station_collision_before_writes(tmp_path):
    batch = Occam1DBatch(
        [Site("A/B"), Site("A B")], tmp_path, _config()
    )

    with pytest.raises(ValueError, match="collision"):
        batch.build_all()
    assert not (tmp_path / "A_B").exists()


def test_batch_can_continue_and_report_station_failure(tmp_path):
    batch = Occam1DBatch(
        [Site("good"), Site("bad", valid=False)], tmp_path, _config()
    ).build_all(continue_on_error=True)

    assert not batch.is_ready
    assert batch.build_state is Occam1DBuildState.FAILED
    assert len(batch.builders) == 1
    assert "bad" in batch.failures
    assert (tmp_path / "good" / "Startup").is_file()
    assert not (tmp_path / "bad").exists()


def test_batch_export_requires_completed_builders(tmp_path):
    batch = Occam1DBatch([], tmp_path, _config())

    with pytest.raises(RuntimeError, match="no completed"):
        batch.export_images()

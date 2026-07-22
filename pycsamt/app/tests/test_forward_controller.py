# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ForwardController — no Qt required."""

from __future__ import annotations

import json

import numpy as np
import pytest

from pycsamt.app.desktop.controllers.forward_controller import (
    GEOLOGY_PRESET_NAMES,
    ForwardController,
)

_MOD = "pycsamt.app.desktop.controllers.forward_controller"


@pytest.fixture(autouse=True)
def isolated_library(monkeypatch, tmp_path):
    """Redirect the on-disk library path away from the real user home."""
    monkeypatch.setattr(f"{_MOD}._LIBRARY_PATH", tmp_path / "forward_models.json")
    return tmp_path / "forward_models.json"


def _record(name="m1"):
    return {
        "dim": "1d",
        "resistivity": [100.0, 10.0, 500.0],
        "thickness": [300.0, 800.0],
    }


# ── Construction / load_library ─────────────────────────────────────────────


def test_init_with_no_file_starts_empty(isolated_library):
    assert not isolated_library.exists()
    ctrl = ForwardController()
    assert ctrl.model_names == []


def test_init_loads_existing_file(isolated_library):
    isolated_library.parent.mkdir(parents=True, exist_ok=True)
    isolated_library.write_text(
        json.dumps({"models": [{"name": "existing", "dim": "1d"}]}),
        encoding="utf-8",
    )
    ctrl = ForwardController()
    assert ctrl.model_names == ["existing"]


def test_load_library_missing_models_key_defaults_empty(isolated_library):
    isolated_library.parent.mkdir(parents=True, exist_ok=True)
    isolated_library.write_text(json.dumps({}), encoding="utf-8")
    ctrl = ForwardController()
    assert ctrl.model_names == []


def test_load_library_corrupt_json_starts_empty(isolated_library, caplog):
    isolated_library.parent.mkdir(parents=True, exist_ok=True)
    isolated_library.write_text("{not valid json!!", encoding="utf-8")
    ctrl = ForwardController()
    assert ctrl.model_names == []


def test_load_library_can_be_called_again_to_refresh(isolated_library):
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())

    # Simulate an external process appending another model to the file.
    on_disk = json.loads(isolated_library.read_text(encoding="utf-8"))
    on_disk["models"].append({"name": "m2", "dim": "1d"})
    isolated_library.write_text(json.dumps(on_disk), encoding="utf-8")

    ctrl.load_library()
    assert set(ctrl.model_names) == {"m1", "m2"}


# ── save_library ─────────────────────────────────────────────────────────────


def test_save_library_writes_json_round_trip(isolated_library):
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())

    assert isolated_library.exists()
    on_disk = json.loads(isolated_library.read_text(encoding="utf-8"))
    assert on_disk["models"][0]["name"] == "m1"


def test_save_library_creates_parent_dir(tmp_path, monkeypatch):
    nested = tmp_path / "nested" / "dir" / "forward_models.json"
    monkeypatch.setattr(f"{_MOD}._LIBRARY_PATH", nested)
    assert not nested.parent.exists()

    ctrl = ForwardController()
    ctrl.save_model("m1", _record())
    assert nested.exists()


def test_save_library_swallows_write_errors(monkeypatch, isolated_library):
    ctrl = ForwardController()

    def boom(self, *a, **k):
        raise OSError("disk full")

    # pathlib.Path instances disallow per-instance attribute assignment
    # (Path.write_text is a read-only slot-backed attribute), so patch the
    # method at the class level for the duration of this test instead.
    monkeypatch.setattr("pathlib.Path.write_text", boom)
    # Must not raise, even though the underlying write fails.
    ctrl.save_model("m1", _record())


# ── CRUD: save_model / get_model_record / model_names ───────────────────────


def test_save_model_adds_new_record():
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())
    assert ctrl.model_names == ["m1"]
    rec = ctrl.get_model_record("m1")
    assert rec["name"] == "m1"
    assert rec["dim"] == "1d"
    assert rec["resistivity"] == [100.0, 10.0, 500.0]


def test_save_model_sets_name_key_even_if_absent_in_record():
    ctrl = ForwardController()
    record = {"dim": "1d", "resistivity": [1.0], "thickness": []}
    assert "name" not in record
    ctrl.save_model("solo", record)
    assert ctrl.get_model_record("solo")["name"] == "solo"


def test_save_model_overwrites_existing_by_name():
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())
    new_record = {"dim": "1d", "resistivity": [1.0, 2.0], "thickness": [5.0]}
    ctrl.save_model("m1", new_record)

    assert ctrl.model_names == ["m1"]  # not duplicated
    rec = ctrl.get_model_record("m1")
    assert rec["resistivity"] == [1.0, 2.0]


def test_save_model_does_not_mutate_caller_dict():
    ctrl = ForwardController()
    record = _record()
    ctrl.save_model("m1", record)
    assert "name" not in record  # save_model copies before mutating


def test_save_model_persists_to_disk(isolated_library):
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())

    reloaded = ForwardController()
    assert reloaded.model_names == ["m1"]


def test_model_names_lists_all_in_order():
    ctrl = ForwardController()
    ctrl.save_model("a", _record())
    ctrl.save_model("b", _record())
    ctrl.save_model("c", _record())
    assert ctrl.model_names == ["a", "b", "c"]


def test_get_model_record_missing_returns_none():
    ctrl = ForwardController()
    assert ctrl.get_model_record("nope") is None


# ── delete_model ──────────────────────────────────────────────────────────────


def test_delete_model_removes_existing_returns_true():
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())
    assert ctrl.delete_model("m1") is True
    assert ctrl.model_names == []


def test_delete_model_missing_returns_false():
    ctrl = ForwardController()
    assert ctrl.delete_model("nope") is False


def test_delete_model_persists_removal(isolated_library):
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())
    ctrl.delete_model("m1")

    reloaded = ForwardController()
    assert reloaded.model_names == []


def test_delete_model_only_removes_matching_name():
    ctrl = ForwardController()
    ctrl.save_model("m1", _record())
    ctrl.save_model("m2", _record())
    ctrl.delete_model("m1")
    assert ctrl.model_names == ["m2"]


# ── rename_model ──────────────────────────────────────────────────────────────


def test_rename_model_existing_returns_true():
    ctrl = ForwardController()
    ctrl.save_model("old", _record())
    assert ctrl.rename_model("old", "new") is True
    assert ctrl.model_names == ["new"]
    assert ctrl.get_model_record("old") is None
    assert ctrl.get_model_record("new") is not None


def test_rename_model_skips_non_matching_entries_before_match():
    # Exercises the loop continuing past a non-matching name before it
    # finds and renames the matching one further down the library.
    ctrl = ForwardController()
    ctrl.save_model("first", _record())
    ctrl.save_model("target", _record())
    assert ctrl.rename_model("target", "renamed") is True
    assert ctrl.model_names == ["first", "renamed"]


def test_rename_model_missing_returns_false():
    ctrl = ForwardController()
    assert ctrl.rename_model("nope", "whatever") is False


def test_rename_model_persists(isolated_library):
    ctrl = ForwardController()
    ctrl.save_model("old", _record())
    ctrl.rename_model("old", "new")

    reloaded = ForwardController()
    assert reloaded.model_names == ["new"]


def test_rename_model_to_duplicate_name_creates_two_entries_with_same_name():
    # Documents current (unvalidated) behaviour: rename does not check for
    # collisions with an existing name, so model_names can end up with a
    # duplicate. See report: possible bug, not fixed per task instructions.
    ctrl = ForwardController()
    ctrl.save_model("a", _record())
    ctrl.save_model("b", _record())
    assert ctrl.rename_model("a", "b") is True
    assert ctrl.model_names == ["b", "b"]


# ── model_to_record / record_to_arrays ───────────────────────────────────────


def test_model_to_record_serialises_arrays_to_lists():
    resistivity = np.array([100.0, 10.0, 500.0])
    thickness = np.array([300.0, 800.0])
    record = ForwardController.model_to_record(
        "m1", "1d", resistivity, thickness
    )
    assert record == {
        "name": "m1",
        "dim": "1d",
        "resistivity": [100.0, 10.0, 500.0],
        "thickness": [300.0, 800.0],
    }
    assert isinstance(record["resistivity"], list)
    assert isinstance(record["thickness"], list)


def test_record_to_arrays_returns_float_ndarrays():
    record = {
        "name": "m1",
        "dim": "1d",
        "resistivity": [100, 10, 500],
        "thickness": [300, 800],
    }
    res, thick = ForwardController.record_to_arrays(record)
    assert isinstance(res, np.ndarray)
    assert isinstance(thick, np.ndarray)
    assert res.dtype == float
    assert thick.dtype == float
    assert res.shape == (3,)
    assert thick.shape == (2,)
    np.testing.assert_array_equal(res, [100.0, 10.0, 500.0])
    np.testing.assert_array_equal(thick, [300.0, 800.0])


def test_model_to_record_and_record_to_arrays_round_trip():
    resistivity = np.array([50.0, 250.0, 1000.0, 5.0])
    thickness = np.array([100.0, 200.0, 300.0])
    record = ForwardController.model_to_record(
        "rt", "1d", resistivity, thickness
    )
    res2, thick2 = ForwardController.record_to_arrays(record)
    np.testing.assert_array_equal(res2, resistivity)
    np.testing.assert_array_equal(thick2, thickness)


def test_model_to_record_works_via_save_and_get(isolated_library):
    ctrl = ForwardController()
    resistivity = np.array([1.0, 2.0, 3.0])
    thickness = np.array([10.0, 20.0])
    record = ForwardController.model_to_record(
        "roundtrip", "1d", resistivity, thickness
    )
    ctrl.save_model("roundtrip", record)

    fetched = ctrl.get_model_record("roundtrip")
    res, thick = ForwardController.record_to_arrays(fetched)
    np.testing.assert_array_equal(res, resistivity)
    np.testing.assert_array_equal(thick, thickness)


# ── build_preset_1d ───────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "name", ["sedimentary", "porphyry", "mineralized", "marine"]
)
def test_build_preset_1d_returns_resistivity_and_thickness(name):
    preset = ForwardController.build_preset_1d(name, seed=42)
    assert set(preset.keys()) == {"resistivity", "thickness"}
    assert isinstance(preset["resistivity"], list)
    assert isinstance(preset["thickness"], list)
    assert len(preset["resistivity"]) >= 1
    # halfspace layer has no thickness entry
    assert len(preset["thickness"]) == len(preset["resistivity"]) - 1
    assert all(r > 0 for r in preset["resistivity"])


def test_build_preset_1d_deterministic_with_fixed_seed():
    p1 = ForwardController.build_preset_1d("sedimentary", seed=7)
    p2 = ForwardController.build_preset_1d("sedimentary", seed=7)
    assert p1 == p2


def test_build_preset_1d_varies_without_fixed_seed():
    p1 = ForwardController.build_preset_1d("sedimentary", seed=None)
    p2 = ForwardController.build_preset_1d("sedimentary", seed=None)
    # Astronomically unlikely to collide with independent RNG draws.
    assert p1 != p2


def test_build_preset_1d_unknown_name_raises_key_error():
    with pytest.raises(KeyError):
        ForwardController.build_preset_1d("not-a-real-geology", seed=1)


def test_geology_preset_names_all_buildable():
    for name in GEOLOGY_PRESET_NAMES:
        preset = ForwardController.build_preset_1d(name, seed=1)
        assert len(preset["resistivity"]) == len(preset["thickness"]) + 1

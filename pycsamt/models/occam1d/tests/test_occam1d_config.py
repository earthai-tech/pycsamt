"""Scientific and serialization contracts for Occam1DConfig."""

from types import MappingProxyType

import pytest

from pycsamt.models.occam1d import Occam1DConfig


def test_defaults_are_valid_and_mode_is_canonical():
    config = Occam1DConfig(mode=" TE ")
    assert config.mode == "te"
    assert config.canonical_component == "xy"
    assert config.frequency_bounds == (None, None)
    assert config.lagrange_min < config.lagrange_start < config.lagrange_max


def test_layer_growth_factor_satisfies_depth_equation():
    config = Occam1DConfig(
        n_layers=12,
        first_thickness=4.0,
        depth_max=1500.0,
    )
    scale = config.layer_growth_factor
    total = config.first_thickness * sum(
        scale**power for power in range(config.n_layers - 1)
    )
    assert scale > 1.0
    assert total == pytest.approx(config.depth_max)


def test_uniform_layers_have_unit_growth_factor():
    config = Occam1DConfig(
        n_layers=6,
        first_thickness=10.0,
        depth_max=50.0,
    )
    assert config.layer_growth_factor == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"mode": "scalar"}, "mode"),
        ({"error_floor_rho": 1.1}, "error_floor_rho"),
        ({"error_floor_phase": 91.0}, "error_floor_phase"),
        ({"freq_min": 10.0, "freq_max": 1.0}, "freq_min"),
        ({"n_layers": 1}, "n_layers"),
        (
            {"n_layers": 10, "first_thickness": 10.0, "depth_max": 80.0},
            "depth_max",
        ),
        ({"roughness_type": 3}, "roughness_type"),
        ({"lagrange_min": -1.0}, "lagrange_min"),
        (
            {"lagrange_min": 10.0, "lagrange_max": 1.0},
            "lagrange_min",
        ),
        (
            {"lagrange_start": 5.0, "lagrange_max": 4.0},
            "lagrange_start",
        ),
        ({"data_file": "folder/Data"}, "data_file"),
        ({"model_file": "Startup"}, "filenames"),
    ],
)
def test_invalid_scientific_and_file_settings_are_rejected(changes, message):
    with pytest.raises(ValueError, match=message):
        Occam1DConfig(**changes)


@pytest.mark.parametrize(
    "changes",
    [
        {"n_layers": True},
        {"max_iterations": 1.5},
        {"target_misfit": "1.0"},
        {"binary_name": 1},
    ],
)
def test_wrong_scalar_types_are_rejected(changes):
    with pytest.raises(TypeError):
        Occam1DConfig(**changes)


def test_update_is_atomic_when_proposed_values_are_invalid():
    config = Occam1DConfig()
    original = config.to_dict()
    with pytest.raises(ValueError):
        config.update(freq_min=100.0, freq_max=1.0)
    assert config.to_dict() == original


def test_updated_returns_independent_validated_copy():
    original = Occam1DConfig()
    changed = original.updated(mode="YX", target_misfit=1.2)
    assert changed.mode == "yx"
    assert changed.target_misfit == pytest.approx(1.2)
    assert original.mode == "determinant"


def test_run_paths_are_absolute_read_only_and_do_not_create(tmp_path):
    root = tmp_path / "not-created"
    paths = Occam1DConfig().run_paths(root)
    assert isinstance(paths, MappingProxyType)
    assert paths["workdir"].is_absolute()
    assert paths["data"].name == "Occam1DData"
    assert not root.exists()
    with pytest.raises(TypeError):
        paths["extra"] = root / "extra"


@pytest.mark.parametrize("suffix", ["py", "json", "yml"])
def test_documented_templates_roundtrip(tmp_path, suffix):
    config = Occam1DConfig(mode="yx", n_layers=24)
    path = config.to_template(tmp_path / f"occam1d.{suffix}")
    text = path.read_text(encoding="utf8")
    restored = Occam1DConfig.from_file(path)
    assert "Sounding response" in text
    assert restored == config


def test_strict_config_read_rejects_misspelled_keys(tmp_path):
    path = tmp_path / "bad.py"
    path.write_text(
        "CONFIG = {'mode': 'xy', 'target_misfitt': 1.0}\n",
        encoding="utf8",
    )
    with pytest.raises(ValueError, match="target_misfitt"):
        Occam1DConfig.from_file(path)


def test_non_strict_config_read_ignores_unknown_keys(tmp_path):
    path = tmp_path / "loose.py"
    path.write_text(
        "CONFIG = {'mode': 'xy', 'application_note': 'draft'}\n",
        encoding="utf8",
    )
    config = Occam1DConfig.from_file(path, strict=False)
    assert config.mode == "xy"

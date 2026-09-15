"""Branch coverage for reusable model configuration templates."""

from __future__ import annotations

import builtins
import importlib
import sys
import types
from dataclasses import dataclass

import pytest

from pycsamt.models import config_io


@dataclass
class ExampleConfig:
    count: int = 2
    enabled: bool = True
    ratio: float = 1e-5
    label: str = "demo"
    optional: object = None


SCHEMA = [
    config_io.ConfigParameter("count", "A deliberately long count description " * 5, "Numbers"),
    config_io.ConfigParameter("enabled", "Enable the feature.", "Flags"),
]


def test_internal_format_and_value_helpers(tmp_path):
    cfg = ExampleConfig()
    assert config_io._config_dict(cfg)["count"] == 2
    assert config_io._config_dict({"count": 3}) == {"count": 3}
    with pytest.raises(TypeError, match="dataclass"):
        config_io._field_names(dict)
    assert config_io._field_names(ExampleConfig) >= {"count", "enabled"}
    assert config_io._format_name("x.yaml", None) == "yaml"
    assert config_io._format_name("x.unknown", None) == "py"
    assert config_io._format_name("x", ".YAML") == "yml"
    with pytest.raises(ValueError, match="fmt must"):
        config_io._format_name("x", "toml")
    assert config_io._target_path(tmp_path / "x", "json").suffix == ".json"
    assert config_io._target_path(tmp_path / "x.cfg", "json").suffix == ".cfg"
    assert config_io._comment_lines("", "# ") == []
    assert len(config_io._comment_lines("word " * 50, "# ")) > 1
    assert config_io._py_value([1]) == "[1]"
    assert config_io._yaml_value(None) == "null"
    assert config_io._yaml_value(True) == "true"
    assert config_io._yaml_value(False) == "false"
    assert config_io._yaml_value(3) == "3"
    assert config_io._yaml_value(1e-5) == "1.0e-05"
    assert config_io._yaml_value("é") == '"é"'


@pytest.mark.parametrize("fmt,suffix", [("py", ".py"), ("json", ".json"), ("yaml", ".yml")])
def test_template_roundtrips_all_formats(tmp_path, fmt, suffix):
    pytest.importorskip("yaml") if fmt == "yaml" else None
    path = config_io.write_config_template(
        tmp_path / f"nested/config-{fmt}", ExampleConfig(), SCHEMA, fmt=fmt, title="Example"
    )
    assert path.suffix == suffix and path.exists()
    values = config_io.read_config_file(path, ExampleConfig)
    assert values == ExampleConfig().__dict__


def test_python_reader_assignment_variants_and_failures(tmp_path):
    annotated = tmp_path / "annotated.py"
    annotated.write_text("CONFIG: dict = {'count': 4}\n")
    assert config_io.read_config_file(annotated, ExampleConfig) == {"count": 4}

    for name, text in [
        ("missing.py", "VALUE = {}\n"),
        ("notdict.py", "CONFIG = []\n"),
        ("annotated-notdict.py", "CONFIG: list = []\n"),
    ]:
        path = tmp_path / name
        path.write_text(text)
        with pytest.raises(ValueError, match="CONFIG dictionary"):
            config_io.read_config_file(path, ExampleConfig)


def test_json_yaml_validation_and_unknown_keys(tmp_path):
    raw = tmp_path / "raw.json"
    raw.write_text('{"count": 8, "unknown": 9, "_note": "ignored"}')
    with pytest.raises(ValueError, match="unknown"):
        config_io.read_config_file(raw, ExampleConfig)
    assert config_io.read_config_file(raw, ExampleConfig, strict=False) == {"count": 8}

    array = tmp_path / "array.json"
    array.write_text("[]")
    with pytest.raises(ValueError, match="object"):
        config_io.read_config_file(array, ExampleConfig)

    yaml = pytest.importorskip("yaml")
    empty = tmp_path / "empty.yml"
    empty.write_text("")
    assert config_io.read_config_file(empty, ExampleConfig) == {}
    sequence = tmp_path / "sequence.yaml"
    sequence.write_text("- one\n- two\n")
    with pytest.raises(ValueError, match="mapping"):
        config_io.read_config_file(sequence, ExampleConfig)

    unsupported = tmp_path / "config.toml"
    unsupported.write_text("")
    with pytest.raises(ValueError, match="Unsupported"):
        config_io.read_config_file(unsupported, ExampleConfig)


def test_yaml_missing_dependency_error(tmp_path, monkeypatch):
    path = tmp_path / "config.yml"
    path.write_text("count: 1\n")
    real_import = builtins.__import__

    def reject_yaml(name, *args, **kwargs):
        if name == "yaml":
            raise ImportError("missing")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_yaml)
    with pytest.raises(ImportError, match="requires PyYAML"):
        config_io.read_config_file(path, ExampleConfig)


def test_models_package_tolerates_optional_import_failure(monkeypatch):
    import pycsamt.models as models

    real_import_module = importlib.import_module
    original_occam = sys.modules.get("pycsamt.models.occam2d")

    def fail_occam(name, package=None):
        if name == "pycsamt.models.occam2d":
            raise ImportError("optional dependency unavailable")
        return real_import_module(name, package)

    with monkeypatch.context() as patcher:
        patcher.setattr(importlib, "import_module", fail_occam)
        reloaded = importlib.reload(models)
        assert isinstance(reloaded.occam2d, types.ModuleType)
        assert reloaded.occam2d.__name__ == "pycsamt.models.occam2d"

    # Restore the real package so this test cannot affect later collection.
    if original_occam is not None:
        sys.modules["pycsamt.models.occam2d"] = original_occam
    else:
        sys.modules.pop("pycsamt.models.occam2d", None)
    importlib.reload(models)

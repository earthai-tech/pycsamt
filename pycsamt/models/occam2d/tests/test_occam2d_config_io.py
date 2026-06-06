"""Tests for Occam2D source-of-truth configuration files."""

from __future__ import annotations

import pytest


def test_write_and_read_python_template(tmp_path):
    """Python templates should preserve defaults and comments."""
    from pycsamt.models.occam2d.config import OccamConfig

    path = OccamConfig.write_template(tmp_path / "occam2d_config.py")
    text = path.read_text()

    assert "CONFIG = {" in text
    assert "Electromagnetic modes" in text
    cfg = OccamConfig.from_file(path)
    assert cfg.binary_name == "Occam2D"
    assert cfg.modes == ["TE", "TM"]


def test_write_and_read_json_template(tmp_path):
    """JSON templates should roundtrip through the config block."""
    from pycsamt.models.occam2d.config import OccamConfig

    path = OccamConfig.write_template(tmp_path / "occam2d_config.json")
    text = path.read_text()

    assert '"_schema"' in text
    assert '"config"' in text
    cfg = OccamConfig.from_file(path)
    assert cfg.target_misfit == pytest.approx(1.0)


def test_write_and_read_yaml_template(tmp_path):
    """YAML templates should roundtrip editable values."""
    from pycsamt.models.occam2d.config import OccamConfig

    path = OccamConfig.write_template(tmp_path / "occam2d_config.yml")
    text = path.read_text()

    assert "# ---- Data Options ----" in text
    cfg = OccamConfig.from_file(path)
    assert cfg.n_layers == 30


def test_to_template_uses_instance_values(tmp_path):
    """Instance templates should write customized values."""
    from pycsamt.models.occam2d.config import OccamConfig

    cfg = OccamConfig(modes=["TM"], n_layers=44)
    path = cfg.to_template(tmp_path / "custom.py")
    loaded = OccamConfig.read(path)

    assert loaded.modes == ["TM"]
    assert loaded.n_layers == 44


def test_from_file_rejects_unknown_keys(tmp_path):
    """Strict reading should protect users from misspelled keys."""
    from pycsamt.models.occam2d.config import OccamConfig

    path = tmp_path / "bad.py"
    path.write_text("CONFIG = {'modes': ['TE'], 'bad_key': 1}\n")

    with pytest.raises(ValueError, match="bad_key"):
        OccamConfig.from_file(path)


def test_from_file_can_ignore_unknown_keys(tmp_path):
    """Non-strict reading should ignore application metadata."""
    from pycsamt.models.occam2d.config import OccamConfig

    path = tmp_path / "loose.py"
    path.write_text("CONFIG = {'modes': ['TE'], 'bad_key': 1}\n")
    cfg = OccamConfig.from_file(path, strict=False)

    assert cfg.modes == ["TE"]

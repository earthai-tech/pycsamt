# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for pycsamt.models.modem.control.ModEmControl.

test_control_numeric_format.py already covers the float-formatting
regression on ``write()``. This file fills the remaining gap: construction
from a default/custom ``ModEmConfig``, ``from_config``, the full
``read``/``write`` round trip (including the integer attribute and the
tolerant handling of unknown or malformed lines), ``FileNotFoundError``,
and the verbose-logging branch.
"""

from __future__ import annotations

import logging

import pytest

from pycsamt.models.modem.config import ModEmConfig
from pycsamt.models.modem.control import ModEmControl


def test_init_defaults_from_default_config():
    ctrl = ModEmControl()
    cfg = ModEmConfig()
    assert ctrl.output_stem == cfg.output_stem
    assert ctrl.initial_lambda == cfg.initial_lambda
    assert ctrl.lambda_divisor == cfg.lambda_divisor
    assert ctrl.initial_alpha == cfg.initial_alpha
    assert ctrl.rms_diff_tol == cfg.rms_diff_tol
    assert ctrl.target_rms == cfg.target_rms
    assert ctrl.lambda_exit == cfg.lambda_exit
    assert ctrl.max_iterations == cfg.max_iterations
    assert ctrl.config is cfg or isinstance(ctrl.config, ModEmConfig)


def test_init_from_custom_config():
    cfg = ModEmConfig(max_iterations=42, target_rms=1.1, output_stem="run_x")
    ctrl = ModEmControl(config=cfg)
    assert ctrl.max_iterations == 42
    assert ctrl.target_rms == 1.1
    assert ctrl.output_stem == "run_x"
    assert ctrl.config is cfg


def test_from_config_classmethod():
    cfg = ModEmConfig(max_iterations=50, target_rms=1.03)
    ctrl = ModEmControl.from_config(cfg)
    assert isinstance(ctrl, ModEmControl)
    assert ctrl.max_iterations == 50
    assert ctrl.target_rms == 1.03


def test_from_config_none_uses_default():
    ctrl = ModEmControl.from_config(None)
    assert ctrl.max_iterations == ModEmConfig().max_iterations


def test_write_creates_parent_dirs_and_returns_path(tmp_path):
    ctrl = ModEmControl.from_config(ModEmConfig())
    dest = tmp_path / "nested" / "dir" / "control.inv"
    out = ctrl.write(dest)
    assert out == dest
    assert dest.exists()


def test_write_int_attribute_formatting(tmp_path):
    ctrl = ModEmControl.from_config(ModEmConfig(max_iterations=77))
    path = ctrl.write(tmp_path / "control.inv")
    text = path.read_text()
    line = next(
        l for l in text.splitlines() if "Maximum number of iterations" in l
    )
    value = line.rsplit(":", 1)[-1].strip()
    assert value == "77"


def test_write_string_attribute_formatting(tmp_path):
    ctrl = ModEmControl.from_config(ModEmConfig(output_stem="myrun"))
    path = ctrl.write(tmp_path / "control.inv")
    text = path.read_text()
    line = next(
        l for l in text.splitlines() if "Model and data output" in l
    )
    value = line.rsplit(":", 1)[-1].strip()
    assert value == "myrun"


def test_read_round_trip_all_fields(tmp_path):
    cfg = ModEmConfig(
        output_stem="roundtrip",
        initial_lambda=5.0,
        lambda_divisor=20.0,
        initial_alpha=3.0,
        rms_diff_tol=1.0e-3,
        target_rms=1.1,
        lambda_exit=1.0e-4,
        max_iterations=99,
    )
    ctrl = ModEmControl.from_config(cfg)
    path = ctrl.write(tmp_path / "control.inv")

    loaded = ModEmControl.read(path)
    assert loaded.output_stem == "roundtrip"
    assert loaded.initial_lambda == pytest.approx(5.0)
    assert loaded.lambda_divisor == pytest.approx(20.0)
    assert loaded.initial_alpha == pytest.approx(3.0)
    assert loaded.rms_diff_tol == pytest.approx(1.0e-3)
    assert loaded.target_rms == pytest.approx(1.1)
    assert loaded.lambda_exit == pytest.approx(1.0e-4)
    assert loaded.max_iterations == 99


def test_read_missing_file_raises(tmp_path):
    missing = tmp_path / "nope.inv"
    with pytest.raises(FileNotFoundError, match="ModEM control file not found"):
        ModEmControl.read(missing)


def test_read_skips_lines_without_colon(tmp_path):
    path = tmp_path / "control.inv"
    path.write_text(
        "this line has no colon and should be skipped\n"
        "Maximum number of iterations:                     30\n",
        encoding="utf-8",
    )
    ctrl = ModEmControl.read(path)
    assert ctrl.max_iterations == 30


def test_read_ignores_unrecognized_keys(tmp_path):
    path = tmp_path / "control.inv"
    path.write_text(
        "Some unrelated Fortran comment: 123\n"
        "Maximum number of iterations:                     15\n",
        encoding="utf-8",
    )
    ctrl = ModEmControl.read(path)
    assert ctrl.max_iterations == 15


def test_read_malformed_float_keeps_default(tmp_path):
    default_lambda = ModEmConfig().initial_lambda
    path = tmp_path / "control.inv"
    path.write_text(
        "Initial damping factor lambda:                     not_a_number\n",
        encoding="utf-8",
    )
    ctrl = ModEmControl.read(path)
    assert ctrl.initial_lambda == default_lambda


def test_read_malformed_int_keeps_default(tmp_path):
    default_iters = ModEmConfig().max_iterations
    path = tmp_path / "control.inv"
    path.write_text(
        "Maximum number of iterations:                     not_a_number\n",
        encoding="utf-8",
    )
    ctrl = ModEmControl.read(path)
    assert ctrl.max_iterations == default_iters


def test_read_int_attr_accepts_float_string(tmp_path):
    path = tmp_path / "control.inv"
    path.write_text(
        "Maximum number of iterations:                     42.0\n",
        encoding="utf-8",
    )
    ctrl = ModEmControl.read(path)
    assert ctrl.max_iterations == 42


def test_read_verbose_logs(tmp_path, caplog):
    # The "pycsamt" logger tree is configured with `propagate: no` (see
    # pycsamt/log/p.configlog.yml) so records never reach caplog's
    # handler, which pytest attaches only to the root logger. Attach it
    # directly to this class's own logger (name is deterministic: see
    # ModEmBase.__init__) to observe the record.
    cfg = ModEmConfig()
    ctrl = ModEmControl.from_config(cfg)
    path = ctrl.write(tmp_path / "control.inv")

    target_logger = logging.getLogger(
        "pycsamt.models.modem.control.ModEmControl"
    )
    target_logger.addHandler(caplog.handler)
    target_logger.setLevel(logging.INFO)
    try:
        with caplog.at_level(logging.INFO):
            loaded = ModEmControl.read(path, verbose=1)
    finally:
        target_logger.removeHandler(caplog.handler)

    assert loaded.verbose == 1
    assert any("loaded from" in rec.message for rec in caplog.records)


def test_read_kwargs_forwarded_to_constructor(tmp_path):
    ctrl = ModEmControl.from_config(ModEmConfig())
    path = ctrl.write(tmp_path / "control.inv")
    loaded = ModEmControl.read(path, verbose=2)
    assert loaded.verbose == 2

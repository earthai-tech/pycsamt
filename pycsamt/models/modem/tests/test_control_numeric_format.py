# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regression tests for ModEmControl's/ModEmForwardControl's numeric format.

Fortran's ``G`` edit descriptor on *input* (used by both ``ModEM.inv``
and the forward-control file) requires an explicit decimal point in the
field; without one, the format's own decimal-digit count silently
re-places one, corrupting the value's magnitude (confirmed both by a
real run and by ModEM's own usage-text examples, which always include a
trailing ``.`` even for whole numbers -- see ``UserCtrl.f90``). Python's
``%g`` formatting omits the decimal point for whole numbers (``10.0`` ->
``"10"``) and for many small-magnitude exponentials (``1e-7`` ->
``"1e-07"``), so both writers must not use it for float fields.
"""

from __future__ import annotations

import pytest

from pycsamt.models.modem.config import ModEmConfig
from pycsamt.models.modem.control import _FLOAT_ATTRS, _KEYS, ModEmControl
from pycsamt.models.modem.forward_control import ModEmForwardControl
from pycsamt.models.modem.forward_control import _KEYS as _FWD_KEYS


def _float_values(text: str, keys) -> dict[str, str]:
    """Map each float-attribute's label to its written value string."""
    lines = text.splitlines()
    out: dict[str, str] = {}
    for entry in keys:
        attr, label = entry[0], entry[1]
        is_float = attr in _FLOAT_ATTRS if len(entry) == 2 else not entry[2]
        if not is_float:
            continue
        line = next(l for l in lines if l.startswith(label))
        out[label] = line.rsplit(":", 1)[-1].strip()
    return out


@pytest.mark.parametrize(
    "cfg_kwargs",
    [
        {},  # whole-number defaults: initial_lambda=10.0, lambda_divisor=100.0
        {"initial_lambda": 10.0, "lambda_divisor": 100.0, "initial_alpha": 10.0},
        {"target_rms": 1.0},
        {"rms_diff_tol": 5.0e-4, "lambda_exit": 1.0e-4},
    ],
)
def test_modem_control_write_always_has_decimal_point(tmp_path, cfg_kwargs):
    cfg = ModEmConfig(**cfg_kwargs)
    ctrl = ModEmControl.from_config(cfg)
    path = ctrl.write(tmp_path / "control.inv")
    values = _float_values(path.read_text(), _KEYS)
    missing = {k: v for k, v in values.items() if "." not in v}
    assert missing == {}
    for v in values.values():
        assert float(v) == float(v)  # parses cleanly as a Python float too


def test_modem_control_whole_number_regression(tmp_path):
    # initial_lambda=10.0 -- the exact value that regressed: `.4g` writes
    # "10" (no decimal point), which a real Mod3DMT/Mod2DMT then misreads
    # as if the decimal point were 7 places from the right.
    cfg = ModEmConfig(initial_lambda=10.0)
    ctrl = ModEmControl.from_config(cfg)
    path = ctrl.write(tmp_path / "control.inv")
    line = next(
        l for l in path.read_text().splitlines() if "Initial damping factor" in l
    )
    value = line.rsplit(":", 1)[-1].strip()
    assert "." in value
    assert float(value) == pytest.approx(10.0)


@pytest.mark.parametrize(
    "cfg_kwargs",
    [
        {},  # defaults: tol_em_fwd=1e-7, tol_em_adj=1e-7, tol_divcor=1e-5
        {"tol_em_fwd": 1.0e-7, "tol_em_adj": 1.0e-7, "tol_divcor": 1.0e-5},
    ],
)
def test_modem_forward_control_write_always_has_decimal_point(tmp_path, cfg_kwargs):
    cfg = ModEmConfig(**cfg_kwargs)
    ctrl = ModEmForwardControl.from_config(cfg)
    path = ctrl.write(tmp_path / "fwd.ctrl")
    values = _float_values(path.read_text(), _FWD_KEYS)
    missing = {k: v for k, v in values.items() if "." not in v}
    assert missing == {}


def test_modem_forward_control_tolerance_regression(tmp_path):
    # The exact bug caught live: a written "1e-07" (no decimal point) was
    # read back by a real Mod3DMT as 0.1000000E-13, not 1.0E-7.
    cfg = ModEmConfig(tol_em_fwd=1.0e-7)
    ctrl = ModEmForwardControl.from_config(cfg)
    path = ctrl.write(tmp_path / "fwd.ctrl")
    text = path.read_text()
    line = next(l for l in text.splitlines() if "forward solver" in l)
    value = line.rsplit(":", 1)[-1].strip()
    assert "." in value
    assert float(value) == pytest.approx(1.0e-7)

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Contract tests for :mod:`pycsamt.forward.config3d` (:class:`ForwardConfig3D`).

Mirrors :mod:`pycsamt.forward.tests.test_config2d`: validation branches,
frequency-grid construction, grid assembly for every ``model_type``,
solver kwargs, template file round-trips, and the summary text.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.config3d import ForwardConfig3D
from pycsamt.forward.grid3d import Grid3D

# ─────────────────────────────────────────────────────────────────────────────
# Defaults / validate()
# ─────────────────────────────────────────────────────────────────────────────


def test_default_construction_validates_cleanly():
    cfg = ForwardConfig3D()
    assert cfg.model_type == "halfspace"
    assert cfg.method == "quasi3d"
    cfg.validate()  # must not raise


def test_validate_rejects_non_positive_freq_min():
    cfg = ForwardConfig3D(freq_min=0.0)
    with pytest.raises(ValueError, match="freq_min"):
        cfg.validate()


def test_validate_rejects_freq_max_not_greater_than_freq_min():
    cfg = ForwardConfig3D(freq_min=10.0, freq_max=10.0)
    with pytest.raises(ValueError, match="freq_min"):
        cfg.validate()
    cfg2 = ForwardConfig3D(freq_min=10.0, freq_max=5.0)
    with pytest.raises(ValueError, match="freq_min"):
        cfg2.validate()


def test_validate_rejects_too_few_freqs():
    cfg = ForwardConfig3D(n_freqs=1)
    with pytest.raises(ValueError, match="n_freqs"):
        cfg.validate()


def test_validate_rejects_unknown_method():
    cfg = ForwardConfig3D(method="fullwave3d")
    with pytest.raises(ValueError, match="quasi3d"):
        cfg.validate()


@pytest.mark.parametrize("field", ["nx", "ny", "nz"])
def test_validate_rejects_undersized_grid(field):
    cfg = ForwardConfig3D(**{field: 3})
    with pytest.raises(ValueError, match="nx, ny, nz"):
        cfg.validate()


@pytest.mark.parametrize("field", ["x_max", "y_max", "z_max"])
def test_validate_rejects_non_positive_extents(field):
    cfg = ForwardConfig3D(**{field: 0.0})
    with pytest.raises(ValueError, match="x_max, y_max, z_max"):
        cfg.validate()


def test_validate_rejects_negative_n_pad():
    cfg = ForwardConfig3D(n_pad=-1)
    with pytest.raises(ValueError, match="n_pad"):
        cfg.validate()


def test_validate_rejects_pad_factor_at_or_below_one():
    cfg = ForwardConfig3D(pad_factor=1.0)
    with pytest.raises(ValueError, match="pad_factor"):
        cfg.validate()


def test_validate_rejects_non_positive_bg_rho():
    cfg = ForwardConfig3D(bg_rho=0.0)
    with pytest.raises(ValueError, match="bg_rho"):
        cfg.validate()


def test_validate_rejects_invalid_model_type():
    cfg = ForwardConfig3D(model_type="bogus")
    with pytest.raises(ValueError, match="model_type"):
        cfg.validate()


def test_validate_rejects_non_positive_anomaly_rho_for_block_anomaly():
    cfg = ForwardConfig3D(model_type="block_anomaly", anomaly_rho=0.0)
    with pytest.raises(ValueError, match="anomaly_rho"):
        cfg.validate()


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_validate_rejects_inverted_anomaly_bounds_per_axis(axis):
    kwargs = {
        f"anomaly_{axis}_lo": 5000.0,
        f"anomaly_{axis}_hi": 5000.0,
    }
    cfg = ForwardConfig3D(model_type="block_anomaly", **kwargs)
    with pytest.raises(ValueError, match=f"anomaly_{axis}_lo"):
        cfg.validate()


def test_anomaly_validation_skipped_when_model_type_is_not_block_anomaly():
    cfg = ForwardConfig3D(
        model_type="halfspace",
        anomaly_rho=-1.0,
        anomaly_x_lo=9000.0,
        anomaly_x_hi=1000.0,
        anomaly_y_lo=9000.0,
        anomaly_y_hi=1000.0,
        anomaly_z_lo=9000.0,
        anomaly_z_hi=1000.0,
    )
    cfg.validate()  # must not raise


def test_validate_rejects_zero_stations():
    cfg = ForwardConfig3D(nx_stations=0)
    with pytest.raises(ValueError, match="nx_stations and ny_stations"):
        cfg.validate()
    cfg2 = ForwardConfig3D(ny_stations=0)
    with pytest.raises(ValueError, match="nx_stations and ny_stations"):
        cfg2.validate()


# ─────────────────────────────────────────────────────────────────────────────
# freq_grid()
# ─────────────────────────────────────────────────────────────────────────────


def test_freq_grid_matches_logspace_definition():
    cfg = ForwardConfig3D(freq_min=1e-2, freq_max=1e2, n_freqs=9)
    grid = cfg.freq_grid()
    expected = np.logspace(-2, 2, 9)
    assert grid.shape == (9,)
    assert np.allclose(grid, expected)
    assert np.all(np.diff(grid) > 0)


# ─────────────────────────────────────────────────────────────────────────────
# to_grid()
# ─────────────────────────────────────────────────────────────────────────────


def test_to_grid_halfspace_uses_uniform_bg_rho():
    cfg = ForwardConfig3D(
        model_type="halfspace",
        bg_rho=250.0,
        nx=8,
        ny=8,
        nz=6,
        n_pad=2,
        nx_stations=3,
        ny_stations=3,
    )
    grid = cfg.to_grid()
    assert isinstance(grid, Grid3D)
    assert grid.nx == 8 + 2 * 2
    assert grid.ny == 8 + 2 * 2
    assert grid.nz == 6 + 2
    assert np.all(grid.resistivity == 250.0)
    assert grid.n_stations == 3 * 3


def test_to_grid_block_anomaly_embeds_anomaly_in_background():
    cfg = ForwardConfig3D(
        model_type="block_anomaly",
        bg_rho=100.0,
        anomaly_rho=5.0,
        anomaly_x_lo=2000.0,
        anomaly_x_hi=6000.0,
        anomaly_y_lo=2000.0,
        anomaly_y_hi=6000.0,
        anomaly_z_lo=300.0,
        anomaly_z_hi=1500.0,
    )
    grid = cfg.to_grid()
    assert isinstance(grid, Grid3D)
    assert np.any(grid.resistivity == 100.0)
    assert np.any(grid.resistivity == 5.0)
    assert grid.resistivity.min() == pytest.approx(5.0)


def test_to_grid_random_layered_is_reproducible_with_seed():
    cfg = ForwardConfig3D(model_type="random_layered", nx=8, ny=8, nz=6, n_pad=2)
    grid_a = cfg.to_grid(seed=42)
    grid_b = cfg.to_grid(seed=42)
    grid_c = cfg.to_grid(seed=7)

    assert isinstance(grid_a, Grid3D)
    assert np.array_equal(grid_a.resistivity, grid_b.resistivity)
    assert not np.array_equal(grid_a.resistivity, grid_c.resistivity)
    assert np.all(grid_a.resistivity > 0)


# ─────────────────────────────────────────────────────────────────────────────
# to_solver_kwargs()
# ─────────────────────────────────────────────────────────────────────────────


def test_to_solver_kwargs_shape():
    cfg = ForwardConfig3D(freq_min=1e-2, freq_max=1e2, n_freqs=5, verbose=False)
    kw = cfg.to_solver_kwargs()
    assert set(kw) == {"freqs", "method", "verbose"}
    assert kw["method"] == "quasi3d"
    assert kw["verbose"] is False
    assert np.allclose(kw["freqs"], cfg.freq_grid())


# ─────────────────────────────────────────────────────────────────────────────
# Template / file I/O round-trip
# ─────────────────────────────────────────────────────────────────────────────


def test_template_round_trip_py_format(tmp_path):
    path = tmp_path / "fwd3d.py"
    cfg = ForwardConfig3D(
        freq_min=5e-3,
        freq_max=5e2,
        n_freqs=17,
        nx=24,
        model_type="block_anomaly",
    )
    out = cfg.to_template(path)
    assert out == path
    assert path.exists()
    assert "CONFIG" in path.read_text(encoding="utf-8")

    loaded = ForwardConfig3D.from_file(path)
    assert loaded.freq_min == pytest.approx(5e-3)
    assert loaded.freq_max == pytest.approx(5e2)
    assert loaded.n_freqs == 17
    assert loaded.nx == 24
    assert loaded.model_type == "block_anomaly"
    loaded.validate()


def test_read_alias_matches_from_file(tmp_path):
    path = tmp_path / "fwd3d_alias.py"
    ForwardConfig3D().to_template(path)
    via_read = ForwardConfig3D.read(path)
    via_from_file = ForwardConfig3D.from_file(path)
    assert via_read == via_from_file


def test_write_template_classmethod_uses_defaults(tmp_path):
    path = tmp_path / "fwd3d_defaults.py"
    out = ForwardConfig3D.write_template(path)
    loaded = ForwardConfig3D.from_file(out)
    assert loaded == ForwardConfig3D()


def test_template_round_trip_json_format(tmp_path):
    path = tmp_path / "fwd3d.json"
    cfg = ForwardConfig3D(nx_stations=6, bg_rho=333.0)
    cfg.to_template(path)
    assert path.exists()

    loaded = ForwardConfig3D.from_file(path)
    assert loaded.nx_stations == 6
    assert loaded.bg_rho == pytest.approx(333.0)


def test_from_file_strict_rejects_unknown_key(tmp_path):
    path = tmp_path / "fwd3d_bad.json"
    ForwardConfig3D().to_template(path)
    import json

    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["config"]["not_a_real_field"] = 1
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="Unknown configuration parameter"):
        ForwardConfig3D.from_file(path, strict=True)

    loaded = ForwardConfig3D.from_file(path, strict=False)
    assert not hasattr(loaded, "not_a_real_field")
    loaded.validate()


# ─────────────────────────────────────────────────────────────────────────────
# summary() / __repr__
# ─────────────────────────────────────────────────────────────────────────────


def test_summary_contains_key_fields_for_halfspace():
    cfg = ForwardConfig3D(model_type="halfspace", bg_rho=77.0, nx_stations=4)
    text = cfg.summary()
    assert "ForwardConfig3D" in text
    assert "'halfspace'" in text
    assert "77.0" in text
    assert "quasi3d" in text
    assert "anomaly" not in text


def test_summary_includes_anomaly_line_for_block_anomaly_model_type():
    cfg = ForwardConfig3D(
        model_type="block_anomaly",
        anomaly_rho=3.5,
        anomaly_x_lo=1000.0,
        anomaly_x_hi=2000.0,
        anomaly_y_lo=1000.0,
        anomaly_y_hi=2000.0,
        anomaly_z_lo=100.0,
        anomaly_z_hi=200.0,
    )
    text = cfg.summary()
    assert "anomaly" in text
    assert "3.5" in text
    assert repr(cfg) == text

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Contract tests for :mod:`pycsamt.forward.config2d` (:class:`ForwardConfig2D`).

Covers parameter validation, frequency-grid construction, grid assembly
for every ``model_type`` branch, solver kwargs, template file round-trips,
and the human-readable summary.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.config2d import ForwardConfig2D
from pycsamt.forward.grid2d import Grid2D

# ─────────────────────────────────────────────────────────────────────────────
# Defaults / validate()
# ─────────────────────────────────────────────────────────────────────────────


def test_default_construction_validates_cleanly():
    cfg = ForwardConfig2D()
    assert cfg.model_type == "halfspace"
    cfg.validate()  # must not raise


@pytest.mark.parametrize(
    "field, value",
    [
        ("freq_min", 0.0),
        ("freq_min", -1.0),
    ],
)
def test_validate_rejects_non_positive_freq_min(field, value):
    cfg = ForwardConfig2D(**{field: value})
    with pytest.raises(ValueError, match="freq_min"):
        cfg.validate()


def test_validate_rejects_freq_max_not_greater_than_freq_min():
    cfg = ForwardConfig2D(freq_min=10.0, freq_max=10.0)
    with pytest.raises(ValueError, match="freq_min"):
        cfg.validate()
    cfg2 = ForwardConfig2D(freq_min=10.0, freq_max=5.0)
    with pytest.raises(ValueError, match="freq_min"):
        cfg2.validate()


def test_validate_rejects_too_few_freqs():
    cfg = ForwardConfig2D(n_freqs=1)
    with pytest.raises(ValueError, match="n_freqs"):
        cfg.validate()


@pytest.mark.parametrize("field", ["nx", "nz"])
def test_validate_rejects_undersized_grid(field):
    cfg = ForwardConfig2D(**{field: 3})
    with pytest.raises(ValueError, match="nx and nz"):
        cfg.validate()


@pytest.mark.parametrize("field", ["x_max", "z_max"])
def test_validate_rejects_non_positive_extents(field):
    cfg = ForwardConfig2D(**{field: 0.0})
    with pytest.raises(ValueError, match="x_max and z_max"):
        cfg.validate()
    cfg2 = ForwardConfig2D(**{field: -100.0})
    with pytest.raises(ValueError, match="x_max and z_max"):
        cfg2.validate()


def test_validate_rejects_negative_n_pad():
    cfg = ForwardConfig2D(n_pad=-1)
    with pytest.raises(ValueError, match="n_pad"):
        cfg.validate()


def test_validate_accepts_zero_n_pad():
    cfg = ForwardConfig2D(n_pad=0)
    cfg.validate()  # zero padding is explicitly allowed


def test_validate_rejects_pad_factor_at_or_below_one():
    cfg = ForwardConfig2D(pad_factor=1.0)
    with pytest.raises(ValueError, match="pad_factor"):
        cfg.validate()
    cfg2 = ForwardConfig2D(pad_factor=0.5)
    with pytest.raises(ValueError, match="pad_factor"):
        cfg2.validate()


def test_validate_rejects_non_positive_bg_rho():
    cfg = ForwardConfig2D(bg_rho=0.0)
    with pytest.raises(ValueError, match="bg_rho"):
        cfg.validate()
    cfg2 = ForwardConfig2D(bg_rho=-5.0)
    with pytest.raises(ValueError, match="bg_rho"):
        cfg2.validate()


def test_validate_rejects_invalid_model_type():
    cfg = ForwardConfig2D(model_type="bogus")
    with pytest.raises(ValueError, match="model_type"):
        cfg.validate()


def test_validate_rejects_non_positive_anomaly_rho_only_for_anomaly_type():
    cfg = ForwardConfig2D(model_type="anomaly", anomaly_rho=0.0)
    with pytest.raises(ValueError, match="anomaly_rho"):
        cfg.validate()


def test_validate_rejects_inverted_anomaly_x_bounds():
    cfg = ForwardConfig2D(
        model_type="anomaly", anomaly_x_lo=5000.0, anomaly_x_hi=5000.0
    )
    with pytest.raises(ValueError, match="anomaly_x_lo"):
        cfg.validate()
    cfg2 = ForwardConfig2D(
        model_type="anomaly", anomaly_x_lo=6000.0, anomaly_x_hi=2000.0
    )
    with pytest.raises(ValueError, match="anomaly_x_lo"):
        cfg2.validate()


def test_validate_rejects_inverted_anomaly_z_bounds():
    cfg = ForwardConfig2D(model_type="anomaly", anomaly_z_lo=1500.0, anomaly_z_hi=300.0)
    with pytest.raises(ValueError, match="anomaly_z_lo"):
        cfg.validate()


def test_anomaly_validation_skipped_when_model_type_is_not_anomaly():
    # These anomaly values are individually invalid, but since model_type
    # is 'halfspace' they must be ignored entirely by validate().
    cfg = ForwardConfig2D(
        model_type="halfspace",
        anomaly_rho=-1.0,
        anomaly_x_lo=9000.0,
        anomaly_x_hi=1000.0,
        anomaly_z_lo=9000.0,
        anomaly_z_hi=1000.0,
    )
    cfg.validate()  # must not raise


def test_validate_rejects_zero_stations():
    cfg = ForwardConfig2D(n_stations=0)
    with pytest.raises(ValueError, match="n_stations"):
        cfg.validate()


def test_validate_accepts_single_station():
    cfg = ForwardConfig2D(n_stations=1)
    cfg.validate()


# ─────────────────────────────────────────────────────────────────────────────
# freq_grid()
# ─────────────────────────────────────────────────────────────────────────────


def test_freq_grid_matches_logspace_definition():
    cfg = ForwardConfig2D(freq_min=1e-2, freq_max=1e2, n_freqs=9)
    grid = cfg.freq_grid()
    expected = np.logspace(-2, 2, 9)
    assert grid.shape == (9,)
    assert np.allclose(grid, expected)
    assert grid[0] == pytest.approx(1e-2)
    assert grid[-1] == pytest.approx(1e2)
    # strictly increasing
    assert np.all(np.diff(grid) > 0)


# ─────────────────────────────────────────────────────────────────────────────
# to_grid()
# ─────────────────────────────────────────────────────────────────────────────


def test_to_grid_halfspace_uses_uniform_bg_rho():
    cfg = ForwardConfig2D(
        model_type="halfspace",
        bg_rho=250.0,
        nx=10,
        nz=8,
        n_pad=3,
        n_stations=4,
    )
    grid = cfg.to_grid()
    assert isinstance(grid, Grid2D)
    # core region matches requested cell counts, padding is added on top
    assert grid.nx == 10 + 2 * 3
    assert grid.nz == 8 + 3
    assert np.all(grid.resistivity == 250.0)
    assert grid.n_stations == 4


def test_to_grid_anomaly_embeds_anomaly_block_in_background():
    cfg = ForwardConfig2D(
        model_type="anomaly",
        bg_rho=100.0,
        anomaly_rho=5.0,
        anomaly_x_lo=2000.0,
        anomaly_x_hi=6000.0,
        anomaly_z_lo=300.0,
        anomaly_z_hi=1500.0,
        nx=40,
        nz=30,
        x_max=10_000.0,
        z_max=6_000.0,
    )
    grid = cfg.to_grid()
    assert isinstance(grid, Grid2D)
    # background value must still be present somewhere
    assert np.any(grid.resistivity == 100.0)
    # anomaly value must appear inside the grid (embedded block)
    assert np.any(grid.resistivity == 5.0)
    # anomaly resistivity strictly lower than background everywhere it appears
    assert grid.resistivity.min() == pytest.approx(5.0)


def test_to_grid_random_builds_grid_without_station_x_max_kwarg():
    """``Grid2D.random()`` has no ``station_x_max`` parameter.

    ``to_grid()`` must not forward it for the 'random' branch (it is only
    meaningful for 'halfspace'/'anomaly'). Regression test for a bug where
    the shared ``kw`` dict was passed through unfiltered and raised
    ``TypeError: random() got an unexpected keyword argument
    'station_x_max'`` for every ``model_type="random"`` config.
    """
    cfg = ForwardConfig2D(
        model_type="random", nx=12, nz=10, n_pad=2, station_x_max=5000.0
    )
    grid = cfg.to_grid(seed=42)
    assert isinstance(grid, Grid2D)
    assert grid.nx == 12 + 2 * 2
    assert grid.nz == 10 + 2

    grid_same_seed = cfg.to_grid(seed=42)
    assert np.array_equal(grid.resistivity, grid_same_seed.resistivity)
    grid_other_seed = cfg.to_grid(seed=43)
    assert not np.array_equal(grid.resistivity, grid_other_seed.resistivity)


# ─────────────────────────────────────────────────────────────────────────────
# to_solver_kwargs()
# ─────────────────────────────────────────────────────────────────────────────


def test_to_solver_kwargs_shape():
    cfg = ForwardConfig2D(freq_min=1e-2, freq_max=1e2, n_freqs=5, verbose=False)
    kw = cfg.to_solver_kwargs()
    assert set(kw) == {"freqs", "verbose"}
    assert kw["verbose"] is False
    assert np.allclose(kw["freqs"], cfg.freq_grid())


# ─────────────────────────────────────────────────────────────────────────────
# Template / file I/O round-trip
# ─────────────────────────────────────────────────────────────────────────────


def test_template_round_trip_py_format(tmp_path):
    path = tmp_path / "fwd2d.py"
    cfg = ForwardConfig2D(
        freq_min=5e-3, freq_max=5e2, n_freqs=17, nx=44, model_type="anomaly"
    )
    out = cfg.to_template(path)
    assert out == path
    assert path.exists()
    assert "CONFIG" in path.read_text()

    loaded = ForwardConfig2D.from_file(path)
    assert loaded.freq_min == pytest.approx(5e-3)
    assert loaded.freq_max == pytest.approx(5e2)
    assert loaded.n_freqs == 17
    assert loaded.nx == 44
    assert loaded.model_type == "anomaly"
    loaded.validate()


def test_read_alias_matches_from_file(tmp_path):
    path = tmp_path / "fwd2d_alias.py"
    ForwardConfig2D().to_template(path)
    via_read = ForwardConfig2D.read(path)
    via_from_file = ForwardConfig2D.from_file(path)
    assert via_read == via_from_file


def test_write_template_classmethod_uses_defaults(tmp_path):
    path = tmp_path / "fwd2d_defaults.py"
    out = ForwardConfig2D.write_template(path)
    loaded = ForwardConfig2D.from_file(out)
    assert loaded == ForwardConfig2D()


@pytest.mark.parametrize("suffix,fmt", [(".json", None), (".yml", None)])
def test_template_round_trip_json_and_yaml(tmp_path, suffix, fmt):
    path = tmp_path / f"fwd2d{suffix}"
    cfg = ForwardConfig2D(n_stations=6, bg_rho=333.0)
    cfg.to_template(path, fmt=fmt)
    assert path.exists()

    loaded = ForwardConfig2D.from_file(path)
    assert loaded.n_stations == 6
    assert loaded.bg_rho == pytest.approx(333.0)


def test_from_file_strict_rejects_unknown_key(tmp_path):
    path = tmp_path / "fwd2d_bad.py"
    ForwardConfig2D().to_template(path)
    text = path.read_text()
    # inject a bogus key into CONFIG dict
    text = text.replace("CONFIG = {", "CONFIG = {\n    'not_a_real_field': 1,")
    path.write_text(text)

    with pytest.raises(ValueError, match="Unknown configuration parameter"):
        ForwardConfig2D.from_file(path, strict=True)

    # non-strict mode silently drops the unknown key
    loaded = ForwardConfig2D.from_file(path, strict=False)
    assert not hasattr(loaded, "not_a_real_field")
    loaded.validate()


# ─────────────────────────────────────────────────────────────────────────────
# summary() / __repr__
# ─────────────────────────────────────────────────────────────────────────────


def test_summary_contains_key_fields_for_halfspace():
    cfg = ForwardConfig2D(model_type="halfspace", bg_rho=77.0, n_stations=9)
    text = cfg.summary()
    assert "ForwardConfig2D" in text
    assert "'halfspace'" in text
    assert "77.0" in text
    assert "9" in text
    # halfspace models must not print the anomaly-specific line
    assert "anomaly" not in text


def test_summary_includes_anomaly_line_for_anomaly_model_type():
    cfg = ForwardConfig2D(
        model_type="anomaly",
        anomaly_rho=3.5,
        anomaly_x_lo=1000.0,
        anomaly_x_hi=2000.0,
        anomaly_z_lo=100.0,
        anomaly_z_hi=200.0,
    )
    text = cfg.summary()
    assert "anomaly" in text
    assert "3.5" in text
    assert "1000.0" in text and "2000.0" in text
    assert repr(cfg) == text

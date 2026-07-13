"""Tests for OccamStartup.read() and OccamIter.read()."""

import math
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"
STARTUP_FILE = DATA_DIR / "Startup"
ITER_FILE = DATA_DIR / "ITER17.iter"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(), reason="sample occam2D data not available"
)


# -----------------------------------------------------------------------
# OccamStartup.read()
# -----------------------------------------------------------------------
@pytest.fixture(scope="module")
def startup():
    from pycsamt.models.occam2d.startup import OccamStartup

    return OccamStartup.read(STARTUP_FILE)


def test_startup_format(startup):
    assert startup.format_str == "OCCAMITER_FLEX"


def test_startup_iteration_zero(startup):
    assert startup.iteration == 0


def test_startup_model_file(startup):
    assert startup.model_file == "Occam2DModel"


def test_startup_data_file(startup):
    assert startup.data_file == "OccamDataFile.dat"


def test_startup_max_iterations(startup):
    assert startup.max_iterations == 100


def test_startup_target_misfit(startup):
    assert math.isclose(startup.target_misfit, 1.0, rel_tol=1e-5)


def test_startup_roughness_type(startup):
    assert startup.roughness_type == 1


def test_startup_n_params(startup):
    assert startup.n_params == 3752


def test_startup_param_values_length(startup):
    assert len(startup.param_values) == 3752


def test_startup_param_values_uniform(startup):
    # Startup uses a uniform halfspace: all log10(rho) values should be 2.0
    assert np.allclose(startup.param_values, 2.0, atol=1e-4)


def test_startup_misfit_not_reached(startup):
    assert startup.misfit_reached is False


def test_startup_lagrange_value(startup):
    assert math.isclose(startup.lagrange_value, 5.0, rel_tol=1e-4)


def test_startup_stepsize_cut_count(startup):
    assert startup.stepsize_cut_count == 8


def test_startup_debug_level(startup):
    assert startup.debug_level == 1


# -----------------------------------------------------------------------
# OccamIter.read()
# -----------------------------------------------------------------------
@pytest.fixture(scope="module")
def iter17():
    from pycsamt.models.occam2d.startup import OccamIter

    return OccamIter.read(ITER_FILE)


def test_iter_format(iter17):
    assert iter17.format_str == "OCCAMITER_FLEX"


def test_iter_iteration_number(iter17):
    assert iter17.iteration == 17


def test_iter_misfit_value(iter17):
    assert math.isclose(iter17.misfit_value, 0.9977012, rel_tol=1e-4)


def test_iter_misfit_reached(iter17):
    assert iter17.misfit_reached is True


def test_iter_lagrange_value(iter17):
    assert math.isclose(iter17.lagrange_value, 0.5204059, rel_tol=1e-4)


def test_iter_roughness_value(iter17):
    assert math.isclose(iter17.roughness_value, 240.7858, rel_tol=1e-4)


def test_iter_n_params(iter17):
    assert iter17.n_params == 3752


def test_iter_param_values_length(iter17):
    assert len(iter17.param_values) == 3752


def test_iter_param_values_not_uniform(iter17):
    # Should not be all 2.0 (unlike Startup)
    assert not np.allclose(iter17.param_values, 2.0, atol=0.5)


def test_iter_first_param(iter17):
    # First value from ITER17: 0.9180370
    assert math.isclose(iter17.param_values[0], 0.9180370, rel_tol=1e-4)


def test_iter_to_resistivity_shape(iter17):
    rho = iter17.to_resistivity()
    assert rho.shape == (3752,)


def test_iter_to_resistivity_values(iter17):
    rho = iter17.to_resistivity()
    # log10(rho[0]) = 0.9180370 → rho ≈ 8.28 Ω·m
    assert math.isclose(rho[0], 10**0.9180370, rel_tol=1e-4)
    # All resistivities should be positive
    assert np.all(rho > 0)


def test_iter_log10_rho_stats(iter17):
    stats = iter17.log10_rho_stats
    assert "min" in stats and "max" in stats
    assert stats["min"] < stats["max"]
    # log10(rho) can be negative (resistivity < 1 Ω·m is physically valid)
    assert stats["max"] > stats["min"]


def test_iter_model_file(iter17):
    assert iter17.model_file == "Occam2DModel"


def test_iter_data_file(iter17):
    assert iter17.data_file == "OccamDataFile.dat"


# -----------------------------------------------------------------------
# Cross-read guards
# -----------------------------------------------------------------------
def test_startup_read_raises_on_iter_file():
    """OccamStartup.read() should reject a non-zero iteration file."""
    from pycsamt.models.occam2d.startup import OccamStartup

    with pytest.raises(ValueError, match="Iteration: 0"):
        OccamStartup.read(ITER_FILE)


def test_iter_read_raises_on_startup_file():
    """OccamIter.read() should reject a Startup (iteration == 0) file."""
    from pycsamt.models.occam2d.startup import OccamIter

    with pytest.raises(ValueError, match="Startup"):
        OccamIter.read(STARTUP_FILE)


# -----------------------------------------------------------------------
# Error handling
# -----------------------------------------------------------------------
def test_missing_file_raises_startup():
    from pycsamt.models.occam2d.startup import OccamStartup

    with pytest.raises(FileNotFoundError):
        OccamStartup.read("/nonexistent/Startup")


def test_missing_file_raises_iter():
    from pycsamt.models.occam2d.startup import OccamIter

    with pytest.raises(FileNotFoundError):
        OccamIter.read("/nonexistent/ITER1.iter")


def test_invalid_format_raises(tmp_path):
    from pycsamt.models.occam2d.startup import OccamIter

    bad = tmp_path / "bad.iter"
    bad.write_text("Format: SOMETHING_ELSE\nIteration: 1\n")
    with pytest.raises(ValueError, match="Expected format"):
        OccamIter.read(bad)

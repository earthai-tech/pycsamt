"""Tests for OccamLog.read()."""

import math
from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"
LOG_FILE = DATA_DIR / "LogFile.logfile"

pytestmark = pytest.mark.skipif(
    not LOG_FILE.exists(), reason="sample log file not available"
)


@pytest.fixture(scope="module")
def log():
    from pycsamt.models.occam2d.log import OccamLog

    return OccamLog.read(LOG_FILE)


# ------------------------------------------------------------------
# Shape / size
# ------------------------------------------------------------------
def test_n_iter(log):
    assert log.n_iter == 17


def test_all_arrays_same_length(log):
    n = log.n_iter
    assert len(log.rms) == n
    assert len(log.roughness) == n
    assert len(log.stepsize) == n
    assert len(log.lagrange) == n


# ------------------------------------------------------------------
# Iteration numbers
# ------------------------------------------------------------------
def test_iteration_numbers_monotonic(log):
    assert list(log.iterations) == list(range(1, 18))


def test_first_iteration(log):
    assert log.iterations[0] == 1


def test_last_iteration(log):
    assert log.iterations[-1] == 17


# ------------------------------------------------------------------
# RMS values (spot-checks from the sample file)
# ------------------------------------------------------------------
def test_rms_iter1(log):
    # Iteration 1 achieves "AND IS = 1.452771" → ROUGHNESS IS = 4.610795
    assert math.isclose(log.rms[0], 1.452771, rel_tol=1e-4)


def test_rms_iter14(log):
    # Last "AND IS" before ROUGHNESS for iter 14 is 0.9996567
    assert math.isclose(log.rms[13], 0.9996567, rel_tol=1e-4)


def test_rms_iter16(log):
    # Iter 16 final "AND IS" = 0.9977012
    assert math.isclose(log.rms[15], 0.9977012, rel_tol=1e-4)


def test_rms_decreasing_overall(log):
    # First finite RMS should be higher than the minimum
    finite = log.rms[np.isfinite(log.rms)]
    assert finite[0] > finite.min()


# ------------------------------------------------------------------
# Roughness and stepsize (spot-checks)
# ------------------------------------------------------------------
def test_roughness_iter1(log):
    assert math.isclose(log.roughness[0], 4.610795, rel_tol=1e-4)


def test_stepsize_iter1(log):
    assert math.isclose(log.stepsize[0], 0.5271110, rel_tol=1e-4)


def test_last_iter_roughness_nan(log):
    # Iteration 17 ended with convergence problems — no ROUGHNESS IS written
    assert math.isnan(log.roughness[-1])


def test_last_iter_stepsize_nan(log):
    assert math.isnan(log.stepsize[-1])


# ------------------------------------------------------------------
# Lagrange
# ------------------------------------------------------------------
def test_lagrange_iter1(log):
    # "MINIMUM TOL FROM fminocc IS AT MU = 2.381966"
    assert math.isclose(log.lagrange[0], 2.381966, rel_tol=1e-4)


# ------------------------------------------------------------------
# Derived properties
# ------------------------------------------------------------------
def test_converged(log):
    # Min RMS across all iters is < 1.0, so converged = True
    assert log.converged is True


def test_best_iteration(log):
    # Iter 16 achieves rms=0.9977012, the minimum of all finite values
    assert log.best_iteration == 16


def test_summary_string(log):
    s = log.summary()
    assert "17 iterations" in s
    assert "converged: True" in s


# ------------------------------------------------------------------
# Error handling
# ------------------------------------------------------------------
def test_missing_file_raises():
    from pycsamt.models.occam2d.log import OccamLog

    with pytest.raises(FileNotFoundError):
        OccamLog.read("/nonexistent/path/log.logfile")


def test_empty_log(tmp_path):
    """An empty file should parse to zero iterations without crashing."""
    from pycsamt.models.occam2d.log import OccamLog

    f = tmp_path / "empty.logfile"
    f.write_text("")
    log = OccamLog.read(f)
    assert log.n_iter == 0
    assert log.converged is False

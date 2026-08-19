"""Focused tests for robust Occam1D convergence-log handling."""

import math

import numpy as np
import pytest

from pycsamt.models.occam1d import Occam1DLog, Occam1DLogRecord


def test_record_is_immutable_and_preserves_missing_values():
    record = Occam1DLogRecord(iteration=3, rms=1.2)

    assert math.isnan(record.roughness)
    assert not record.is_complete
    with pytest.raises(AttributeError):
        record.rms = 0.9


@pytest.mark.parametrize("field", ["rms", "roughness", "lagrange", "stepsize"])
def test_record_rejects_negative_and_infinite_metrics(field):
    with pytest.raises(ValueError):
        Occam1DLogRecord(iteration=1, **{field: -1.0})
    with pytest.raises(ValueError):
        Occam1DLogRecord(iteration=1, **{field: np.inf})


def test_read_normalizes_dialects_and_fortran_exponents(tmp_path):
    path = tmp_path / "occam.log"
    path.write_text(
        "ITERATION 1\n"
        "R.M.S. = 2.5D+00\n"
        "ROUGHNESS IS: 1.2D+02\n"
        "MU = 4.0D-01\n"
        "STEP SIZE = 5.0D-01\n"
        "Iteration: 2 RMS: 9.5E-01 ROUGHNESS: 130\n",
        encoding="utf8",
    )

    log = Occam1DLog.read(path)

    assert log.n_iterations == 2
    assert log.best_iteration == 2
    assert log.has_converged
    assert log.records[0].is_complete
    assert log.records[0].lagrange == pytest.approx(0.4)
    assert log.records[0].stepsize == pytest.approx(0.5)


def test_read_strict_and_permissive_empty_log(tmp_path):
    path = tmp_path / "empty.log"
    path.write_text("solver banner only\n", encoding="utf8")

    with pytest.raises(ValueError, match="No Occam iteration"):
        Occam1DLog.read(path)

    log = Occam1DLog.read(path, strict=False)
    assert log.n_iterations == 0
    assert log.path == path.resolve()
    assert log.warnings


def test_restarted_and_incomplete_histories_are_diagnosed(tmp_path):
    path = tmp_path / "restart.log"
    path.write_text(
        "ITERATION 1\nRMS = 2\nROUGHNESS = 10\n"
        "ITERATION 2\nRMS = 1.5\n"
        "ITERATION 1\nRMS = 1.1\nROUGHNESS = 12\n",
        encoding="utf8",
    )

    log = Occam1DLog.read(path)

    assert len(log.warnings) == 2
    assert log.get_iteration(1).rms == pytest.approx(1.1)
    assert len(log.get_iteration(1, all_matches=True)) == 2


def test_queries_and_exports_remain_aligned(tmp_path):
    log = Occam1DLog(
        [
            {
                "iteration": 1,
                "rms": 3.0,
                "roughness": 10.0,
                "lagrange": 0.5,
                "stepsize": 1.0,
            },
            {"iteration": 2, "rms": 0.8},
        ],
        target_misfit=0.9,
    )

    assert log.initial_rms == 3.0
    assert log.final_rms == 0.8
    assert log.best_index == 1
    assert log.converged(1.0)
    assert log.to_array().shape == (2, 5)
    assert log.diagnostics()["has_converged"] is True
    assert "best_iteration=2" in log.summary()

    complete = log.to_csv(tmp_path / "complete.csv", include_missing=False)
    assert len(complete.read_text(encoding="utf8").splitlines()) == 2


def test_invalid_public_query_arguments_are_rejected():
    log = Occam1DLog([{"iteration": 1, "rms": 1, "roughness": 1}])

    with pytest.raises((TypeError, ValueError)):
        log.converged(0)
    with pytest.raises((TypeError, ValueError)):
        log.get_iteration(-1)
    with pytest.raises(KeyError):
        log.get_iteration(99)

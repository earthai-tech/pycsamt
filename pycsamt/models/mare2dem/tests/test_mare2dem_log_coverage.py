# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for models/mare2dem/log.py.

Uses small hand-written ``OccamLog.2012.0``-style logs and group-RMS CSV
logs so the test does not depend on the (gitignored, not present in CI)
bundled ``data/mare2dem`` example directory.
"""

from __future__ import annotations

from pycsamt.models.mare2dem.log import (
    GroupRMSLog,
    IterationRecord,
    Mare2DEMLog,
    read_group_rms_log,
)

_TWO_ITER_LOG = """\
Format:     OccamLog.2012.0
Start time: 12/21/2022 15:36:27.271

** Iteration     1 **

 Occam iteration completed successfully:

                 Target Misfit:        1.000
                  Model Misfit:        6.875
                     Roughness:        1.407
                    Optimal Mu:        4.000
            Convergence Status:        0

** Iteration     2 **

 Occam iteration completed successfully:

                 Target Misfit:        1.000
                  Model Misfit:        0.998
                     Roughness:        37.803
                    Optimal Mu:        3.000
            Convergence Status:        1

 Model roughness can not be decreased further, stopping.
 End time: 12/21/2022 15:40:12.284
"""

_TEXT_CONVERGED_LOG = """\
** Iteration     1 **
            Model Misfit:        1.234
               Roughness:        2.000
              Optimal Mu:        1.000
 Target misfit achieved, stopping inversion.
"""

_MISFIT_REACHED_LOG = """\
** Iteration     1 **
            Model Misfit:        0.999
 Misfit reached target, stopping.
"""

_MISSING_ROUGHNESS_MU_LOG = """\
** Iteration     1 **
            Model Misfit:        5.0
** Iteration     2 **
            Model Misfit:        3.0
               Roughness:        2.0
"""


def test_mare2dem_log_parses_two_iterations_and_converges(tmp_path):
    path = tmp_path / "demo.logfile"
    path.write_text(_TWO_ITER_LOG)

    log = Mare2DEMLog(path)

    assert log.n_iterations == 2
    assert log.converged is True
    assert log.final_rms == 0.998
    assert log.rms_history() == [6.875, 0.998]
    rec = log.iterations[0]
    assert isinstance(rec, IterationRecord)
    assert rec.iteration == 1
    assert rec.rms == 6.875
    assert rec.roughness == 1.407
    assert rec.lambda_ == 4.000
    assert "n_iterations=2" in repr(log)
    assert "converged=True" in repr(log)


def test_mare2dem_log_missing_file_leaves_empty(tmp_path):
    log = Mare2DEMLog(tmp_path / "does_not_exist.logfile")
    assert log.n_iterations == 0
    assert log.converged is False
    assert log.final_rms is None
    assert log.rms_history() == []


def test_mare2dem_log_target_misfit_achieved_text_sets_converged(tmp_path):
    path = tmp_path / "target.logfile"
    path.write_text(_TEXT_CONVERGED_LOG)
    log = Mare2DEMLog(path)
    assert log.converged is True
    assert log.n_iterations == 1


def test_mare2dem_log_misfit_reached_text_sets_converged(tmp_path):
    path = tmp_path / "misfit_reached.logfile"
    path.write_text(_MISFIT_REACHED_LOG)
    log = Mare2DEMLog(path)
    assert log.converged is True


def test_mare2dem_log_missing_roughness_and_mu_default_to_zero(tmp_path):
    path = tmp_path / "partial.logfile"
    path.write_text(_MISSING_ROUGHNESS_MU_LOG)
    log = Mare2DEMLog(path)

    assert log.n_iterations == 2
    first, second = log.iterations
    assert first.roughness == 0.0
    assert first.lambda_ == 0.0
    assert second.roughness == 2.0
    assert second.lambda_ == 0.0


def test_mare2dem_log_resolves_relative_path(tmp_path, monkeypatch):
    path = tmp_path / "rel.logfile"
    path.write_text(_TWO_ITER_LOG)
    monkeypatch.chdir(tmp_path)
    log = Mare2DEMLog("rel.logfile")
    assert log.n_iterations == 2
    assert log.path.is_absolute()


# ---------------------------------------------------------------------------
# GroupRMSLog / read_group_rms_log
# ---------------------------------------------------------------------------

_GROUP_RMS_TEXT = """\
  Iteration,                Total RMS,                    CSEM,                      MT
           1                   29.738                   48.837                   12.235
           2                   22.610                   34.567                   13.287
           3                    1.001                    0.910                    1.043
"""


def test_read_group_rms_log(tmp_path):
    path = tmp_path / "demo.group_rms.log"
    path.write_text(_GROUP_RMS_TEXT)

    log = read_group_rms_log(path)

    assert isinstance(log, GroupRMSLog)
    assert log.headers == ["Iteration", "Total RMS", "CSEM", "MT"]
    assert log.n_groups == 4
    assert log.n_iterations == 3
    assert log.rms_log.shape == (3, 4)
    assert log.rms_log[0, 0] == 1.0
    assert log.rms_log[-1, 1] == 1.001
    assert "n_iterations=3" in repr(log)


def test_read_group_rms_log_missing_file_raises(tmp_path):
    import pytest

    with pytest.raises(FileNotFoundError):
        read_group_rms_log(tmp_path / "missing.log")


def test_read_group_rms_log_empty_file(tmp_path):
    path = tmp_path / "empty.log"
    path.write_text("")
    log = read_group_rms_log(path)
    assert log.headers == []
    assert log.n_iterations == 0


def test_read_group_rms_log_header_only_no_data_rows(tmp_path):
    path = tmp_path / "header_only.log"
    path.write_text("Iteration, Total RMS\n")
    log = read_group_rms_log(path)
    assert log.headers == ["Iteration", "Total RMS"]
    assert log.n_iterations == 0
    assert log.rms_log.shape == (0, 0)


def test_read_group_rms_log_partial_trailing_row_is_dropped(tmp_path):
    # Header names 2 columns but only a single trailing numeric value is
    # present -> n_rows == 0, so no partial row is materialized.
    path = tmp_path / "partial.log"
    path.write_text("Iteration, Total RMS\n1\n")
    log = read_group_rms_log(path)
    assert log.headers == ["Iteration", "Total RMS"]
    assert log.n_iterations == 0
    assert log.rms_log.shape == (0, 0)


def test_read_group_rms_log_reexported_from_log_module():
    # Mare2DEMLog re-exports GroupRMSLog/read_group_rms_log from iotools.
    from pycsamt.models.mare2dem.iotools.group_rms import (
        GroupRMSLog as _Direct,
        read_group_rms_log as _direct_reader,
    )

    assert GroupRMSLog is _Direct
    assert read_group_rms_log is _direct_reader

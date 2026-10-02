# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Self-contained coverage tests for ``InversionResult``.

Unlike ``test_mare2dem_results.py`` (which relies on the gitignored
``data/mare2dem/`` example directory and is skipped in CI), these
tests build minimal synthetic output files directly in ``tmp_path``
so they contribute real coverage in every environment.
"""

from __future__ import annotations

from pathlib import Path

from pycsamt.models.mare2dem import InversionResult
from pycsamt.models.mare2dem.iotools.emdata import EMDataFile, write_emdata
from pycsamt.models.mare2dem.iotools.resistivity import (
    ResistivityFile,
    write_resistivity,
)

_LOG_TEXT = """\
Format:     OccamLog.2012.0
** Iteration     1 **
                  Model Misfit:        6.875
                     Roughness:        1.407
                    Optimal Mu:        4.000
            Convergence Status:        0
** Iteration     2 **
                  Model Misfit:        0.987
                     Roughness:        1.900
                    Optimal Mu:        2.000
            Convergence Status:        1
"""


def _write_synthetic_run(workdir: Path) -> None:
    (workdir / "run.logfile").write_text(_LOG_TEXT)
    write_resistivity(ResistivityFile(), workdir / "run.resistivity")
    write_emdata(EMDataFile(), workdir / "run.emdata")
    write_emdata(EMDataFile(), workdir / "run_MARE2DEM.emdata")


def test_scan_missing_workdir_verbose_logs_warning(tmp_path):
    result = InversionResult(tmp_path / "does_not_exist", verbose=1)
    assert result.log is None
    assert result.model is None
    assert result.data is None
    assert result.response is None
    assert result.converged is False
    assert result.final_rms is None
    assert result.n_iterations == 0


def test_scan_populates_all_outputs_and_properties(tmp_path):
    _write_synthetic_run(tmp_path)

    result = InversionResult(tmp_path, verbose=1)

    assert result.log is not None
    assert result.model is not None
    assert result.data is not None
    assert result.response is not None

    assert result.converged is True
    assert result.final_rms == 0.987
    assert result.n_iterations == 2


def test_summary_and_print_summary(tmp_path, capsys):
    _write_synthetic_run(tmp_path)
    result = InversionResult(tmp_path)

    summary = result.summary()
    assert str(result.workdir) in summary
    assert "converged   : True" in summary
    assert "final RMS   : 0.987" in summary
    assert "n_iterations: 2" in summary
    assert "model" in summary
    assert "data" in summary
    assert "response" in summary

    result.print_summary()
    captured = capsys.readouterr()
    assert captured.out.strip() == summary


def test_repr_reports_convergence_and_rms(tmp_path):
    _write_synthetic_run(tmp_path)
    result = InversionResult(tmp_path)
    r = repr(result)
    assert "converged=True" in r
    assert "final_rms=0.987" in r


def test_empty_dir_leaves_all_outputs_none(tmp_path):
    result = InversionResult(tmp_path)
    assert result.converged is False
    assert result.final_rms is None
    assert result.n_iterations == 0
    assert "converged   : False" in result.summary()

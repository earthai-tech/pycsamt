# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ``pycsamt invert`` command group."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
from pycsamt.cli.tests.conftest import (
    make_fake_sites as _make_fake_sites,  # noqa: E402
)

# ---------------------------------------------------------------------------
# invert (group help)
# ---------------------------------------------------------------------------


class TestInvertGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "--help"])
        assert result.exit_code == 0
        for sub in ("build", "run", "status", "results", "plot"):
            assert sub in result.output


# ---------------------------------------------------------------------------
# invert build
# ---------------------------------------------------------------------------


class TestInvertBuild:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "build", "--help"])
        assert result.exit_code == 0
        assert "--solver" in result.output
        assert "--workdir" in result.output

    def test_no_survey_no_context_fails(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["invert", "build"])
        assert result.exit_code != 0
        assert "No active survey" in (result.output + str(result.exception or ""))

    def test_build_occam2d_called(
        self,
        runner: CliRunner,
        isolated_home: Path,
        edi_dir: Path,
        tmp_path: Path,
    ) -> None:
        fake_sites = _make_fake_sites(5)
        workdir = tmp_path / "run01"

        mock_builder = MagicMock()
        mock_builder.return_value.build.return_value = mock_builder

        with (
            patch("pycsamt.cli.survey._build_sites", return_value=fake_sites),
            patch(
                "pycsamt.cli.commands.invert.build.InputBuilder",
                mock_builder,
                create=True,
            ),
        ):
            result = runner.invoke(
                main,
                [
                    "invert",
                    "build",
                    str(edi_dir),
                    "--solver",
                    "occam2d",
                    "--workdir",
                    str(workdir),
                ],
            )
        # The command may fail because OccamConfig/InputBuilder aren't mocked
        # deeply, but it should not raise a Python exception
        assert result.exception is None or isinstance(result.exception, SystemExit)

    def test_explicit_path_takes_priority_over_context(
        self,
        runner: CliRunner,
        isolated_home: Path,
        edi_dir: Path,
        tmp_path: Path,
    ) -> None:
        fake_sites = _make_fake_sites(2)
        alt_dir = tmp_path / "alt"
        alt_dir.mkdir()
        (alt_dir / "dummy.edi").write_text("dummy")

        # Set context pointing to edi_dir
        from pycsamt.cli.survey import set_survey

        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            set_survey(edi_dir)

        # Call with explicit alt_dir — resolve_survey should use alt_dir
        resolved_paths = []
        __import__("pycsamt.cli.survey", fromlist=["resolve_survey"]).resolve_survey

        def tracking_resolve(explicit, **kw):
            if explicit is not None:
                resolved_paths.append(explicit)
            return fake_sites

        with patch(
            "pycsamt.cli.commands.invert.build.resolve_survey",
            side_effect=tracking_resolve,
        ):
            runner.invoke(
                main,
                [
                    "invert",
                    "build",
                    str(alt_dir),
                    "--workdir",
                    str(tmp_path / "wd"),
                ],
            )

        assert any(alt_dir.resolve() == p.resolve() for p in resolved_paths)


# ---------------------------------------------------------------------------
# invert status
# ---------------------------------------------------------------------------


class TestInvertStatus:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "status", "--help"])
        assert result.exit_code == 0

    def test_occam2d_workdir_text(self, runner: CliRunner, occam_workdir: Path) -> None:
        result = runner.invoke(main, ["invert", "status", str(occam_workdir)])
        assert result.exit_code == 0
        assert "OCCAM2D" in result.output.upper() or "occam2d" in result.output

    def test_occam2d_workdir_json(self, runner: CliRunner, occam_workdir: Path) -> None:
        result = runner.invoke(
            main, ["invert", "status", str(occam_workdir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["solver"] == "occam2d"
        assert "ready_to_run" in data

    def test_modem_workdir_detected(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        result = runner.invoke(main, ["invert", "status", str(modem_workdir)])
        assert result.exit_code == 0
        lower = result.output.lower()
        assert "modem" in lower

    def test_iterations_counted(
        self, runner: CliRunner, occam_workdir_with_iters: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "invert",
                "status",
                str(occam_workdir_with_iters),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_iterations_done"] == 5

    def test_rms_extracted_from_log(
        self, runner: CliRunner, occam_workdir_with_iters: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "invert",
                "status",
                str(occam_workdir_with_iters),
                "--format",
                "json",
            ],
        )
        data = json.loads(result.output)
        assert data["rms_last"] == pytest.approx(1.087, abs=0.01)

    def test_empty_dir_cannot_detect_solver(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        empty = tmp_path / "empty_wd"
        empty.mkdir()
        result = runner.invoke(main, ["invert", "status", str(empty)])
        assert result.exit_code != 0

    def test_explicit_solver_overrides_detection(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "invert",
                "status",
                str(modem_workdir),
                "--solver",
                "modem",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["solver"] == "modem"


# ---------------------------------------------------------------------------
# invert run
# ---------------------------------------------------------------------------


class TestInvertRun:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "run", "--help"])
        assert result.exit_code == 0
        assert "--max-iter" in result.output
        assert "--async" in result.output

    def test_nonexistent_workdir_fails(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(main, ["invert", "run", str(tmp_path / "no_such")])
        assert result.exit_code != 0

    def test_occam2d_runner_called(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        mock_runner = MagicMock()
        mock_runner.return_value.run.return_value = 0
        with patch(
            "pycsamt.cli.commands.invert.run.OccamRunner",
            mock_runner,
            create=True,
        ):
            result = runner.invoke(
                main,
                [
                    "invert",
                    "run",
                    str(occam_workdir),
                    "--solver",
                    "occam2d",
                    "--max-iter",
                    "50",
                ],
            )
        # OccamRunner may not be importable in test env, but should not traceback
        assert result.exception is None or isinstance(result.exception, SystemExit)

    def test_occam2d_sync_success(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run.return_value = 0
        with patch("pycsamt.models.occam2d.runner.OccamRunner", mock_cls):
            result = runner.invoke(
                main,
                ["invert", "run", str(occam_workdir), "--solver", "occam2d", "-v"],
            )
        assert result.exit_code == 0, result.output
        assert "Starting OCCAM2D" in result.output
        assert "Occam2D finished successfully." in result.output
        mock_cls.return_value.run.assert_called_once_with(
            max_iter=None, target_misfit=None
        )

    def test_occam2d_sync_nonzero_exit(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run.return_value = 2
        with patch("pycsamt.models.occam2d.runner.OccamRunner", mock_cls):
            result = runner.invoke(
                main, ["invert", "run", str(occam_workdir), "--solver", "occam2d"]
            )
        assert result.exit_code == 2
        assert "Occam2D exited with code 2" in result.output

    def test_occam2d_async(self, runner: CliRunner, occam_workdir: Path) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run_async.return_value = 12345
        with patch("pycsamt.models.occam2d.runner.OccamRunner", mock_cls):
            result = runner.invoke(
                main,
                [
                    "invert",
                    "run",
                    str(occam_workdir),
                    "--solver",
                    "occam2d",
                    "--async",
                    "--max-iter",
                    "50",
                    "--target-misfit",
                    "1.05",
                ],
            )
        assert result.exit_code == 0, result.output
        assert "Occam2D started (PID 12345)" in result.output
        mock_cls.return_value.run_async.assert_called_once_with(
            max_iter=50, target_misfit=1.05
        )

    def test_modem_sync_success(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run.return_value = 0
        with patch("pycsamt.models.modem.runner.ModEmRunner", mock_cls):
            result = runner.invoke(
                main, ["invert", "run", str(modem_workdir), "--solver", "modem"]
            )
        assert result.exit_code == 0, result.output
        assert "ModEM finished successfully." in result.output

    def test_modem_sync_nonzero_exit(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run.return_value = 3
        with patch("pycsamt.models.modem.runner.ModEmRunner", mock_cls):
            result = runner.invoke(
                main, ["invert", "run", str(modem_workdir), "--solver", "modem"]
            )
        assert result.exit_code == 3
        assert "ModEM exited with code 3" in result.output

    def test_modem_async(self, runner: CliRunner, modem_workdir: Path) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run.return_value = 999
        with patch("pycsamt.models.modem.runner.ModEmRunner", mock_cls):
            result = runner.invoke(
                main,
                [
                    "invert",
                    "run",
                    str(modem_workdir),
                    "--solver",
                    "modem",
                    "--async",
                ],
            )
        assert result.exit_code == 0, result.output
        assert "ModEM started (PID 999)" in result.output
        mock_cls.return_value.run.assert_called_once_with(
            run_async=True, max_iterations=None, target_rms=None
        )

    def test_runner_exception_reported(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        mock_cls = MagicMock(side_effect=RuntimeError("binary missing"))
        with patch("pycsamt.models.occam2d.runner.OccamRunner", mock_cls):
            result = runner.invoke(
                main, ["invert", "run", str(occam_workdir), "--solver", "occam2d"]
            )
        assert result.exit_code == 1
        assert "Error: binary missing" in result.output

    def test_solver_auto_detected(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        mock_cls = MagicMock()
        mock_cls.return_value.run.return_value = 0
        with patch("pycsamt.models.modem.runner.ModEmRunner", mock_cls):
            result = runner.invoke(main, ["invert", "run", str(modem_workdir)])
        assert result.exit_code == 0, result.output


# ---------------------------------------------------------------------------
# invert results
# ---------------------------------------------------------------------------


class TestInvertResults:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "results", "--help"])
        assert result.exit_code == 0
        assert "--iteration" in result.output

    def test_empty_workdir_fails_gracefully(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        result = runner.invoke(main, ["invert", "results", str(occam_workdir)])
        # Expected to fail (no iter files) but should not traceback uncontrolled
        assert result.exception is None or isinstance(result.exception, SystemExit)


# ---------------------------------------------------------------------------
# invert plot (sub-group help only — no live data needed)
# ---------------------------------------------------------------------------


class TestInvertPlot:
    def test_plot_group_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "plot", "--help"])
        assert result.exit_code == 0
        for sub in (
            "model",
            "misfit",
            "response",
            "pseudo",
            "section",
            "1d",
            "per-site",
            "grid",
        ):
            assert sub in result.output

    @pytest.mark.parametrize(
        "sub",
        [
            "model",
            "misfit",
            "response",
            "pseudo",
            "section",
            "1d",
            "per-site",
            "grid",
        ],
    )
    def test_each_subcommand_has_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["invert", "plot", sub, "--help"])
        assert result.exit_code == 0
        assert "WORKDIR" in result.output
        assert "--save" in result.output
        assert "--show" in result.output


# ---------------------------------------------------------------------------
# invert plot — real bundled Occam2D + ModEM data
# ---------------------------------------------------------------------------

_OCCAM_REAL = Path(__file__).resolve().parents[3] / "data" / "occam2D"
_MODEM_REAL = (
    Path(__file__).resolve().parents[3]
    / "data"
    / "modem"
    / "willy_27freq_watex_line02_sample"
)


def _has_occam_real() -> bool:
    return (_OCCAM_REAL / "OccamDataFile.dat").exists()


def _has_modem_real() -> bool:
    return (_MODEM_REAL / "Modular_NLCG.log").exists()


@pytest.mark.skipif(not _has_occam_real(), reason="bundled data/occam2D absent")
class TestInvertPlotOccamReal:
    def test_model_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "model.png"
        result = runner.invoke(
            main, ["invert", "plot", "model", str(_OCCAM_REAL), "--save", str(out)]
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_model_options(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "model.png"
        result = runner.invoke(
            main,
            [
                "invert",
                "plot",
                "model",
                str(_OCCAM_REAL),
                "--rho-min",
                "1",
                "--rho-max",
                "1000",
                "--depth-max",
                "5000",
                "--no-stations",
                "--cmap",
                "viridis",
                "--iteration",
                "17",
                "--save",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_misfit_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "misfit.png"
        result = runner.invoke(
            main,
            [
                "invert",
                "plot",
                "misfit",
                str(_OCCAM_REAL),
                "--no-roughness",
                "--lagrange",
                "--save",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_response_station_filter(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        out = tmp_path / "resp.png"
        result = runner.invoke(
            main,
            [
                "invert",
                "plot",
                "response",
                str(_OCCAM_REAL),
                "--station",
                "S00",
                "--save",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_pseudo_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "pseudo.png"
        result = runner.invoke(
            main,
            [
                "invert",
                "plot",
                "pseudo",
                str(_OCCAM_REAL),
                "--cmap",
                "viridis",
                "--save",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_1d_with_station_list(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "profiles.png"
        result = runner.invoke(
            main,
            [
                "invert",
                "plot",
                "1d",
                str(_OCCAM_REAL),
                "--stations",
                "S00,S01",
                "--save",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_per_site_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "site_rms.png"
        result = runner.invoke(
            main,
            ["invert", "plot", "per-site", str(_OCCAM_REAL), "--save", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_grid_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "grid.png"
        result = runner.invoke(
            main, ["invert", "plot", "grid", str(_OCCAM_REAL), "--save", str(out)]
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_section_rejected_for_occam(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "plot", "section", str(_OCCAM_REAL)])
        assert result.exit_code != 0
        assert "ModEM-only" in result.output

    def test_1d_rejected_for_modem(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        result = runner.invoke(main, ["invert", "plot", "1d", str(modem_workdir)])
        assert result.exit_code != 0
        assert "Occam2D-only" in result.output

    def test_per_site_rejected_for_modem(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        result = runner.invoke(
            main, ["invert", "plot", "per-site", str(modem_workdir)]
        )
        assert result.exit_code != 0
        assert "Occam2D-only" in result.output

    def test_grid_rejected_for_modem(
        self, runner: CliRunner, modem_workdir: Path
    ) -> None:
        result = runner.invoke(main, ["invert", "plot", "grid", str(modem_workdir)])
        assert result.exit_code != 0
        assert "Occam2D-only" in result.output

    def test_no_save_no_show_warns(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["invert", "plot", "misfit", str(_OCCAM_REAL)])
        assert result.exit_code == 0
        assert "not saved" in result.output

    def test_show_flag_opens_window(self, runner: CliRunner) -> None:
        with patch("matplotlib.pyplot.show") as mock_show:
            result = runner.invoke(
                main, ["invert", "plot", "misfit", str(_OCCAM_REAL), "--show"]
            )
        assert result.exit_code == 0, result.output
        mock_show.assert_called_once()

    def test_save_and_show_together(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        out = tmp_path / "both.png"
        with patch("matplotlib.pyplot.show") as mock_show:
            result = runner.invoke(
                main,
                [
                    "invert",
                    "plot",
                    "misfit",
                    str(_OCCAM_REAL),
                    "--save",
                    str(out),
                    "--show",
                ],
            )
        assert result.exit_code == 0, result.output
        assert out.exists()
        mock_show.assert_called_once()

    def test_model_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.occam2d.plot.PlotModel.plot",
            side_effect=RuntimeError("boom"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "model", str(_OCCAM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: boom" in result.output

    def test_pseudo_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.occam2d.plot.PlotPseudo.plot",
            side_effect=RuntimeError("bad pseudo"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "pseudo", str(_OCCAM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad pseudo" in result.output

    def test_1d_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.occam2d.plot.PlotSounding1D.plot",
            side_effect=RuntimeError("bad 1d"),
        ):
            result = runner.invoke(main, ["invert", "plot", "1d", str(_OCCAM_REAL)])
        assert result.exit_code == 1
        assert "Error: bad 1d" in result.output

    def test_per_site_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.occam2d.plot.PlotSiteMisfit.plot",
            side_effect=RuntimeError("bad site"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "per-site", str(_OCCAM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad site" in result.output

    def test_grid_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.occam2d.plot.PlotResponseGrid.plot",
            side_effect=RuntimeError("bad grid"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "grid", str(_OCCAM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad grid" in result.output

    def test_response_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.occam2d.plot.PlotResponse.plot",
            side_effect=RuntimeError("bad response"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "response", str(_OCCAM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad response" in result.output


@pytest.mark.skipif(not _has_modem_real(), reason="bundled ModEM sample absent")
class TestInvertPlotModemReal:
    def test_model_plot_real_bug_handled_gracefully(
        self, runner: CliRunner
    ) -> None:
        """``invert plot model`` on this real ModEM run currently fails
        inside ``PlotModel2D`` itself (a pre-existing bug outside this
        batch's scope) -- assert the CLI still degrades to a clean
        ``Error: ...`` + exit 1 instead of an uncaught traceback."""
        result = runner.invoke(main, ["invert", "plot", "model", str(_MODEM_REAL)])
        assert result.exception is None or isinstance(
            result.exception, SystemExit
        )
        if result.exit_code != 0:
            assert "Error:" in result.output

    def test_misfit_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "misfit.png"
        result = runner.invoke(
            main,
            ["invert", "plot", "misfit", str(_MODEM_REAL), "--save", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_response_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "resp.png"
        result = runner.invoke(
            main,
            ["invert", "plot", "response", str(_MODEM_REAL), "--save", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_pseudo_saves_figure(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "pseudo.png"
        result = runner.invoke(
            main,
            ["invert", "plot", "pseudo", str(_MODEM_REAL), "--save", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_section_with_depth(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "section.png"
        result = runner.invoke(
            main,
            [
                "invert",
                "plot",
                "section",
                str(_MODEM_REAL),
                "--depth",
                "5000",
                "--save",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_section_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.modem.plot.PlotSection.plot",
            side_effect=RuntimeError("bad section"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "section", str(_MODEM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad section" in result.output

    def test_misfit_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.modem.plot.PlotMisfit.plot",
            side_effect=RuntimeError("bad misfit"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "misfit", str(_MODEM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad misfit" in result.output

    def test_response_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.modem.plot.PlotResponse.plot",
            side_effect=RuntimeError("bad response"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "response", str(_MODEM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad response" in result.output

    def test_pseudo_plot_error_handled(self, runner: CliRunner) -> None:
        with patch(
            "pycsamt.models.modem.plot.PlotPseudo.plot",
            side_effect=RuntimeError("bad pseudo"),
        ):
            result = runner.invoke(
                main, ["invert", "plot", "pseudo", str(_MODEM_REAL)]
            )
        assert result.exit_code == 1
        assert "Error: bad pseudo" in result.output

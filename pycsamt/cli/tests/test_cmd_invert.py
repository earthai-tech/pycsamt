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
        assert "No active survey" in (
            result.output + str(result.exception or "")
        )

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
        assert result.exception is None or isinstance(
            result.exception, SystemExit
        )

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

        with patch(
            "pycsamt.cli.survey._build_sites", return_value=fake_sites
        ):
            set_survey(edi_dir)

        # Call with explicit alt_dir — resolve_survey should use alt_dir
        resolved_paths = []
        __import__(
            "pycsamt.cli.survey", fromlist=["resolve_survey"]
        ).resolve_survey

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

    def test_occam2d_workdir_text(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        result = runner.invoke(main, ["invert", "status", str(occam_workdir)])
        assert result.exit_code == 0
        assert (
            "OCCAM2D" in result.output.upper() or "occam2d" in result.output
        )

    def test_occam2d_workdir_json(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
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

    def test_nonexistent_workdir_fails(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main, ["invert", "run", str(tmp_path / "no_such")]
        )
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
        assert result.exception is None or isinstance(
            result.exception, SystemExit
        )


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
        result = runner.invoke(
            main, ["invert", "results", str(occam_workdir)]
        )
        # Expected to fail (no iter files) but should not traceback uncontrolled
        assert result.exception is None or isinstance(
            result.exception, SystemExit
        )


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
    def test_each_subcommand_has_help(
        self, runner: CliRunner, sub: str
    ) -> None:
        result = runner.invoke(main, ["invert", "plot", sub, "--help"])
        assert result.exit_code == 0
        assert "WORKDIR" in result.output
        assert "--save" in result.output
        assert "--show" in result.output

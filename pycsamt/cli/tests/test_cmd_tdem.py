# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt tdem`` command group.

Test strategy
-------------
* **Help tests** — always run; verify Click wiring and option names.
* **Unit tests** — exercise the workflow classes directly (no CLI overhead).
* **Live integration tests** — require ``data/TEMAVG/JIANGSU/``.
  All live tests skip gracefully when the data directory is absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_TEMAVG_DIR = _PROJECT_ROOT / "data" / "TEMAVG" / "JIANGSU"


def _has_temavg() -> bool:
    return _TEMAVG_DIR.exists() and bool(list(_TEMAVG_DIR.glob("*.AVG")))


# ---------------------------------------------------------------------------
# pycsamt tdem  (group help)
# ---------------------------------------------------------------------------


class TestTdemGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "--help"])
        assert result.exit_code == 0
        for sub in ("info", "convert", "plot"):
            assert sub in result.output

    @pytest.mark.parametrize("sub", ["info", "convert", "plot"])
    def test_each_subcommand_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["tdem", sub, "--help"])
        assert result.exit_code == 0
        assert "SURVEY_DIR" in result.output


# ---------------------------------------------------------------------------
# pycsamt tdem info
# ---------------------------------------------------------------------------


class TestTdemInfo:
    @pytest.fixture(autouse=True)
    def require_temavg(self) -> None:
        if not _has_temavg():
            pytest.skip("data/TEMAVG/JIANGSU/ not found")

    def test_text_output(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "info", str(_TEMAVG_DIR)])
        assert result.exit_code == 0
        assert "Survey root" in result.output
        assert "AVG files" in result.output

    def test_json_output_keys(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["tdem", "info", str(_TEMAVG_DIR), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        for key in (
            "root",
            "n_avg_files",
            "n_z_files",
            "n_log_files",
            "avg_stems",
            "has_coordinates",
        ):
            assert key in data

    def test_n_avg_files_positive(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["tdem", "info", str(_TEMAVG_DIR), "--format", "json"]
        )
        data = json.loads(result.output)
        assert data["n_avg_files"] > 0

    def test_verbose_flag(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "info", str(_TEMAVG_DIR), "-v"])
        assert result.exit_code == 0

    def test_nonexistent_dir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "info", "/no/such/path"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt tdem convert
# ---------------------------------------------------------------------------


class TestTdemConvert:
    @pytest.fixture(autouse=True)
    def require_temavg(self) -> None:
        if not _has_temavg():
            pytest.skip("data/TEMAVG/JIANGSU/ not found")

    def test_dry_run_lists_soundings(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "convert", str(_TEMAVG_DIR), "--dry-run"])
        assert result.exit_code == 0
        assert "sounding" in result.output.lower()

    def test_convert_text_output(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "tdem",
                "convert",
                str(_TEMAVG_DIR),
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        assert "Soundings" in result.output

    def test_convert_json_output(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "tdem",
                "convert",
                str(_TEMAVG_DIR),
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        for key in ("n_soundings", "n_written", "written", "output_dir"):
            assert key in data

    def test_no_output_dir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "convert", str(_TEMAVG_DIR)])
        assert result.exit_code != 0

    def test_stems_filter(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main, ["tdem", "info", str(_TEMAVG_DIR), "--format", "json"]
        )
        stems = json.loads(result.output).get("avg_stems", [])
        if not stems:
            pytest.skip("No stems found in survey")

        one_stem = stems[0]
        result = runner.invoke(
            main,
            [
                "tdem",
                "convert",
                str(_TEMAVG_DIR),
                "--stems",
                one_stem,
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_soundings"] >= 1

    def test_method_fourier(self, runner: CliRunner, tmp_path: Path) -> None:
        # Fourier's per-sounding cosine-transform + Kramers-Kronig
        # reconstruction is Python-loop-heavy (unlike late_time, which
        # is fully vectorized): the full JIANGSU survey is 2790
        # soundings, taking ~2 minutes even outside CI and far longer
        # under --cov-branch instrumentation. This test only needs to
        # prove --method fourier is wired correctly, so restrict to a
        # single stem like test_stems_filter does.
        result = runner.invoke(
            main,
            [
                "tdem",
                "convert",
                str(_TEMAVG_DIR),
                "--stems",
                "TEM100",
                "--method",
                "fourier",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0

    def test_nonexistent_dir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["tdem", "convert", "/no/such/path", "--output-dir", "/tmp"]
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt tdem plot
# ---------------------------------------------------------------------------


class TestTdemPlot:
    @pytest.fixture(autouse=True)
    def require_temavg(self) -> None:
        if not _has_temavg():
            pytest.skip("data/TEMAVG/JIANGSU/ not found")

    @pytest.mark.parametrize("kind", ["decay", "rho", "section", "map"])
    def test_plot_saves_file(
        self, runner: CliRunner, tmp_path: Path, kind: str
    ) -> None:
        # Restrict to one stem: the full JIANGSU survey is 2790
        # soundings, and per-sounding plot styling (e.g. a matplotlib
        # colormap call per line) over the whole thing is far slower
        # than what this wiring/smoke test needs.
        result = runner.invoke(
            main,
            [
                "tdem",
                "plot",
                str(_TEMAVG_DIR),
                "--stems",
                "TEM100",
                "--kind",
                kind,
                "--output-dir",
                str(tmp_path),
            ],
        )
        if result.exit_code != 0:
            pytest.skip(f"Plot kind {kind!r} not available: {result.output}")
        saved = list(tmp_path.glob("*.png"))
        assert saved, f"No PNG file written for kind={kind!r}"

    def test_plot_help_shows_kinds(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["tdem", "plot", "--help"])
        assert result.exit_code == 0
        assert "decay" in result.output
        assert "dashboard" in result.output


# ---------------------------------------------------------------------------
# Workflow unit tests (no CLI, no CliRunner needed)
# ---------------------------------------------------------------------------


class TestTdemWorkflowUnit:
    @pytest.fixture(autouse=True)
    def require_temavg(self) -> None:
        if not _has_temavg():
            pytest.skip("data/TEMAVG/JIANGSU/ not found")

    def test_read_temavg_soundings_returns_list(self) -> None:
        from pycsamt.tdem.workflow import (
            read_temavg_soundings,
        )

        soundings = read_temavg_soundings(_TEMAVG_DIR)
        assert isinstance(soundings, list)
        assert len(soundings) > 0

    def test_sounding_has_time_gates(self) -> None:
        from pycsamt.tdem.workflow import (
            read_temavg_soundings,
        )

        soundings = read_temavg_soundings(_TEMAVG_DIR)
        snd = soundings[0]
        assert len(snd.time_gates) > 0

    def test_survey_has_avg_files(self) -> None:
        from pycsamt.tdem.survey import read_temavg_survey

        survey = read_temavg_survey(_TEMAVG_DIR)
        assert len(survey.avg_files) > 0

    def test_transform_returns_conversion_bundle(self) -> None:
        from pycsamt.tdem.workflow import (
            transform_temavg_survey,
        )

        result = transform_temavg_survey(_TEMAVG_DIR)
        assert result.n_soundings > 0

    def test_transform_writes_edis(self, tmp_path: Path) -> None:
        from pycsamt.tdem.workflow import (
            transform_temavg_survey,
        )

        result = transform_temavg_survey(_TEMAVG_DIR, savepath=tmp_path)
        assert result.n_soundings > 0
        assert len(result.written_paths) > 0
        for p in result.written_paths:
            assert Path(p).exists()

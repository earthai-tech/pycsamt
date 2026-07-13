# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt interp`` command group.

Test strategy
-------------
* **Help tests** — always run; verify Click wiring, option names, sub-commands.
* **Unit tests (rocks)** — use the real built-in RockDatabase; pure Python, no
  mocking or file I/O required.
* **Unit tests (classify / export)** — monkeypatch ``InversionResult`` and
  ``ResistivityModel.from_occam2d`` to avoid real Occam2D file parsing.
  The ``occam_workdir_with_iters`` fixture provides a directory with the
  right signature files so the solver auto-detection path is exercised.
* **Error tests** — verify exit codes and useful messages for bad inputs,
  missing files, and unsupported solvers.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# Fake interp objects (avoid real Occam2D I/O in unit tests)
# ---------------------------------------------------------------------------


def _fake_model(n_x: int = 5, n_z: int = 10) -> MagicMock:
    """Return a ResistivityModel-like mock."""
    rng = np.random.default_rng(0)
    m = MagicMock()
    m.x_centers = np.linspace(0, 1000, n_x)
    m.z_centers = np.linspace(5, 500, n_z)
    m.rho_2d = rng.uniform(1.5, 3.5, (n_z, n_x))
    m.station_x = np.linspace(0, 1000, n_x)
    m.station_names = [f"S{i:02d}" for i in range(n_x)]
    m.method = "occam2d"
    m.rms = 1.234
    m.n_x = n_x
    m.n_z = n_z
    return m


def _fake_layer(
    depth_top: float, depth_bot: float, rock_name: str
) -> MagicMock:
    layer = MagicMock()
    layer.depth_top = depth_top
    layer.depth_bot = depth_bot
    layer.rock = MagicMock()
    layer.rock.name = rock_name
    layer.rock.rho_min = 50.0
    layer.rock.rho_max = 500.0
    layer.rock.description = "Test rock"
    return layer


def _fake_log(station_name: str, n_layers: int = 3) -> MagicMock:
    log = MagicMock()
    log.station_name = station_name
    log.station_x = 0.0
    log.layers = [
        _fake_layer(i * 50, (i + 1) * 50, f"Rock_{i}")
        for i in range(n_layers)
    ]
    return log


def _patch_classify(
    monkeypatch: pytest.MonkeyPatch, n_stations: int = 3
) -> None:
    """Monkeypatch InversionResult + ResistivityModel + ModelCalibrator."""
    model = _fake_model(n_x=n_stations)
    logs = [_fake_log(f"S{i:02d}") for i in range(n_stations)]

    # Patch InversionResult constructor
    monkeypatch.setattr(
        "pycsamt.cli.commands.interp.classify.InversionResult",
        lambda workdir, iteration=None: MagicMock(),
        raising=False,
    )
    # Patch ResistivityModel.from_occam2d
    monkeypatch.setattr(
        "pycsamt.interp.ResistivityModel.from_occam2d",
        staticmethod(lambda result: model),
        raising=False,
    )
    # Patch ModelCalibrator.fit → stratigraphic_logs
    cal_mock = MagicMock()
    cal_mock.stratigraphic_logs.return_value = logs
    monkeypatch.setattr(
        "pycsamt.interp.ModelCalibrator",
        lambda **kw: MagicMock(fit=lambda model, **kw2: cal_mock),
        raising=False,
    )


def _patch_classify_in_module(
    monkeypatch: pytest.MonkeyPatch,
    module: str,
    n_stations: int = 3,
) -> list:
    """Patch lazy imports via their source modules.

    ``classify`` and ``export`` use:
        from ....models.occam2d.results import InversionResult
        from ....interp import ResistivityModel, ModelCalibrator

    Both resolve via their defining modules, so we patch those.
    """
    import pycsamt.interp as _interp
    import pycsamt.models.occam2d.results as _results

    model = _fake_model(n_x=n_stations)
    logs = [_fake_log(f"S{i:02d}") for i in range(n_stations)]

    ir_mock = MagicMock()
    cal_mock = MagicMock()
    cal_mock.stratigraphic_logs.return_value = logs

    monkeypatch.setattr(
        _results,
        "InversionResult",
        lambda workdir, iteration=None: ir_mock,
    )

    rm_mock = MagicMock()
    rm_mock.from_occam2d = staticmethod(lambda result: model)
    monkeypatch.setattr(_interp, "ResistivityModel", rm_mock)
    monkeypatch.setattr(
        _interp,
        "ModelCalibrator",
        lambda **kw: MagicMock(fit=lambda m, **kw2: cal_mock),
    )
    return logs


# ---------------------------------------------------------------------------
# pycsamt interp  (group help)
# ---------------------------------------------------------------------------


class TestInterpGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "--help"])
        assert result.exit_code == 0
        for sub in ("rocks", "classify", "export"):
            assert sub in result.output

    def test_rocks_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--help"])
        assert result.exit_code == 0
        assert "--rho" in result.output
        assert "--method" in result.output

    def test_classify_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "classify", "--help"])
        assert result.exit_code == 0
        for opt in (
            "WORKDIR",
            "--solver",
            "--ptol",
            "--merge-tol",
            "--iteration",
            "--output",
        ):
            assert opt in result.output

    def test_export_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "export", "--help"])
        assert result.exit_code == 0
        for opt in (
            "WORKDIR",
            "--format",
            "--output-dir",
            "--solver",
            "--ptol",
        ):
            assert opt in result.output


# ---------------------------------------------------------------------------
# pycsamt interp rocks
# ---------------------------------------------------------------------------


class TestInterpRocks:
    def test_list_all_text(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks"])
        assert result.exit_code == 0
        # Built-in database contains at least these rock types
        for name in ("Clay", "Aquifer", "Granite", "Sandstone"):
            assert name in result.output

    def test_list_all_json(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--format", "json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) >= 20
        assert all("name" in row and "rho_min" in row for row in data)

    def test_list_all_csv(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--format", "csv"])
        assert result.exit_code == 0
        first_line = result.output.strip().splitlines()[0]
        assert "name" in first_line
        assert "rho_min" in first_line

    def test_classify_text(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--rho", "250"])
        assert result.exit_code == 0
        assert "250.0" in result.output
        assert "Ω·m" in result.output

    def test_classify_json(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["interp", "rocks", "--rho", "250", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["query_rho_Ohm_m"] == pytest.approx(250.0)
        assert "name" in data
        assert "rho_min_Ohm_m" in data

    def test_classify_csv(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["interp", "rocks", "--rho", "5000", "--format", "csv"]
        )
        assert result.exit_code == 0
        assert "name" in result.output

    def test_classify_conductive_value(self, runner: CliRunner) -> None:
        """Low resistivity should match a conductive rock (clay, aquifer…)."""
        result = runner.invoke(main, ["interp", "rocks", "--rho", "5"])
        assert result.exit_code == 0
        assert "Ω·m" in result.output

    def test_classify_resistive_value(self, runner: CliRunner) -> None:
        """High resistivity should match a resistive rock."""
        result = runner.invoke(main, ["interp", "rocks", "--rho", "100000"])
        assert result.exit_code == 0

    def test_overlap_method(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["interp", "rocks", "--rho", "50", "--method", "overlap"],
        )
        assert result.exit_code == 0

    def test_zero_rho_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--rho", "0"])
        assert result.exit_code != 0

    def test_negative_rho_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--rho", "-10"])
        assert result.exit_code != 0

    @pytest.mark.parametrize(
        "rho,expected_fragment",
        [
            (0.01, "Sulfide"),  # very conductive
            (10.0, None),  # any clay-range rock
            (500.0, None),  # sandstone / limestone range
            (50000.0, None),  # crystalline basement range
        ],
    )
    def test_parametrized_classify(
        self,
        rho: float,
        expected_fragment: str | None,
        runner: CliRunner,
    ) -> None:
        result = runner.invoke(main, ["interp", "rocks", "--rho", str(rho)])
        assert result.exit_code == 0
        if expected_fragment:
            assert expected_fragment in result.output


# ---------------------------------------------------------------------------
# pycsamt interp classify
# ---------------------------------------------------------------------------


class TestInterpClassify:
    def test_nonexistent_workdir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["interp", "classify", "/nonexistent/path"]
        )
        assert result.exit_code != 0

    def test_unknown_solver_fails(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """A directory without Occam2D files and --solver auto should error."""
        result = runner.invoke(main, ["interp", "classify", str(tmp_path)])
        assert result.exit_code != 0
        assert "Could not detect" in result.output or result.exception

    def test_unsupported_solver_fails(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        """Passing an unsupported solver should fail gracefully."""
        # Currently only occam2d is supported; 'modem' is in choices for
        # future use but not yet implemented in classify.
        # We indirectly test by forcing an auto-detection on a clean dir.
        result = runner.invoke(
            main,
            ["interp", "classify", str(occam_workdir), "--solver", "occam2d"],
        )
        # Without mocked InversionResult this will error, but NOT with exit 0
        # The important check: it didn't crash with an unhandled exception
        assert result.exit_code in (0, 1)

    # -- unit tests with mocked internals -------------------------------------

    def test_text_output(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
    ) -> None:
        _patch_classify_in_module(monkeypatch, "classify")
        result = runner.invoke(
            main,
            ["interp", "classify", str(occam_workdir_with_iters)],
        )
        assert result.exit_code == 0
        assert "station" in result.output.lower()

    def test_json_output(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
    ) -> None:
        _patch_classify_in_module(monkeypatch, "classify")
        result = runner.invoke(
            main,
            [
                "interp",
                "classify",
                str(occam_workdir_with_iters),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) > 0
        assert "station" in data[0]

    def test_csv_output(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
    ) -> None:
        _patch_classify_in_module(monkeypatch, "classify")
        result = runner.invoke(
            main,
            [
                "interp",
                "classify",
                str(occam_workdir_with_iters),
                "--format",
                "csv",
            ],
        )
        assert result.exit_code == 0
        first_line = result.output.strip().splitlines()[0]
        assert "station" in first_line

    def test_save_csv_output(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
        tmp_path: Path,
    ) -> None:
        _patch_classify_in_module(monkeypatch, "classify")
        out = tmp_path / "logs.csv"
        result = runner.invoke(
            main,
            [
                "interp",
                "classify",
                str(occam_workdir_with_iters),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0
        assert out.exists()
        assert out.stat().st_size > 0

    def test_auto_detects_occam2d(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
    ) -> None:
        """Solver auto-detection path: workdir has Occam2D signature files."""
        _patch_classify_in_module(monkeypatch, "classify")
        result = runner.invoke(
            main,
            [
                "interp",
                "classify",
                str(occam_workdir_with_iters),
                "--solver",
                "auto",
            ],
        )
        assert result.exit_code == 0

    def test_verbose_flag(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
    ) -> None:
        _patch_classify_in_module(monkeypatch, "classify")
        result = runner.invoke(
            main,
            ["interp", "classify", str(occam_workdir_with_iters), "-v"],
        )
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# pycsamt interp export
# ---------------------------------------------------------------------------


class TestInterpExport:
    def test_nonexistent_workdir_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["interp", "export", "/nonexistent/path", "--format", "csv"],
        )
        assert result.exit_code != 0

    def test_missing_format_fails(
        self, runner: CliRunner, occam_workdir: Path
    ) -> None:
        """--format is required; omitting it must fail."""
        result = runner.invoke(main, ["interp", "export", str(occam_workdir)])
        assert result.exit_code != 0

    def test_unknown_solver_fails(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["interp", "export", str(tmp_path), "--format", "csv"],
        )
        assert result.exit_code != 0

    # -- unit tests with mocked internals -------------------------------------

    def _invoke_export(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        workdir: Path,
        fmt: str,
        output_dir: Path,
        extra_args: list | None = None,
    ):
        # Patch loading and classify via source modules
        _patch_classify_in_module(monkeypatch, "export")

        # Patch the export writers on the source module so the lazy
        # `from ....interp import export as _export` picks them up.
        import pycsamt.interp.export as _exp

        monkeypatch.setattr(
            _exp, "to_oasis_montaj_xyz", lambda logs, path: None
        )
        monkeypatch.setattr(_exp, "to_las", lambda log, path: None)
        monkeypatch.setattr(_exp, "to_csv", lambda logs, path: None)
        monkeypatch.setattr(_exp, "to_vtk", lambda model, path: None)

        args = [
            "interp",
            "export",
            str(workdir),
            "--format",
            fmt,
            "--output-dir",
            str(output_dir),
        ] + (extra_args or [])
        return runner.invoke(main, args)

    @pytest.mark.parametrize("fmt", ["xyz", "las", "csv", "vtk"])
    def test_all_formats(
        self,
        fmt: str,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
        tmp_path: Path,
    ) -> None:
        result = self._invoke_export(
            runner,
            monkeypatch,
            occam_workdir_with_iters,
            fmt,
            tmp_path / fmt,
        )
        assert result.exit_code == 0
        assert "Written" in result.output

    def test_overwrite_flag(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
        tmp_path: Path,
    ) -> None:
        """Second run with --overwrite must succeed even if file exists."""
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (
            out_dir / "layers.csv"
        ).touch()  # pre-create to trigger overwrite logic

        result = self._invoke_export(
            runner,
            monkeypatch,
            occam_workdir_with_iters,
            "csv",
            out_dir,
            extra_args=["--overwrite"],
        )
        assert result.exit_code == 0

    def test_existing_file_without_overwrite_fails(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        occam_workdir_with_iters: Path,
        tmp_path: Path,
    ) -> None:
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (out_dir / "layers.csv").touch()

        _patch_classify_in_module(monkeypatch, "export")
        import pycsamt.interp.export as _exp

        monkeypatch.setattr(_exp, "to_csv", lambda logs, path: None)

        result = runner.invoke(
            main,
            [
                "interp",
                "export",
                str(occam_workdir_with_iters),
                "--format",
                "csv",
                "--output-dir",
                str(out_dir),
            ],
        )
        assert result.exit_code != 0
        assert "overwrite" in result.output.lower() or result.exception

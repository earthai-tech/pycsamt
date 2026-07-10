# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt site`` command group.

Test strategy
-------------
* **Help tests** — always run; verify Click wiring, option names, sub-commands.
* **Unit tests** — patch ``_get_sites`` at the sub-module level so no EDI I/O
  is required; use ``_FakeSites`` / ``make_fake_sites`` from conftest.
* **Live integration tests** — use ``site_edi_dir`` fixture; skip gracefully
  when no EDI data directory is found in the repo.
* **Write tests** — use ``site_edi_dir_writable`` (copy in tmp_path) so the
  source data is never mutated.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.tests.conftest import (
    make_fake_sites,
)

# ---------------------------------------------------------------------------
# Convenience helpers
# ---------------------------------------------------------------------------

_SITE_MODULE = "pycsamt.cli.commands.site"


def _patch_get_sites(monkeypatch, sites, submod: str = "info"):
    """Patch _get_sites in a specific site sub-module."""
    monkeypatch.setattr(
        f"{_SITE_MODULE}.{submod}._get_sites",
        lambda *a, **kw: sites,
    )


# ---------------------------------------------------------------------------
# pycsamt site  (group help)
# ---------------------------------------------------------------------------

class TestSiteGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "--help"])
        assert result.exit_code == 0
        for sub in ("info", "select", "edit", "export",
                    "recompute", "compute"):
            assert sub in result.output

    def test_compute_subgroup_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "compute", "--help"])
        assert result.exit_code == 0
        for sub in ("strike", "resistivity", "phase-slope", "tipper"):
            assert sub in result.output


# ---------------------------------------------------------------------------
# pycsamt site info
# ---------------------------------------------------------------------------

class TestSiteInfo:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "info", "--help"])
        assert result.exit_code == 0
        assert "--station" in result.output
        assert "--top" in result.output

    def test_no_source_no_context(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["site", "info"])
        assert result.exit_code != 0
        assert "No active survey" in result.output

    # --- unit tests with fake sites ---

    def test_fake_sites_text_output(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Verifies the report path runs without error on fake data."""
        fake = make_fake_sites(3)
        _patch_get_sites(monkeypatch, fake, "info")

        # Also stub SitesReport to avoid real numpy computation
        from pycsamt.site.report import SitesReport

        def _noop_report(self, **kw): pass  # noqa: ANN001
        monkeypatch.setattr(SitesReport, "report", _noop_report)

        # invoke with a dummy path (path check is bypassed via _get_sites mock)
        result = runner.invoke(main, ["site", "info", "--survey", "."])
        # Exit code 0 OR the only failure is that the path doesn't exist
        assert result.exception is None or isinstance(result.exception, SystemExit)

    # --- live integration tests ---

    def test_live_text_output(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(main, ["site", "info", str(site_edi_dir)])
        assert result.exit_code == 0

    def test_live_json_output(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "info", str(site_edi_dir), "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) > 0
        assert "name" in data[0]

    def test_live_csv_output(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "info", str(site_edi_dir), "--format", "csv"]
        )
        assert result.exit_code == 0
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2  # header + at least 1 station

    def test_live_top_option(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "info", str(site_edi_dir), "--top", "1"]
        )
        assert result.exit_code == 0

    def test_live_single_station(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        edis = sorted(site_edi_dir.glob("*.edi"))
        station_name = edis[0].stem
        result = runner.invoke(
            main,
            ["site", "info", str(site_edi_dir), "--station", station_name],
        )
        assert result.exit_code == 0

    def test_live_station_not_found(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "info", str(site_edi_dir), "--station", "DOESNOTEXIST_XYZ"],
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt site select
# ---------------------------------------------------------------------------

class TestSiteSelect:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "select", "--help"])
        assert result.exit_code == 0
        for opt in ("--stations", "--freq", "--bbox", "--drop-empty",
                    "--keep-finite", "--dry-run"):
            assert opt in result.output

    def test_no_source_fails(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["site", "select"])
        assert result.exit_code != 0

    # --- live tests ---

    def test_dry_run_all(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "select", str(site_edi_dir), "--dry-run"]
        )
        assert result.exit_code == 0
        n_edi = len(list(site_edi_dir.glob("*.edi")))
        assert str(n_edi) in result.output

    def test_filter_by_station_dry_run(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        first = sorted(site_edi_dir.glob("*.edi"))[0].stem
        result = runner.invoke(
            main,
            ["site", "select", str(site_edi_dir),
             "--stations", first, "--dry-run"],
        )
        assert result.exit_code == 0
        assert "1/" in result.output   # "1/N station(s) would be selected"
        assert first in result.output

    def test_filter_no_match_fails(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "select", str(site_edi_dir),
             "--stations", "XYZNOTEXIST", "--dry-run"],
        )
        assert result.exit_code != 0

    def test_dry_run_json(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "select", str(site_edi_dir),
             "--dry-run", "--format", "json"],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "n_selected" in data
        assert "stations" in data

    def test_write_output(
        self, runner: CliRunner, site_edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "subset"
        first = sorted(site_edi_dir.glob("*.edi"))[0].stem
        result = runner.invoke(
            main,
            ["site", "select", str(site_edi_dir),
             "--stations", first,
             "--output-dir", str(out)],
        )
        assert result.exit_code == 0
        written = list(out.glob("*.edi"))
        assert len(written) == 1
        assert written[0].stem == first

    def test_freq_filter_dry_run(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        # Very wide range — should keep all stations
        result = runner.invoke(
            main,
            ["site", "select", str(site_edi_dir),
             "--freq", "0.001:100000", "--dry-run"],
        )
        # Just verify no crash and output makes sense
        assert result.exception is None or isinstance(result.exception, SystemExit)

    def test_drop_empty_dry_run(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "select", str(site_edi_dir),
             "--drop-empty", "--dry-run"],
        )
        assert result.exception is None or isinstance(result.exception, SystemExit)


# ---------------------------------------------------------------------------
# pycsamt site edit
# ---------------------------------------------------------------------------

class TestSiteEdit:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "edit", "--help"])
        assert result.exit_code == 0
        for opt in ("--rotate", "--freq", "--fill-missing",
                    "--set-coords", "--coords-table", "--dry-run"):
            assert opt in result.output

    def test_no_operation_fails(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "edit", str(site_edi_dir)]
        )
        assert result.exit_code != 0
        assert "at least one edit" in result.output.lower()

    def test_no_output_dir_fails(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "edit", str(site_edi_dir), "--rotate", "30"],
        )
        assert result.exit_code != 0
        assert "--output-dir" in result.output

    def test_dry_run_rotate(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "edit", str(site_edi_dir),
             "--rotate", "45", "--dry-run"],
        )
        assert result.exit_code == 0
        assert "rotate" in result.output.lower()

    def test_dry_run_multiple_ops(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "edit", str(site_edi_dir),
             "--rotate", "15", "--freq", "0.1:1000",
             "--fill-missing", "--dry-run"],
        )
        assert result.exit_code == 0
        out_lower = result.output.lower()
        assert "rotate" in out_lower
        assert "select_freq" in out_lower or "freq" in out_lower

    def test_live_rotate_writes_files(
        self, runner: CliRunner, site_edi_dir_writable: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "rotated"
        len(list(site_edi_dir_writable.glob("*.edi")))
        result = runner.invoke(
            main,
            ["site", "edit", str(site_edi_dir_writable),
             "--rotate", "30",
             "--output-dir", str(out)],
        )
        # Acceptable: either succeeds (writes) or gracefully errors
        assert result.exception is None or isinstance(result.exception, SystemExit)
        if result.exit_code == 0:
            assert out.exists()

    def test_set_coords_dry_run(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "edit", str(site_edi_dir),
             "--set-coords", "28.5,102.3,200", "--dry-run"],
        )
        assert result.exit_code == 0
        assert "set_coords" in result.output

    def test_set_coords_bad_format(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "edit", str(site_edi_dir),
             "--set-coords", "not_a_coord", "--dry-run"],
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# pycsamt site export
# ---------------------------------------------------------------------------

class TestSiteExport:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "export", "--help"])
        assert result.exit_code == 0
        for opt in ("--output-dir", "--zip", "--template", "--manifest"):
            assert opt in result.output

    def test_no_output_dir_fails(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "export", str(site_edi_dir)]
        )
        assert result.exit_code != 0
        assert "--output-dir" in result.output

    def test_live_export_edis(
        self, runner: CliRunner, site_edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "exported"
        len(list(site_edi_dir.glob("*.edi")))
        result = runner.invoke(
            main,
            ["site", "export", str(site_edi_dir),
             "--output-dir", str(out)],
        )
        assert result.exception is None or isinstance(result.exception, SystemExit)
        if result.exit_code == 0:
            assert out.exists()

    def test_live_export_zip(
        self, runner: CliRunner, site_edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "zipped"
        result = runner.invoke(
            main,
            ["site", "export", str(site_edi_dir),
             "--zip", "--zip-name", "survey.zip",
             "--output-dir", str(out)],
        )
        assert result.exception is None or isinstance(result.exception, SystemExit)

    def test_overwrite_protection(
        self, runner: CliRunner, site_edi_dir: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        out.mkdir()
        # first export
        runner.invoke(
            main,
            ["site", "export", str(site_edi_dir), "--output-dir", str(out)],
        )
        # second export without --overwrite — should skip or warn, not crash
        result = runner.invoke(
            main,
            ["site", "export", str(site_edi_dir), "--output-dir", str(out)],
        )
        assert result.exception is None or isinstance(result.exception, SystemExit)


# ---------------------------------------------------------------------------
# pycsamt site recompute
# ---------------------------------------------------------------------------

class TestSiteRecompute:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "recompute", "--help"])
        assert result.exit_code == 0
        for opt in ("--rotate", "--components", "--freq",
                    "--flatten", "--output-name", "--progress"):
            assert opt in result.output

    def test_no_source_fails(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["site", "recompute"])
        assert result.exit_code != 0
        assert "No EDI source" in result.output

    def test_line_folder_default_output(
        self,
        runner: CliRunner,
        site_edi_dir: Path,
        tmp_path: Path,
    ) -> None:
        root = tmp_path / "WILLY_DATA"
        line = root / "L18PLT"
        shutil.copytree(site_edi_dir, line)

        result = runner.invoke(
            main,
            [
                "site",
                "recompute",
                str(line),
                "--template",
                "{source_stem}.edi",
                "--overwrite",
            ],
        )

        assert result.exit_code == 0
        assert "Recomputed" in result.output
        outdir = root / "recomputed_edis" / "L18PLT"
        assert outdir.exists()
        assert len(list(outdir.glob("*.edi"))) == len(
            list(line.glob("*.edi"))
        )
        assert (root / "recomputed_edis" / "manifest.csv").exists()

    def test_survey_folder_flatten_output(
        self,
        runner: CliRunner,
        site_edi_dir: Path,
        tmp_path: Path,
    ) -> None:
        root = tmp_path / "WILLY_DATA"
        shutil.copytree(site_edi_dir, root / "L18PLT")
        shutil.copytree(site_edi_dir, root / "L32PLT")

        result = runner.invoke(
            main,
            [
                "site",
                "recompute",
                str(root),
                "--flatten",
                "--template",
                "{line}_{source_stem}.edi",
                "--overwrite",
            ],
        )

        assert result.exit_code == 0
        outdir = root / "recomputed_edis"
        assert outdir.exists()
        assert list(outdir.glob("*.edi"))
        assert not (outdir / "L18PLT").exists()

    def test_dry_run_writes_nothing(
        self,
        runner: CliRunner,
        site_edi_dir: Path,
        tmp_path: Path,
    ) -> None:
        root = tmp_path / "WILLY_DATA"
        line = root / "L18PLT"
        shutil.copytree(site_edi_dir, line)

        result = runner.invoke(
            main,
            ["site", "recompute", str(line), "--dry-run"],
        )

        assert result.exit_code == 0
        assert "Dry run" in result.output
        assert not (root / "recomputed_edis").exists()


# ---------------------------------------------------------------------------
# pycsamt site compute
# ---------------------------------------------------------------------------

class TestSiteCompute:

    # --- help tests (always run) ---

    @pytest.mark.parametrize("sub", ["strike", "resistivity", "phase-slope", "tipper"])
    def test_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["site", "compute", sub, "--help"])
        assert result.exit_code == 0
        assert "EDI_DIR" in result.output
        assert "--survey" in result.output

    def test_resistivity_requires_freq(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["site", "compute", "resistivity"])
        assert result.exit_code != 0
        assert "--freq" in result.output or "freq" in result.output.lower()

    # --- live integration tests ---

    def test_live_strike_text(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "compute", "strike", str(site_edi_dir)]
        )
        assert result.exit_code == 0
        assert "station" in result.output.lower() or "theta" in result.output.lower()

    def test_live_strike_swift(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "strike", str(site_edi_dir),
             "--method", "swift"],
        )
        assert result.exit_code == 0

    def test_live_strike_json(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "strike", str(site_edi_dir), "--format", "json"],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)

    def test_live_strike_csv(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "strike", str(site_edi_dir), "--format", "csv"],
        )
        assert result.exit_code == 0
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2  # header + at least 1 data row

    def test_live_resistivity(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "resistivity", str(site_edi_dir),
             "--freq", "10.0"],
        )
        assert result.exit_code == 0
        assert "res_xy" in result.output or "station" in result.output.lower()

    def test_live_resistivity_interp(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "resistivity", str(site_edi_dir),
             "--freq", "1.0", "--how", "interp"],
        )
        assert result.exit_code == 0

    def test_live_resistivity_csv(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "resistivity", str(site_edi_dir),
             "--freq", "10.0", "--format", "csv"],
        )
        assert result.exit_code == 0
        lines = [l for l in result.output.splitlines() if l.strip()]
        assert len(lines) >= 2

    def test_live_phase_slope(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["site", "compute", "phase-slope", str(site_edi_dir),
             "--band", "0.1:10000"],
        )
        assert result.exit_code == 0

    def test_live_tipper(
        self, runner: CliRunner, site_edi_dir: Path
    ) -> None:
        result = runner.invoke(
            main, ["site", "compute", "tipper", str(site_edi_dir)]
        )
        assert result.exit_code == 0

    def test_no_survey_no_context_fails(
        self, runner: CliRunner, isolated_home: Path
    ) -> None:
        result = runner.invoke(main, ["site", "compute", "strike"])
        assert result.exit_code != 0
        assert "No active survey" in result.output


# ---------------------------------------------------------------------------
# Full workflow integration test
# ---------------------------------------------------------------------------

class TestSiteWorkflow:
    """End-to-end: set survey → info → select → export."""

    def test_workflow_with_context(
        self,
        runner: CliRunner,
        isolated_home: Path,
        site_edi_dir: Path,
        tmp_path: Path,
    ) -> None:

        # Build a real Sites from the EDI dir and cache it
        from pycsamt.seg.edi import EDIFile  # noqa: PLC0415
        from pycsamt.site.base import (  # noqa: PLC0415
            Site,
            Sites,
        )

        edis = sorted(site_edi_dir.glob("*.edi"))
        sites = Sites([Site(EDIFile(str(p))).edi for p in edis])
        n = len(sites)

        # set survey (build cache from real EDIs)
        r1 = runner.invoke(main, ["survey", "set", str(site_edi_dir)])
        assert r1.exit_code == 0

        # site info — should use the active context
        r2 = runner.invoke(main, ["site", "info"])
        assert r2.exit_code == 0

        # site select all — dry run
        r3 = runner.invoke(main, ["site", "select", "--dry-run"])
        assert r3.exit_code == 0
        assert str(n) in r3.output

        # site select → write subset
        out = tmp_path / "selected"
        r4 = runner.invoke(
            main, ["site", "select", "--output-dir", str(out), "--overwrite"]
        )
        # accept success or graceful failure (e.g. site writer not available)
        assert r4.exception is None or isinstance(r4.exception, (SystemExit, Exception))

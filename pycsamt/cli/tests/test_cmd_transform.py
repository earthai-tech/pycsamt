# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt transform`` command group.

Test strategy
-------------
* **Help tests** — always run; verify Click wiring and option names.
* **Unit tests** — use ``SpectraToEDI`` directly (no CLI overhead)
  to verify the transformer class itself.
* **Live integration tests** — require ``data/AMT/SPECTRA/`` (spectra EDIs),
  ``data/avg/`` (AVG files), and Jones J files from ``edi_out/`` or
  ``data/``.  All live tests skip gracefully when data is absent.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main

# ---------------------------------------------------------------------------
# Project-level data paths (resolved from conftest._PROJECT_ROOT)
# ---------------------------------------------------------------------------

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_SPECTRA_DIR = _PROJECT_ROOT / "data" / "AMT" / "SPECTRA"
_AVG_DIR = _PROJECT_ROOT / "data" / "avg"
_J_DIR = _PROJECT_ROOT / "edi_out"  # has test_j_kb0-s001.edi


def _has_spectra() -> bool:
    return _SPECTRA_DIR.exists() and bool(list(_SPECTRA_DIR.glob("*.edi")))


def _has_avg() -> bool:
    return _AVG_DIR.exists() and bool(list(_AVG_DIR.glob("*.avg")))


def _has_j() -> bool:
    return _J_DIR.exists() and bool(list(_J_DIR.glob("*.edi")))


# ---------------------------------------------------------------------------
# pycsamt transform  (group help)
# ---------------------------------------------------------------------------


class TestTransformGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["transform", "--help"])
        assert result.exit_code == 0
        for sub in ("spectra", "avg", "j"):
            assert sub in result.output

    @pytest.mark.parametrize("sub", ["spectra", "avg", "j"])
    def test_each_subcommand_help(self, runner: CliRunner, sub: str) -> None:
        result = runner.invoke(main, ["transform", sub, "--help"])
        assert result.exit_code == 0
        assert "FILE_OR_DIR" in result.output
        assert "--output-dir" in result.output


# ---------------------------------------------------------------------------
# SpectraToEDI unit tests (no CLI, no CliRunner needed)
# ---------------------------------------------------------------------------


class TestSpectraToEDIUnit:
    """Unit tests that exercise the transformer class directly."""

    @pytest.fixture(autouse=True)
    def require_spectra(self) -> None:
        if not _has_spectra():
            pytest.skip(
                "data/AMT/SPECTRA/ not found — skipping spectra unit tests"
            )

    def test_single_file_returns_collection(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI().transform(_SPECTRA_DIR / "HBH03.edi")
        assert len(col) == 1

    def test_single_file_has_z(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI().transform(_SPECTRA_DIR / "HBH03.edi")
        ed = col[0]
        assert ed.Z is not None
        assert ed.Z.n_freq > 0

    def test_single_file_has_tipper(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI().transform(_SPECTRA_DIR / "HBH03.edi")
        ed = col[0]
        assert ed.has_tipper

    def test_directory_processes_all_files(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI().transform(_SPECTRA_DIR)
        n_edi = len(list(_SPECTRA_DIR.glob("*.edi")))
        assert len(col) == n_edi

    def test_station_suffix_applied(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI(station_suffix="_IMP").transform(
            _SPECTRA_DIR / "HBH03.edi"
        )
        assert col[0].station.endswith("_IMP")

    def test_station_name_override(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI().transform(
            _SPECTRA_DIR / "HBH03.edi",
            station_name="CUSTOM_01",
        )
        assert col[0].station == "CUSTOM_01"

    def test_estimate_error(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI(estimate_error=True).transform(
            _SPECTRA_DIR / "HBH03.edi"
        )
        assert col[0].Z.n_freq > 0

    def test_write_to_output_dir(self, tmp_path: Path) -> None:
        from pycsamt.transformers import SpectraToEDI

        SpectraToEDI().transform(
            _SPECTRA_DIR / "HBH03.edi",
            output_dir=tmp_path,
        )
        written = list(tmp_path.glob("*.edi"))
        assert len(written) >= 1

    def test_batch_result_contains_failures(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        t = SpectraToEDI(skip_errors=True)
        # Provide a corrupt file that will fail
        with tempfile.NamedTemporaryFile(suffix=".edi", delete=False) as fp:
            fp.write(b">HEAD\n>END\n")
            bad = Path(fp.name)
        try:
            result = t.transform_batch(bad)
            # Should fail gracefully with 1 failure, 0 ok
            assert result.n_fail == 1
            assert result.n_ok == 0
        finally:
            bad.unlink(missing_ok=True)

    def test_skip_errors_false_raises(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        t = SpectraToEDI(skip_errors=False)
        with tempfile.NamedTemporaryFile(suffix=".edi", delete=False) as fp:
            fp.write(b">HEAD\n>END\n")
            bad = Path(fp.name)
        try:
            with pytest.raises(
                RuntimeError, match="[Cc]onversion failed|failed to convert"
            ):
                t.transform(bad)
        finally:
            bad.unlink(missing_ok=True)

    def test_transform_result_repr(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        result = SpectraToEDI().transform_batch(_SPECTRA_DIR)
        assert "TransformResult" in repr(result)
        assert str(result.n_ok) in repr(result)

    def test_resolve_list_input(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        files = list(_SPECTRA_DIR.glob("*.edi"))
        col = SpectraToEDI().transform(files)
        assert len(col) == len(files)

    def test_with_errors_classmethod(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        t = SpectraToEDI.with_errors()
        assert t.estimate_error is True

    def test_with_tipper_suffix_classmethod(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        t = SpectraToEDI.with_tipper_suffix()
        assert t.station_suffix == "_IMP"

    def test_invalid_source_type_raises(self) -> None:
        from pycsamt.transformers import SpectraToEDI

        with pytest.raises(TypeError):
            SpectraToEDI()._resolve_sources(12345)

    def test_z_shape(self) -> None:
        import numpy as np

        from pycsamt.transformers import SpectraToEDI

        col = SpectraToEDI().transform(_SPECTRA_DIR / "HBH03.edi")
        z = col[0].Z.z
        assert z.ndim == 3
        assert z.shape[1:] == (2, 2)
        assert np.any(np.isfinite(z))


# ---------------------------------------------------------------------------
# CLI integration tests — pycsamt transform spectra
# ---------------------------------------------------------------------------


class TestTransformSpectraCLI:
    @pytest.fixture(autouse=True)
    def require_spectra(self) -> None:
        if not _has_spectra():
            pytest.skip("data/AMT/SPECTRA/ not found")

    def test_dry_run_no_output(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main,
            ["transform", "spectra", str(_SPECTRA_DIR), "--dry-run"],
        )
        assert result.exit_code == 0
        # Nothing written
        assert not list(tmp_path.glob("*.edi"))

    def test_dry_run_lists_files(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["transform", "spectra", str(_SPECTRA_DIR), "--dry-run"],
        )
        n = len(list(_SPECTRA_DIR.glob("*.edi")))
        assert result.exit_code == 0
        assert f"{n} file(s)" in result.output

    def test_single_file_text_output(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = sorted(_SPECTRA_DIR.glob("*.edi"))[0]
        result = runner.invoke(
            main,
            ["transform", "spectra", str(src), "--output-dir", str(tmp_path)],
        )
        assert result.exit_code == 0
        assert "Converted: 1/1" in result.output

    def test_single_file_json_output(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = sorted(_SPECTRA_DIR.glob("*.edi"))[0]
        result = runner.invoke(
            main,
            [
                "transform",
                "spectra",
                str(src),
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["n_ok"] == 1
        assert data["n_fail"] == 0
        assert data["converted"][0]["has_tipper"] is True

    def test_directory_converts_all(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        n = len(list(_SPECTRA_DIR.glob("*.edi")))
        result = runner.invoke(
            main,
            [
                "transform",
                "spectra",
                str(_SPECTRA_DIR),
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0
        assert f"Converted: {n}/{n}" in result.output

    def test_station_suffix_in_output(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        src = sorted(_SPECTRA_DIR.glob("*.edi"))[0]
        result = runner.invoke(
            main,
            [
                "transform",
                "spectra",
                str(src),
                "--station-suffix",
                "_IMP",
                "--output-dir",
                str(tmp_path),
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["converted"][0]["station"].endswith("_IMP")

    def test_estimate_error_flag(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "transform",
                "spectra",
                str(_SPECTRA_DIR),
                "--estimate-error",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code == 0

    def test_verbose_output(self, runner: CliRunner, tmp_path: Path) -> None:
        result = runner.invoke(
            main,
            [
                "transform",
                "spectra",
                str(_SPECTRA_DIR),
                "--output-dir",
                str(tmp_path),
                "-v",
            ],
        )
        assert result.exit_code == 0
        assert "freq" in result.output  # verbose line shows freq count

    def test_missing_output_dir_fails(self, runner: CliRunner) -> None:
        src = sorted(_SPECTRA_DIR.glob("*.edi"))[0]
        result = runner.invoke(main, ["transform", "spectra", str(src)])
        assert result.exit_code != 0
        assert "--output-dir" in result.output

    def test_nonexistent_source_fails(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        result = runner.invoke(
            main,
            [
                "transform",
                "spectra",
                "/no/such/path.edi",
                "--output-dir",
                str(tmp_path),
            ],
        )
        assert result.exit_code != 0

    def test_files_actually_written(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        n = len(list(_SPECTRA_DIR.glob("*.edi")))
        runner.invoke(
            main,
            [
                "transform",
                "spectra",
                str(_SPECTRA_DIR),
                "--output-dir",
                str(tmp_path),
            ],
        )
        written = list(tmp_path.glob("*.edi"))
        assert len(written) == n


# ---------------------------------------------------------------------------
# CLI integration tests — pycsamt transform avg
# ---------------------------------------------------------------------------


class TestTransformAvgCLI:
    @pytest.fixture(autouse=True)
    def require_avg(self) -> None:
        if not _has_avg():
            pytest.skip("data/avg/ not found")

    def test_dry_run(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["transform", "avg", str(_AVG_DIR), "--dry-run"],
        )
        assert result.exit_code == 0

    def test_single_file(self, runner: CliRunner, tmp_path: Path) -> None:
        avgs = sorted(_AVG_DIR.glob("*.avg"))
        result = runner.invoke(
            main,
            ["transform", "avg", str(avgs[0]), "--output-dir", str(tmp_path)],
        )
        assert result.exception is None or isinstance(
            result.exception, SystemExit
        )

    def test_no_output_dir_fails(self, runner: CliRunner) -> None:
        avgs = sorted(_AVG_DIR.glob("*.avg"))
        result = runner.invoke(main, ["transform", "avg", str(avgs[0])])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# CLI integration tests — pycsamt transform j
# ---------------------------------------------------------------------------


class TestTransformJCLI:
    @pytest.fixture(autouse=True)
    def require_j(self) -> None:
        # Look for .j files in data/
        j_files = list((_PROJECT_ROOT / "data").rglob("*.j"))
        if not j_files:
            pytest.skip("No .j files found — skipping J transform tests")
        self.j_dir = j_files[0].parent

    def test_dry_run(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["transform", "j", str(self.j_dir), "--dry-run"],
        )
        assert result.exit_code == 0

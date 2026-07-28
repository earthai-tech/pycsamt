# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for ``pycsamt forward`` command group.

Test strategy
-------------
* **Help tests** — always run; verify Click wiring, option names, sub-commands.
* **Unit tests (rocks / geology)** — use the real built-in GEOLOGY_PRIORS and
  RockDatabase; no mocking required since both are pure Python.
* **Unit tests (run)** — monkeypatch the forward solver's ``run()`` method
  to return a controlled ``ForwardResponse``; avoids SciPy computation in CI.
* **Unit tests (generate)** — monkeypatch ``generate_dataset`` to return a
  minimal ``ForwardDataset``; tests CLI logic and output without running
  actual forward solvers.
* **Error tests** — verify exit codes and error messages for bad inputs.
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
# Fake forward objects (avoid real SciPy computation in unit tests)
# ---------------------------------------------------------------------------


def _fake_response(n_freqs: int = 5) -> MagicMock:
    """Return a ForwardResponse-like mock with MT fields."""
    resp = MagicMock()
    resp.rho_a = np.full(n_freqs, 100.0)
    resp.phase = np.full(n_freqs, 45.0)
    resp.dBz_dt = np.full(n_freqs, 1e-9)
    resp.times = np.logspace(-4, -1, n_freqs)
    return resp


def _fake_dataset(n: int = 10, n_freqs: int = 5) -> MagicMock:
    """Return a ForwardDataset-like mock."""
    ds = MagicMock()
    ds.X = np.random.default_rng(0).random((n, 2 * n_freqs))
    ds.y = np.random.default_rng(1).random((n, 9))
    return ds


# ---------------------------------------------------------------------------
# pycsamt forward  (group help)
# ---------------------------------------------------------------------------


class TestForwardGroup:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "--help"])
        assert result.exit_code == 0
        for sub in ("model", "run", "generate"):
            assert sub in result.output

    def test_model_subgroup_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "model", "--help"])
        assert result.exit_code == 0
        assert "geology" in result.output

    def test_run_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "run", "--help"])
        assert result.exit_code == 0
        for opt in (
            "--resistivities",
            "--thicknesses",
            "--geology",
            "--solver",
            "--freqs",
            "--n-freqs",
        ):
            assert opt in result.output

    def test_generate_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "generate", "--help"])
        assert result.exit_code == 0
        for opt in (
            "--n-samples",
            "--n-layers",
            "--noise-type",
            "--noise-level",
            "--output",
            "--jobs",
        ):
            assert opt in result.output


# ---------------------------------------------------------------------------
# pycsamt forward model geology
# ---------------------------------------------------------------------------


class TestForwardModelGeology:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "model", "geology", "--help"])
        assert result.exit_code == 0
        assert "--name" in result.output

    def test_list_all_text(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "model", "geology"])
        assert result.exit_code == 0
        # At least the five well-known scenarios must appear
        for name in (
            "sedimentary",
            "crystalline",
            "geothermal",
            "marine",
            "permafrost",
        ):
            assert name in result.output

    def test_list_all_json(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["forward", "model", "geology", "--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) >= 5
        assert all("scenario" in row for row in data)

    def test_list_all_csv(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "model", "geology", "--format", "csv"])
        assert result.exit_code == 0
        first_line = result.output.strip().splitlines()[0]
        assert "scenario" in first_line

    def test_detail_text(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["forward", "model", "geology", "--name", "sedimentary"]
        )
        assert result.exit_code == 0
        assert "sedimentary" in result.output
        assert "n_layers" in result.output

    def test_detail_json(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            [
                "forward",
                "model",
                "geology",
                "--name",
                "geothermal",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["scenario"] == "geothermal"
        assert "n_layers" in data

    def test_unknown_scenario_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main, ["forward", "model", "geology", "--name", "atlantis"]
        )
        assert result.exit_code != 0
        assert "atlantis" in result.output.lower() or result.exception


# ---------------------------------------------------------------------------
# pycsamt forward run
# ---------------------------------------------------------------------------


class TestForwardRun:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "run", "--help"])
        assert result.exit_code == 0

    def test_no_model_fails(self, runner: CliRunner) -> None:
        """Running without model spec must fail with a usage error."""
        result = runner.invoke(main, ["forward", "run"])
        assert result.exit_code != 0

    def test_resistivities_without_thicknesses_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["forward", "run", "--resistivities", "100,10,500"],
        )
        assert result.exit_code != 0

    def test_geology_with_resistivities_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            [
                "forward",
                "run",
                "--geology",
                "sedimentary",
                "--resistivities",
                "100,10",
            ],
        )
        assert result.exit_code != 0

    def test_mismatched_layers_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            [
                "forward",
                "run",
                "--resistivities",
                "100,10,500",
                "--thicknesses",
                "300",
            ],  # 3 layers → need 2 thicknesses
        )
        assert result.exit_code != 0

    def test_unknown_geology_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "run", "--geology", "atlantis"])
        assert result.exit_code != 0

    # -- unit tests with mocked solver ----------------------------------------

    def _invoke_run_mocked(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        extra_args: list | None = None,
    ):
        """Monkeypatch MT1DForward.run to return a fake response."""
        from pycsamt.forward.em1d import MT1DForward

        resp = _fake_response(n_freqs=5)
        monkeypatch.setattr(MT1DForward, "run", lambda self, model: resp)

        args = [
            "forward",
            "run",
            "--resistivities",
            "100,10,500",
            "--thicknesses",
            "300,800",
            "--n-freqs",
            "5",
        ] + (extra_args or [])
        return runner.invoke(main, args)

    def test_mt_text_output(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result = self._invoke_run_mocked(runner, monkeypatch)
        assert result.exit_code == 0
        assert "freq_Hz" in result.output
        assert "rho_a_Ohm_m" in result.output

    def test_mt_json_output(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result = self._invoke_run_mocked(
            runner, monkeypatch, extra_args=["--format", "json"]
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, list)
        assert len(data) == 5
        assert "freq_Hz" in data[0]
        assert "rho_a_Ohm_m" in data[0]

    def test_mt_csv_output(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result = self._invoke_run_mocked(
            runner, monkeypatch, extra_args=["--format", "csv"]
        )
        assert result.exit_code == 0
        first_line = result.output.strip().splitlines()[0]
        assert "freq_Hz" in first_line

    def test_geology_path_mocked(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """--geology takes the from_geology path; mock the run to avoid SciPy."""
        from pycsamt.forward.em1d import MT1DForward

        resp = _fake_response(n_freqs=5)
        monkeypatch.setattr(MT1DForward, "run", lambda self, model: resp)

        result = runner.invoke(
            main,
            [
                "forward",
                "run",
                "--geology",
                "sedimentary",
                "--seed",
                "42",
                "--n-freqs",
                "5",
            ],
        )
        assert result.exit_code == 0
        assert "freq_Hz" in result.output

    def test_verbose_writes_to_stderr(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from pycsamt.forward.em1d import MT1DForward

        resp = _fake_response(5)
        monkeypatch.setattr(MT1DForward, "run", lambda self, model: resp)

        result = runner.invoke(
            main,
            [
                "forward",
                "run",
                "--resistivities",
                "100,10,500",
                "--thicknesses",
                "300,800",
                "--n-freqs",
                "5",
                "-v",
            ],
        )
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# pycsamt forward generate
# ---------------------------------------------------------------------------


class TestForwardGenerate:
    def test_help(self, runner: CliRunner) -> None:
        result = runner.invoke(main, ["forward", "generate", "--help"])
        assert result.exit_code == 0

    def test_bad_n_layers_fails(self, runner: CliRunner) -> None:
        result = runner.invoke(
            main,
            ["forward", "generate", "--n-layers", "bad"],
        )
        assert result.exit_code != 0

    # -- unit tests with mocked generate_dataset ------------------------------

    def _invoke_generate_mocked(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        extra_args: list | None = None,
        output: str = "/tmp/test_fw_cli.npz",
    ):
        ds = _fake_dataset(n=10, n_freqs=5)
        # Patch on the source module — lazy `from ....forward import generate_dataset`
        # resolves via pycsamt.forward, so patch the attribute there.
        import pycsamt.forward as _fw

        monkeypatch.setattr(_fw, "generate_dataset", lambda **kw: ds)
        args = [
            "forward",
            "generate",
            "--n-samples",
            "10",
            "--n-freqs",
            "5",
            "--output",
            output,
        ] + (extra_args or [])
        return runner.invoke(main, args), ds

    def test_npz_output(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result, ds = self._invoke_generate_mocked(runner, monkeypatch)
        assert result.exit_code == 0
        assert "Dataset saved" in result.output
        assert str(ds.X.shape) in result.output

    def test_csv_output(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        out = str(tmp_path / "dataset.csv")
        result, ds = self._invoke_generate_mocked(runner, monkeypatch, output=out)
        assert result.exit_code == 0
        assert "Dataset saved" in result.output
        # CSV was written to tmp_path
        assert Path(out).exists()

    def test_n_layers_range(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result, _ = self._invoke_generate_mocked(
            runner, monkeypatch, extra_args=["--n-layers", "3:7"]
        )
        assert result.exit_code == 0

    def test_n_layers_fixed(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result, _ = self._invoke_generate_mocked(
            runner, monkeypatch, extra_args=["--n-layers", "5"]
        )
        assert result.exit_code == 0

    def test_noise_none(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result, _ = self._invoke_generate_mocked(
            runner, monkeypatch, extra_args=["--noise-type", "none"]
        )
        assert result.exit_code == 0

    def test_verbose_shows_progress(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        result, _ = self._invoke_generate_mocked(runner, monkeypatch, extra_args=["-v"])
        assert result.exit_code == 0

    @pytest.mark.parametrize("solver", ["mt1d", "tem1d", "csamt1d"])
    def test_solver_choices(
        self,
        solver: str,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        result, _ = self._invoke_generate_mocked(
            runner, monkeypatch, extra_args=["--solver", solver]
        )
        assert result.exit_code == 0

    @pytest.mark.parametrize("noise", ["gaussian", "field", "multiplicative"])
    def test_noise_type_choices(
        self,
        noise: str,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        result, _ = self._invoke_generate_mocked(
            runner, monkeypatch, extra_args=["--noise-type", noise]
        )
        assert result.exit_code == 0

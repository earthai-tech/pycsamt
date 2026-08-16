# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the ``pycsamt pipe`` command group.

Test strategy
-------------
* **Help tests** — always run; verify Click wiring, option names, sub-commands.
* **Unit tests** — mock ``_resolve_sites`` / ``Pipeline.run`` so no EDI I/O or
  real MT processing is required.  Uses ``make_fake_sites`` + ``_FakeSites``
  from conftest.
* **Param tests** — validate the ``PipeStepList`` Click param type in
  isolation before exercising it through the CLI.
* **Integration tests** — use real WILLY EDI data; skipped automatically
  when ``data/AMT/WILLY_DATA/`` is absent.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.tests.conftest import (
    _FakeSites,
    make_fake_sites,
)

# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

_PIPE_MOD = "pycsamt.cli.commands.pipe"
_RUN_MOD = f"{_PIPE_MOD}.run"

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_DATA_WILLY = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"


def _patch_resolve_sites(monkeypatch: pytest.MonkeyPatch, sites: _FakeSites) -> None:
    """Replace ``_resolve_sites`` in the run sub-module with a lambda."""
    monkeypatch.setattr(
        f"{_RUN_MOD}._resolve_sites",
        lambda *a, **kw: sites,
    )


def _make_fake_result(fake_sites: _FakeSites, outdir: Path | None = None):
    """Build a minimal :class:`~pycsamt.pipeline.PipelineResult` for mocking."""
    from pycsamt.pipeline import (  # noqa: PLC0415
        PipelineResult,
        StepResult,
    )

    sr = StepResult(
        step_idx=1,
        step_name="notch",
        step_code="NR001",
        step_label="Power-line Harmonic Notch",
        params={"mains_hz": 50},
        elapsed_sec=0.05,
        plots=[],
        n_sites_in=len(fake_sites),
        n_sites_out=len(fake_sites),
    )
    return PipelineResult(
        sites_in=fake_sites,
        sites_out=fake_sites,
        step_results=[sr],
        outdir=outdir,
        elapsed_sec=0.08,
        processed_paths=[],
        pipeline_name="mock_pipe",
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def fake_sites() -> _FakeSites:
    return make_fake_sites(4)


@pytest.fixture()
def patched_run(monkeypatch: pytest.MonkeyPatch, fake_sites: _FakeSites):
    """Patch _resolve_sites + Pipeline.run for unit-testing pipe run.

    Returns the fake PipelineResult that Pipeline.run will return.
    """
    from pycsamt.pipeline import Pipeline  # noqa: PLC0415

    _patch_resolve_sites(monkeypatch, fake_sites)
    result = _make_fake_result(fake_sites)
    monkeypatch.setattr(Pipeline, "run", lambda self, *a, **kw: result)
    return result


@pytest.fixture()
def yaml_config_file(tmp_path: Path) -> Path:
    """Write a minimal YAML pipeline config to a temp file."""
    from pycsamt.pipeline import Pipeline  # noqa: PLC0415

    pipe = Pipeline.from_preset("basic_qc")
    p = tmp_path / "workflow.yaml"
    pipe.to_yaml(p)
    return p


# ============================================================================
# 1 · pycsamt pipe  (group help)
# ============================================================================


class TestPipeGroup:
    def test_help_lists_subcommands(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "--help"])
        assert r.exit_code == 0
        for sub in ("run", "steps", "presets", "plugins", "init", "show", "history"):
            assert sub in r.output

    def test_help_lists_opt_in_flags(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "--help"])
        assert r.exit_code == 0
        assert "--with-plugins" in r.output
        assert "--with-ai-steps" in r.output

    def test_each_subcommand_has_help(self, runner: CliRunner) -> None:
        for sub in ("run", "steps", "presets", "plugins", "init", "show", "history"):
            r = runner.invoke(main, ["pipe", sub, "--help"])
            assert r.exit_code == 0, f"pipe {sub} --help failed: {r.output}"
            assert "--help" in r.output

    def test_run_help_lists_key_options(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "run", "--help"])
        for opt in (
            "--config",
            "--preset",
            "--steps",
            "--dry-run",
            "--from-step",
            "--until-step",
            "--n-steps",
            "--on-error",
            "--no-plots",
            "--no-edi",
            "--no-report",
            "--dpi",
            "--plot-fmt",
            "--out",
        ):
            assert opt in r.output, f"Missing option {opt!r} in pipe run --help"


# ============================================================================
# 2 · PipeStepList param type
# ============================================================================


class TestPipeStepList:
    def test_valid_codes_parsed(self) -> None:
        from pycsamt.api.cli.params import PipeStepList

        result = PipeStepList().convert("NR001,FREQ002,FREQ004", None, None)
        assert result == ["NR001", "FREQ002", "FREQ004"]

    def test_valid_names_resolved_to_codes(self) -> None:
        from pycsamt.api.cli.params import PipeStepList

        result = PipeStepList().convert("notch_powerline,drop_duplicates", None, None)
        assert result == ["NR001", "FREQ002"]

    def test_mixed_codes_and_names(self) -> None:
        from pycsamt.api.cli.params import PipeStepList

        result = PipeStepList().convert("NR001,align_grid,SS001", None, None)
        assert result == ["NR001", "FREQ004", "SS001"]

    def test_unknown_code_raises_bad_param(self) -> None:
        from click import BadParameter

        from pycsamt.api.cli.params import PipeStepList

        with pytest.raises(BadParameter, match="Unknown step"):
            PipeStepList().convert("NR001,DOES_NOT_EXIST", None, None)

    def test_empty_string_raises(self) -> None:
        from click import BadParameter

        from pycsamt.api.cli.params import PipeStepList

        with pytest.raises(BadParameter):
            PipeStepList().convert("", None, None)

    def test_unknown_step_surfaces_in_cli(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--steps",
                "BADCODE999",
                "--dry-run",
                "--survey",
                str(_DATA_WILLY),
            ],
        )
        assert "Unknown step" in r.output or r.exit_code != 0


# ============================================================================
# 3 · pycsamt pipe steps
# ============================================================================


class TestPipeSteps:
    def test_all_steps_text(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps"])
        assert r.exit_code == 0
        # All 8 category headers should appear
        for cat in (
            "FREQUENCY",
            "NOISE_REMOVAL",
            "STATIC_SHIFT",
            "TENSOR",
            "DIMENSIONALITY",
            "SKEW",
            "SOURCE_EFFECTS",
            "QC",
        ):
            assert cat in r.output.upper(), f"Category {cat} missing from steps output"

    def test_category_filter_noise_removal(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--category", "noise_removal"])
        assert r.exit_code == 0
        for code in ("NR001", "NR004", "NR010"):
            assert code in r.output
        # Other categories must not appear
        assert "SS001" not in r.output
        assert "TZ001" not in r.output

    def test_info_by_code(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--info", "NR001"])
        assert r.exit_code == 0
        assert "NR001" in r.output
        assert "notch_powerline" in r.output
        assert "nr_qc_harmonic_waterfall" in r.output

    def test_info_by_name(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--info", "notch_powerline"])
        assert r.exit_code == 0
        assert "NR001" in r.output

    def test_info_unknown_raises_error(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--info", "UNKNOWN_CODE"])
        assert r.exit_code != 0
        assert "Unknown step" in r.output or "Error" in r.output

    def test_codes_only_text(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--codes-only"])
        assert r.exit_code == 0
        codes = [line.strip() for line in r.output.splitlines() if line.strip()]
        assert "NR001" in codes
        assert "SS001" in codes
        assert len(codes) == 55

    def test_codes_only_json(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--codes-only", "--format", "json"])
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert isinstance(data, list)
        assert "NR001" in data
        assert len(data) == 55

    def test_format_json_full(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "steps",
                "--category",
                "static_shift",
                "--format",
                "json",
            ],
        )
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert "static_shift" in data
        codes = [s["code"] for s in data["static_shift"]]
        assert "SS001" in codes
        assert "SS002" in codes
        assert "SS003" in codes

    def test_format_csv(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "steps", "--format", "csv"])
        assert r.exit_code == 0
        lines = r.output.strip().splitlines()
        assert lines[0] == "code,name,label,category,returns_sites"
        assert any(l.startswith("NR001,") for l in lines)


# ============================================================================
# 3b · pycsamt pipe plugins
# ============================================================================


@dataclasses.dataclass
class _FakeEntryPoint:
    name: str
    _target: object

    def load(self):
        return self._target


def _register_demo_cli_plugin() -> None:
    from pycsamt.pipeline import StepSpec, register_step  # noqa: PLC0415

    register_step(
        StepSpec(
            code="DEMOCLI001",
            name="demo_cli_step",
            label="Demo CLI Plugin Step",
            category="plugin_demo",
            override_fn=lambda sites, **kw: sites,
        )
    )


def _broken_cli_plugin() -> None:
    raise RuntimeError("deliberately broken CLI plugin")


@pytest.fixture()
def fake_entry_points(monkeypatch: pytest.MonkeyPatch):
    """Patch pycsamt.pipeline discovery to yield a given set of fake plugins,
    and unregister DEMOCLI001 afterwards regardless of what a test did."""
    from pycsamt.pipeline import _plugins as pipeline_plugins_mod
    from pycsamt.pipeline import unregister_step

    def _install(*entry_points: _FakeEntryPoint) -> None:
        monkeypatch.setattr(
            pipeline_plugins_mod, "_iter_entry_points", lambda group: entry_points
        )

    yield _install
    unregister_step("DEMOCLI001", missing_ok=True)


class TestPipePlugins:
    def test_no_plugins_installed_text(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points()  # simulate an environment with none installed
        r = runner.invoke(main, ["pipe", "plugins"])
        assert r.exit_code == 0
        assert "No pipeline plugins found." in r.output
        assert "No plugin steps registered." in r.output

    def test_no_plugins_installed_json(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points()
        r = runner.invoke(main, ["pipe", "plugins", "--format", "json"])
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert data == {"plugins": [], "steps": []}

    def test_no_plugins_installed_csv(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points()
        r = runner.invoke(main, ["pipe", "plugins", "--format", "csv"])
        assert r.exit_code == 0
        assert "plugin_name,ok,error" in r.output
        assert "code,name,label,category" in r.output

    def test_discovers_and_lists_a_registered_plugin(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points(_FakeEntryPoint("demo", _register_demo_cli_plugin))

        r = runner.invoke(main, ["pipe", "plugins"])
        assert r.exit_code == 0
        assert "demo" in r.output
        assert "DEMOCLI001" in r.output
        assert "demo_cli_step" in r.output

    def test_discovers_and_lists_a_registered_plugin_json(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points(_FakeEntryPoint("demo", _register_demo_cli_plugin))

        r = runner.invoke(main, ["pipe", "plugins", "--format", "json"])
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert data["plugins"] == [{"name": "demo", "ok": True, "error": None}]
        assert data["steps"][0]["code"] == "DEMOCLI001"

    def test_strict_propagates_plugin_failure(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points(_FakeEntryPoint("broken", _broken_cli_plugin))

        r = runner.invoke(main, ["pipe", "plugins", "--strict"])
        assert r.exit_code != 0

    def test_default_mode_reports_failure_without_crashing(
        self, runner: CliRunner, fake_entry_points
    ) -> None:
        fake_entry_points(_FakeEntryPoint("broken", _broken_cli_plugin))

        r = runner.invoke(main, ["pipe", "plugins"])
        assert r.exit_code == 0
        assert "broken" in r.output
        assert "FAILED" in r.output


# ============================================================================
# 4 · pycsamt pipe presets
# ============================================================================


class TestPipePresets:
    def test_list_all_text(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "presets"])
        assert r.exit_code == 0
        for name in (
            "basic_qc",
            "noise_reduction",
            "full_processing",
            "tensor_analysis",
            "dimensionality_filter",
            "publication_ready",
        ):
            assert name in r.output

    def test_list_all_json(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "presets", "--format", "json"])
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert isinstance(data, list)
        # the preset catalogue can grow; require the documented core set
        assert len(data) >= 6
        names = [p["name"] for p in data]
        assert "full_processing" in names

    def test_list_all_csv(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "presets", "--format", "csv"])
        assert r.exit_code == 0
        lines = r.output.strip().splitlines()
        assert lines[0].startswith("name,")
        assert any("basic_qc" in l for l in lines)

    def test_expand_preset_text(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "presets", "--expand", "basic_qc"])
        assert r.exit_code == 0
        # basic_qc has NR001, FREQ002, FREQ001, FREQ004, QC001
        for code in ("NR001", "FREQ002", "FREQ001", "FREQ004"):
            assert code in r.output

    def test_expand_preset_json(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "presets",
                "--expand",
                "full_processing",
                "--format",
                "json",
            ],
        )
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert data["name"] == "full_processing"
        assert "steps" in data
        codes = [s["code"] for s in data["steps"]]
        assert "NR001" in codes
        assert "SS001" in codes

    def test_expand_unknown_preset_error(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "presets", "--expand", "nonexistent_preset"])
        assert r.exit_code != 0
        assert "Unknown preset" in r.output or "Error" in r.output


# ============================================================================
# 5 · pycsamt pipe init
# ============================================================================


class TestPipeInit:
    def test_print_yaml(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "init", "--print"])
        assert r.exit_code == 0
        assert "name:" in r.output
        assert "steps:" in r.output

    def test_print_json(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "init", "--print", "--format", "json"])
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert "name" in data
        assert "steps" in data

    def test_print_py(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "init", "--print", "--format", "py"])
        assert r.exit_code == 0
        assert "pipeline_config" in r.output
        # Must be valid Python
        import ast

        ast.parse(r.output)

    def test_print_with_preset(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main, ["pipe", "init", "--print", "--preset", "tensor_analysis"]
        )
        assert r.exit_code == 0
        # Tensor preset codes should be active (uncommented)
        assert "TZ001" in r.output

    def test_write_default_yaml(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "my_workflow.yaml"
        r = runner.invoke(main, ["pipe", "init", "-o", str(out)])
        assert r.exit_code == 0
        assert out.exists()
        assert out.stat().st_size > 0

    def test_write_with_name_and_format(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        out = tmp_path / "willy.py"
        r = runner.invoke(
            main,
            [
                "pipe",
                "init",
                "--name",
                "willy_survey",
                "--format",
                "py",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0
        assert out.exists()
        # generated templates are UTF-8 (contain e.g. Ω·m); the default
        # Windows codec cannot read them back
        assert "willy_survey" in out.read_text(encoding="utf-8")

    def test_write_to_directory_infers_filename(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "init",
                "--name",
                "my_pipe",
                "--format",
                "json",
                "-o",
                str(tmp_path),
            ],
        )
        assert r.exit_code == 0
        expected = tmp_path / "my_pipe.json"
        assert expected.exists()

    def test_yaml_roundtrip(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "wf.yaml"
        r = runner.invoke(
            main,
            [
                "pipe",
                "init",
                "--name",
                "roundtrip_test",
                "--preset",
                "basic_qc",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0
        from pycsamt.pipeline import Pipeline  # noqa: PLC0415

        loaded = Pipeline.from_yaml(out)
        assert loaded.name == "roundtrip_test"
        assert len(loaded) == 5  # basic_qc has 5 steps

    def test_py_roundtrip(self, runner: CliRunner, tmp_path: Path) -> None:
        out = tmp_path / "wf.py"
        r = runner.invoke(
            main,
            [
                "pipe",
                "init",
                "--name",
                "py_test",
                "--format",
                "py",
                "--preset",
                "tensor_analysis",
                "-o",
                str(out),
            ],
        )
        assert r.exit_code == 0
        from pycsamt.pipeline import Pipeline  # noqa: PLC0415

        loaded = Pipeline.from_py(out)
        assert loaded.name == "py_test"
        codes = [s.spec.code for _, s in loaded._steps]
        assert "TZ001" in codes

    def test_invalid_format_raises(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "init", "--print", "--format", "xml"])
        assert r.exit_code != 0


# ============================================================================
# 6 · pycsamt pipe show
# ============================================================================


class TestPipeShow:
    def test_show_preset_text(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "show", "--preset", "tensor_analysis"])
        assert r.exit_code == 0
        for code in ("TZ001", "TZ002", "TZ003", "TZ004"):
            assert code in r.output

    def test_show_preset_json(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            ["pipe", "show", "--preset", "basic_qc", "--format", "json"],
        )
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert data["name"] == "basic_qc"
        assert "steps" in data
        codes = [s["code"] for s in data["steps"]]
        assert "NR001" in codes

    def test_show_preset_csv(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            ["pipe", "show", "--preset", "basic_qc", "--format", "csv"],
        )
        assert r.exit_code == 0
        lines = r.output.strip().splitlines()
        assert lines[0].startswith("idx,")
        assert any("NR001" in l for l in lines)

    def test_show_from_file(self, runner: CliRunner, yaml_config_file: Path) -> None:
        r = runner.invoke(main, ["pipe", "show", str(yaml_config_file)])
        assert r.exit_code == 0
        assert "NR001" in r.output

    def test_show_n_steps_slices(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            ["pipe", "show", "--preset", "full_processing", "--n-steps", "3"],
        )
        assert r.exit_code == 0
        assert "3 steps" in r.output
        # SK001 (step 5) must not appear
        assert "SK001" not in r.output

    def test_show_from_step(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "show",
                "--preset",
                "full_processing",
                "--from-step",
                "mask_skew",
            ],
        )
        assert r.exit_code == 0
        # Steps before mask_skew (NR001, FREQ002, etc.) must not appear
        assert "NR001" not in r.output
        assert "SK001" in r.output

    def test_show_no_source_raises(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "show"])
        assert r.exit_code != 0
        assert "CONFIG_FILE" in r.output or "preset" in r.output.lower()


# ============================================================================
# 6b · pycsamt pipe history
# ============================================================================


class TestPipeHistory:
    def test_no_runs_logged_text(self, runner: CliRunner, tmp_path: Path) -> None:
        r = runner.invoke(
            main, ["pipe", "history", "--file", str(tmp_path / "missing.jsonl")]
        )
        assert r.exit_code == 0
        assert "No pipeline runs logged yet." in r.output

    def test_no_runs_logged_json(self, runner: CliRunner, tmp_path: Path) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "history",
                "--file",
                str(tmp_path / "missing.jsonl"),
                "--format",
                "json",
            ],
        )
        assert r.exit_code == 0
        assert json.loads(r.output) == []

    def test_lists_a_real_logged_run(self, runner: CliRunner, tmp_path: Path) -> None:
        from pycsamt.pipeline import Pipeline

        hist_file = tmp_path / "hist.jsonl"
        pipe = Pipeline.from_preset("basic_qc")
        pipe.name = "cli_history_test"

        class _Sites:
            def __len__(self):
                return 3

            def __iter__(self):
                return iter(())

        # Real run against an empty-but-well-shaped sites stand-in, so this
        # exercises the actual Pipeline.run(history=...) write path rather
        # than mocking it — TestPipeRunUnit already covers the CLI wiring
        # via mocks; this test covers the on-disk round trip specifically.
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pipe.run(
                _Sites(),
                outdir=None,
                save_report=False,
                save_plots=False,
                history=hist_file,
            )

        r = runner.invoke(main, ["pipe", "history", "--file", str(hist_file)])
        assert r.exit_code == 0
        assert "cli_history_test" in r.output

    def test_last_option(self, runner: CliRunner, tmp_path: Path) -> None:
        hist_file = tmp_path / "hist.jsonl"
        hist_file.write_text(
            "\n".join(
                json.dumps(
                    {
                        "timestamp": f"2026-01-0{i}T00:00:00Z",
                        "pipeline_name": f"run{i}",
                        "ok": True,
                        "n_errors": 0,
                        "elapsed_sec": 1.0,
                        "n_sites_in": 1,
                        "n_sites_out": 1,
                        "outdir": None,
                        "steps": [],
                    }
                )
                for i in range(1, 4)
            )
            + "\n",
            encoding="utf-8",
        )
        r = runner.invoke(
            main, ["pipe", "history", "--file", str(hist_file), "--last", "1"]
        )
        assert r.exit_code == 0
        assert "run3" in r.output
        assert "run1" not in r.output

    def test_csv_format(self, runner: CliRunner, tmp_path: Path) -> None:
        hist_file = tmp_path / "hist.jsonl"
        hist_file.write_text(
            json.dumps(
                {
                    "timestamp": "2026-01-01T00:00:00Z",
                    "pipeline_name": "run1",
                    "ok": True,
                    "n_errors": 0,
                    "elapsed_sec": 1.0,
                    "n_sites_in": 1,
                    "n_sites_out": 1,
                    "outdir": None,
                    "steps": [],
                }
            )
            + "\n",
            encoding="utf-8",
        )
        r = runner.invoke(
            main, ["pipe", "history", "--file", str(hist_file), "--format", "csv"]
        )
        assert r.exit_code == 0
        lines = r.output.strip().splitlines()
        assert lines[0].startswith("timestamp,pipeline_name,ok,")
        assert "run1" in lines[1]


# ============================================================================
# 7 · pycsamt pipe run  (unit — mocked)
# ============================================================================


class TestPipeRunUnit:
    # ── help ────────────────────────────────────────────────────────────────

    def test_help(self, runner: CliRunner) -> None:
        r = runner.invoke(main, ["pipe", "run", "--help"])
        assert r.exit_code == 0

    # ── dry-run: no processing ──────────────────────────────────────────────

    def test_dry_run_preset_shows_pipeline(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert "NR001" in r.output
        assert "Dry run" in r.output

    def test_dry_run_steps_builds_pipeline(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--steps",
                "NR001,FREQ002,FREQ004",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert "NR001" in r.output
        assert "FREQ002" in r.output
        assert "FREQ004" in r.output

    def test_steps_flag_rejects_unknown_plugin_code_by_default(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
        fake_entry_points,
    ) -> None:
        """Without --with-plugins, a plugin step code is never discovered —
        this is the opt-in guarantee, and the error should say how to opt in."""
        _patch_resolve_sites(monkeypatch, fake_sites)
        fake_entry_points(_FakeEntryPoint("demo", _register_demo_cli_plugin))
        r = runner.invoke(
            main,
            ["pipe", "run", "--steps", "DEMOCLI001,NR001", "--survey", ".", "--dry-run"],
        )
        assert r.exit_code != 0
        assert "--with-plugins" in r.output

    def test_with_plugins_flag_makes_plugin_step_usable_in_one_shot(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
        fake_entry_points,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        fake_entry_points(_FakeEntryPoint("demo", _register_demo_cli_plugin))
        r = runner.invoke(
            main,
            [
                "pipe",
                "--with-plugins",
                "run",
                "--steps",
                "DEMOCLI001,NR001",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert "DEMOCLI001" in r.output

    def test_export_steps_need_no_opt_in(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--steps",
                "NR001,EXPORT001",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert "EXPORT001" in r.output

    def test_ai_step_rejected_without_with_ai_steps(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            ["pipe", "run", "--steps", "AI001", "--survey", ".", "--dry-run"],
        )
        assert r.exit_code != 0
        assert "Unknown step" in r.output

    def test_with_ai_steps_flag_makes_ai001_usable_in_one_shot(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        from pycsamt.pipeline import unregister_step

        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "--with-ai-steps",
                "run",
                "--steps",
                "AI001",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        unregister_step("AI001", missing_ok=True)
        assert r.exit_code == 0
        assert "AI001" in r.output

    def test_dry_run_shows_site_count(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert str(len(fake_sites)) in r.output

    def test_dry_run_shows_cache_disabled_by_default(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            ["pipe", "run", "--preset", "basic_qc", "--survey", ".", "--dry-run"],
        )
        assert r.exit_code == 0
        assert "Cache   : disabled" in r.output

    def test_dry_run_shows_cache_dir_when_given(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
        tmp_path: Path,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--cache-dir",
                str(tmp_path),
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert str(tmp_path) in r.output

    def test_dry_run_shows_history_disabled_by_default(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            ["pipe", "run", "--preset", "basic_qc", "--survey", ".", "--dry-run"],
        )
        assert r.exit_code == 0
        assert "History : disabled" in r.output

    def test_dry_run_shows_history_file_when_given(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
        tmp_path: Path,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        hist_file = tmp_path / "hist.jsonl"
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--history-file",
                str(hist_file),
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert str(hist_file) in r.output

    def test_dry_run_n_steps_slices(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "full_processing",
                "--survey",
                ".",
                "--dry-run",
                "--n-steps",
                "2",
            ],
        )
        assert r.exit_code == 0
        assert "2 steps" in r.output

    def test_dry_run_from_step(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "full_processing",
                "--survey",
                ".",
                "--dry-run",
                "--from-step",
                "mask_skew",
            ],
        )
        assert r.exit_code == 0
        # NR001 was before mask_skew and should be sliced
        assert "NR001" not in r.output
        assert "SK001" in r.output

    def test_dry_run_until_step(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "full_processing",
                "--survey",
                ".",
                "--dry-run",
                "--until-step",
                "align_grid",
            ],
        )
        assert r.exit_code == 0
        # SS001 (static shift) comes after align_grid and must be sliced out
        assert "SS001" not in r.output

    # ── error paths ──────────────────────────────────────────────────────────

    def test_no_pipeline_specified_raises(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(main, ["pipe", "run", "--survey", "."])
        assert r.exit_code != 0
        assert "--config" in r.output or "--preset" in r.output or "--steps" in r.output

    def test_unknown_preset_raises(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "does_not_exist",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code != 0
        assert "Unknown preset" in r.output or "Error" in r.output

    def test_unknown_step_in_steps_flag(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--steps",
                "NR001,BADSTEP",
                "--survey",
                ".",
                "--dry-run",
            ],
        )
        assert r.exit_code != 0
        assert "Unknown step" in r.output or "Error" in r.output

    def test_bad_from_step_raises(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--dry-run",
                "--from-step",
                "ghost_step",
            ],
        )
        assert r.exit_code != 0

    # ── mocked execution ─────────────────────────────────────────────────────

    def test_run_text_output(self, runner: CliRunner, patched_run) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
        )
        assert r.exit_code == 0
        assert "mock_pipe" in r.output or "sites" in r.output.lower()

    def test_cache_flag_forwarded_to_pipeline_run(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        captured = {}

        def _fake_run(self, *a, **kw):
            captured.update(kw)
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
                "--cache",
            ],
        )
        assert r.exit_code == 0
        assert captured["cache"] is True

    def test_cache_dir_flag_forwarded_to_pipeline_run(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites,
        tmp_path: Path,
    ) -> None:
        from pycsamt.pipeline import Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        captured = {}

        def _fake_run(self, *a, **kw):
            captured.update(kw)
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
                "--cache-dir",
                str(tmp_path),
            ],
        )
        assert r.exit_code == 0
        assert captured["cache"] == tmp_path

    def test_no_cache_flag_forwards_false(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        captured = {}

        def _fake_run(self, *a, **kw):
            captured.update(kw)
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
        )
        assert r.exit_code == 0
        assert captured["cache"] is False

    def test_history_flag_forwarded_to_pipeline_run(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        captured = {}

        def _fake_run(self, *a, **kw):
            captured.update(kw)
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
                "--history",
            ],
        )
        assert r.exit_code == 0
        assert captured["history"] is True

    def test_history_file_flag_forwarded_to_pipeline_run(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites,
        tmp_path: Path,
    ) -> None:
        from pycsamt.pipeline import Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        captured = {}

        def _fake_run(self, *a, **kw):
            captured.update(kw)
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        hist_file = tmp_path / "hist.jsonl"
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
                "--history-file",
                str(hist_file),
            ],
        )
        assert r.exit_code == 0
        assert captured["history"] == hist_file

    def test_no_history_flag_forwards_false(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        captured = {}

        def _fake_run(self, *a, **kw):
            captured.update(kw)
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
        )
        assert r.exit_code == 0
        assert captured["history"] is False

    def test_live_flag_sets_rich_progress_style(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import PYCSAMT_PIPE, Pipeline
        from pycsamt.api.pipe import reset_pipe

        _patch_resolve_sites(monkeypatch, fake_sites)
        seen_style = {}

        def _fake_run(self, *a, **kw):
            seen_style["progress_style"] = PYCSAMT_PIPE.progress_style
            seen_style["show_progress"] = PYCSAMT_PIPE.show_progress
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        try:
            r = runner.invoke(
                main,
                [
                    "pipe",
                    "run",
                    "--preset",
                    "basic_qc",
                    "--survey",
                    ".",
                    "--no-plots",
                    "--no-edi",
                    "--no-report",
                    "--live",
                ],
            )
            assert r.exit_code == 0
            assert seen_style["progress_style"] == "rich"
            assert seen_style["show_progress"] is True
        finally:
            # --live sets progress_style="rich" on the global PYCSAMT_PIPE
            # singleton with no matching reset in run.py itself; clean up
            # explicitly so it can't leak into a later test in this file.
            reset_pipe()

    def test_without_live_flag_progress_style_untouched(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import PYCSAMT_PIPE, Pipeline, configure_pipe
        from pycsamt.api.pipe import reset_pipe

        configure_pipe(progress_style="log")
        try:
            _patch_resolve_sites(monkeypatch, fake_sites)
            seen_style = {}

            def _fake_run(self, *a, **kw):
                seen_style["progress_style"] = PYCSAMT_PIPE.progress_style
                return _make_fake_result(fake_sites)

            monkeypatch.setattr(Pipeline, "run", _fake_run)
            r = runner.invoke(
                main,
                [
                    "pipe",
                    "run",
                    "--preset",
                    "basic_qc",
                    "--survey",
                    ".",
                    "--no-plots",
                    "--no-edi",
                    "--no-report",
                ],
            )
            assert r.exit_code == 0
            # --live was NOT passed: a user's own progress_style must survive.
            assert seen_style["progress_style"] == "log"
        finally:
            reset_pipe()  # this test mutates the global PYCSAMT_PIPE singleton;
            # no autouse reset fixture exists in this file, so clean up
            # explicitly to avoid leaking "log" into later tests' state.

    def test_dashboard_flag_adds_dashboard_report_format(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import PYCSAMT_PIPE, Pipeline
        from pycsamt.api.pipe import reset_pipe

        _patch_resolve_sites(monkeypatch, fake_sites)
        seen = {}

        def _fake_run(self, *a, **kw):
            seen["report_formats"] = PYCSAMT_PIPE.report_formats
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        try:
            r = runner.invoke(
                main,
                [
                    "pipe",
                    "run",
                    "--preset",
                    "basic_qc",
                    "--survey",
                    ".",
                    "--no-plots",
                    "--no-edi",
                    "--no-report",
                    "--dashboard",
                ],
            )
            assert r.exit_code == 0
            assert "dashboard" in seen["report_formats"]
        finally:
            reset_pipe()

    def test_without_dashboard_flag_report_formats_untouched(
        self, runner: CliRunner, monkeypatch: pytest.MonkeyPatch, fake_sites
    ) -> None:
        from pycsamt.pipeline import PYCSAMT_PIPE, Pipeline

        _patch_resolve_sites(monkeypatch, fake_sites)
        seen = {}

        def _fake_run(self, *a, **kw):
            seen["report_formats"] = PYCSAMT_PIPE.report_formats
            return _make_fake_result(fake_sites)

        monkeypatch.setattr(Pipeline, "run", _fake_run)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
        )
        assert r.exit_code == 0
        assert "dashboard" not in seen["report_formats"]

    def test_dry_run_shows_dashboard_disabled_by_default(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            ["pipe", "run", "--preset", "basic_qc", "--survey", ".", "--dry-run"],
        )
        assert r.exit_code == 0
        assert "Dashboard : disabled" in r.output

    def test_dry_run_shows_dashboard_enabled_when_flag_given(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        _patch_resolve_sites(monkeypatch, fake_sites)
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--dashboard",
                "--dry-run",
            ],
        )
        assert r.exit_code == 0
        assert "Dashboard : enabled" in r.output

    def test_run_json_output(self, runner: CliRunner, patched_run) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
                "--format",
                "json",
            ],
        )
        assert r.exit_code == 0
        data = json.loads(r.output)
        assert data["ok"] is True
        assert "steps" in data
        assert data["pipeline"] == "mock_pipe"

    def test_run_csv_output(self, runner: CliRunner, patched_run) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
                "--format",
                "csv",
            ],
        )
        assert r.exit_code == 0
        lines = r.output.strip().splitlines()
        assert lines[0].startswith("idx,")
        assert any("NR001" in l for l in lines)

    def test_run_from_config_file(
        self,
        runner: CliRunner,
        patched_run,
        yaml_config_file: Path,
    ) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--config",
                str(yaml_config_file),
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
        )
        assert r.exit_code == 0

    def test_run_error_pipeline_exits_1(
        self,
        runner: CliRunner,
        monkeypatch: pytest.MonkeyPatch,
        fake_sites: _FakeSites,
    ) -> None:
        """A pipeline with at least one failed step should exit with code 1."""
        from pycsamt.pipeline import (  # noqa: PLC0415
            Pipeline,
            PipelineResult,
            StepResult,
        )

        _patch_resolve_sites(monkeypatch, fake_sites)

        bad_sr = StepResult(
            step_idx=1,
            step_name="notch",
            step_code="NR001",
            step_label="Power-line Harmonic Notch",
            params={},
            elapsed_sec=0.01,
            plots=[],
            n_sites_in=4,
            n_sites_out=4,
            error=RuntimeError("deliberate"),
        )
        bad_result = PipelineResult(
            sites_in=fake_sites,
            sites_out=fake_sites,
            step_results=[bad_sr],
            outdir=None,
            elapsed_sec=0.01,
            processed_paths=[],
            pipeline_name="err_pipe",
        )
        monkeypatch.setattr(Pipeline, "run", lambda self, *a, **kw: bad_result)

        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                ".",
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
        )
        assert r.exit_code == 1


# ============================================================================
# 8 · pycsamt pipe run  (integration — real WILLY data)
# ============================================================================


@pytest.mark.skipif(
    not _DATA_WILLY.exists(),
    reason=f"WILLY EDI data not found at {_DATA_WILLY}",
)
class TestPipeRunIntegration:
    """Integration tests that exercise real processing on 5 WILLY stations."""

    def test_dry_run_real_data(self, runner: CliRunner) -> None:
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                str(_DATA_WILLY),
                "--dry-run",
            ],
        )
        assert r.exit_code == 0, r.output
        assert "NR001" in r.output
        assert "Dry run" in r.output

    def test_run_basic_qc_no_disk_writes(self, runner: CliRunner) -> None:
        """Run basic_qc with all disk output disabled — verifies step execution."""
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                str(_DATA_WILLY),
                "--no-plots",
                "--no-edi",
                "--no-report",
            ],
            # No --out → outdir=None → no files written.  (--out "" would
            # resolve to Path(".") and drop pipeline.yaml in the CWD.)
        )
        # Exit 0 means all 5 steps completed without error
        assert r.exit_code == 0, r.output

    def test_run_writes_pipe_results(self, runner: CliRunner, tmp_path: Path) -> None:
        """End-to-end: run + write processed EDIs and reports."""
        out = tmp_path / "pipe_results_cli_test"
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                str(_DATA_WILLY),
                "--out",
                str(out),
                "--no-plots",  # skip plots for speed
                "--dpi",
                "72",
                "-v",
            ],
        )
        assert r.exit_code == 0, r.output
        assert (out / "processed").is_dir()
        assert list((out / "processed").glob("*.edi")), "No EDIs in processed/"
        assert (out / "pipeline.yaml").exists()
        assert (out / "summary.txt").exists()

    def test_run_with_dashboard_flag_writes_dashboard_html(
        self, runner: CliRunner, tmp_path: Path
    ) -> None:
        """End-to-end: --dashboard writes dashboard.html alongside the fast tier."""
        out = tmp_path / "pipe_results_dashboard_cli_test"
        r = runner.invoke(
            main,
            [
                "pipe",
                "run",
                "--preset",
                "basic_qc",
                "--survey",
                str(_DATA_WILLY),
                "--out",
                str(out),
                "--no-plots",
                "--dpi",
                "72",
                "--dashboard",
            ],
        )
        assert r.exit_code == 0, r.output
        assert (out / "report.html").exists()
        assert (out / "summary.txt").exists()
        dashboard_html = out / "dashboard.html"
        assert dashboard_html.exists()
        content = dashboard_html.read_text(encoding="utf-8")
        assert "pycsamt-dashboard" in content
        assert "NR001" in content

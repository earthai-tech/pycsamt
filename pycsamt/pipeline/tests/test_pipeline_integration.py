"""Integration test: real WILLY AMT data → full pipe_results/ output.

This test loads real EDI files from ``data/AMT/WILLY_DATA/L22PLT/``, runs a
multi-step pipeline, and verifies every expected output artefact is produced:

* ``pipe_results/processed/``    — processed EDI files
* ``pipe_results/plots/``        — QC figures, one sub-folder per step
* ``pipe_results/pipeline.yaml`` — reproduced config (reproducibility)
* ``pipe_results/summary.txt``   — plain-text run report
* ``pipe_results/report.html``   — HTML run report

The output directory is **not cleaned up** after the test so the user can
inspect the results directly.  It is listed in ``.gitignore`` to prevent
accidental commits.

Skip condition: the test is skipped automatically when the ``data/`` tree is
not present (e.g. on a fresh CI clone without the data submodule).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

from pycsamt.emtools._core import ensure_sites
from pycsamt.pipeline import (
    Pipeline,
    PipelineResult,
    Step,
    configure_pipe,
    reset_pipe,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_DATA_DIR = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
_PIPE_RESULTS = _PROJECT_ROOT / "pipe_results"

# ── Skip guard ────────────────────────────────────────────────────────────────

pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"WILLY EDI data not found at {_DATA_DIR}",
)

# ── How many stations to use (keep the test fast) ────────────────────────────

_N_STATIONS = 5

# ── Module-level fixtures (run once for the whole module) ─────────────────────


@pytest.fixture(scope="module")
def willy_sites():
    """Load the first *_N_STATIONS* stations from the WILLY L22 profile."""
    paths = sorted(_DATA_DIR.glob("*.edi"))[:_N_STATIONS]
    assert paths, f"No EDI files found in {_DATA_DIR}"
    return ensure_sites([str(p) for p in paths])


@pytest.fixture(scope="module")
def pipeline_result(willy_sites):
    """Run the pipeline once and yield the result.

    The output directory is cleared before the run so each test execution
    starts with a clean slate.  It is *not* removed afterwards so the user
    can inspect ``pipe_results/``.
    """
    configure_pipe(show_progress=False, plot_dpi=72, plot_fmt="png")

    if _PIPE_RESULTS.exists():
        shutil.rmtree(_PIPE_RESULTS)

    pipe = Pipeline(
        [
            ("drop_dup", Step("FREQ002")),
            ("select_band", Step("FREQ001", band_hz=(0.01, 10_000))),
            ("align", Step("FREQ004")),
            ("notch", Step("NR001", mains_hz=50)),
            ("qc_snap", Step("QC001")),
        ],
        name="willy_l22_integration",
    )

    result = pipe.run(
        willy_sites,
        outdir=_PIPE_RESULTS,
        save_plots=True,
        save_edis=True,
        save_report=True,
    )

    yield result

    reset_pipe()


# ── Tests ─────────────────────────────────────────────────────────────────────


class TestPipelineResult:
    """Validate the PipelineResult object returned by the run."""

    def test_result_type(self, pipeline_result):
        assert isinstance(pipeline_result, PipelineResult)

    def test_result_ok(self, pipeline_result):
        assert pipeline_result.ok, (
            f"Pipeline had {pipeline_result.n_errors} error(s):\n"
            + "\n".join(
                f"  [{sr.step_code}] {sr.step_name}: {sr.error}"
                for sr in pipeline_result.step_results
                if not sr.ok
            )
        )

    def test_n_errors_zero(self, pipeline_result):
        assert pipeline_result.n_errors == 0

    def test_step_count(self, pipeline_result):
        assert len(pipeline_result.step_results) == 5

    def test_all_steps_ok(self, pipeline_result):
        for sr in pipeline_result.step_results:
            assert sr.ok, f"Step {sr.step_code!r} failed: {sr.error}"

    def test_elapsed_positive(self, pipeline_result):
        assert pipeline_result.elapsed_sec > 0

    def test_sites_in_count(self, pipeline_result):
        assert len(pipeline_result.sites_in) == _N_STATIONS

    def test_sites_out_count(self, pipeline_result):
        assert len(pipeline_result.sites_out) == _N_STATIONS

    def test_outdir_is_pipe_results(self, pipeline_result):
        assert pipeline_result.outdir == _PIPE_RESULTS

    def test_summary_is_string(self, pipeline_result):
        s = pipeline_result.summary()
        assert isinstance(s, str)
        assert "willy_l22_integration" in s


class TestOutputDirectory:
    """Verify the on-disk artefacts inside pipe_results/."""

    def test_root_exists(self, pipeline_result):
        assert _PIPE_RESULTS.exists()

    def test_processed_subdir_exists(self, pipeline_result):
        assert (_PIPE_RESULTS / "processed").is_dir()

    def test_plots_subdir_exists(self, pipeline_result):
        assert (_PIPE_RESULTS / "plots").is_dir()

    def test_pipeline_yaml_written(self, pipeline_result):
        assert (_PIPE_RESULTS / "pipeline.yaml").is_file()

    def test_summary_txt_written(self, pipeline_result):
        assert (_PIPE_RESULTS / "summary.txt").is_file()

    def test_report_html_written(self, pipeline_result):
        assert (_PIPE_RESULTS / "report.html").is_file()


class TestProcessedEDIs:
    """Verify EDI files were written to processed/."""

    def test_edi_files_exist(self, pipeline_result):
        edi_files = list((_PIPE_RESULTS / "processed").glob("*.edi"))
        assert len(edi_files) > 0, "No EDI files in processed/"

    def test_edi_count_matches_sites(self, pipeline_result):
        edi_files = list((_PIPE_RESULTS / "processed").glob("*.edi"))
        assert len(edi_files) == _N_STATIONS

    def test_edi_files_are_non_empty(self, pipeline_result):
        for edi in (_PIPE_RESULTS / "processed").glob("*.edi"):
            assert edi.stat().st_size > 0, f"{edi.name} is empty"


class TestPlotOutput:
    """Verify QC figures were saved under plots/."""

    def test_plots_has_step_subdirs(self, pipeline_result):
        subdirs = [
            d for d in (_PIPE_RESULTS / "plots").iterdir() if d.is_dir()
        ]
        assert len(subdirs) > 0, "No step subdirectories in plots/"

    def test_step_subdir_names_are_numbered(self, pipeline_result):
        subdirs = sorted(
            d.name for d in (_PIPE_RESULTS / "plots").iterdir() if d.is_dir()
        )
        for name in subdirs:
            prefix = name.split("_")[0]
            assert prefix.isdigit(), f"Expected numeric prefix, got {name!r}"

    def test_at_least_one_png_produced(self, pipeline_result):
        pngs = list((_PIPE_RESULTS / "plots").rglob("*.png"))
        assert len(pngs) > 0, "No PNG files produced under plots/"

    def test_plots_count_matches_result(self, pipeline_result):
        pngs = list((_PIPE_RESULTS / "plots").rglob("*.png"))
        assert len(pngs) == len(pipeline_result.plots)

    def test_plot_paths_point_inside_outdir(self, pipeline_result):
        for p in pipeline_result.plots:
            assert _PIPE_RESULTS in p.parents, f"{p} is outside pipe_results/"


class TestReportContent:
    """Verify the written reports contain expected content."""

    def test_summary_txt_contains_pipeline_name(self, pipeline_result):
        txt = (_PIPE_RESULTS / "summary.txt").read_text(encoding="utf-8")
        assert "willy_l22_integration" in txt

    def test_summary_txt_contains_step_codes(self, pipeline_result):
        txt = (_PIPE_RESULTS / "summary.txt").read_text(encoding="utf-8")
        for code in ("FREQ002", "FREQ001", "FREQ004", "NR001", "QC001"):
            assert code in txt, f"Step code {code!r} missing from summary.txt"

    def test_summary_txt_shows_ok_status(self, pipeline_result):
        txt = (_PIPE_RESULTS / "summary.txt").read_text(encoding="utf-8")
        assert "OK" in txt

    def test_report_html_is_valid_html(self, pipeline_result):
        html = (_PIPE_RESULTS / "report.html").read_text(encoding="utf-8")
        assert "<html" in html
        assert "</html>" in html

    def test_report_html_contains_pipeline_name(self, pipeline_result):
        html = (_PIPE_RESULTS / "report.html").read_text(encoding="utf-8")
        assert "willy_l22_integration" in html

    def test_report_html_contains_all_step_codes(self, pipeline_result):
        html = (_PIPE_RESULTS / "report.html").read_text(encoding="utf-8")
        for code in ("FREQ002", "NR001", "QC001"):
            assert code in html, (
                f"Step code {code!r} missing from report.html"
            )


class TestPipelineYAML:
    """Verify pipeline.yaml is valid and reloadable."""

    def test_yaml_is_valid_yaml(self, pipeline_result):
        import yaml

        text = (_PIPE_RESULTS / "pipeline.yaml").read_text(encoding="utf-8")
        data = yaml.safe_load(text)
        assert isinstance(data, dict)

    def test_yaml_contains_name(self, pipeline_result):
        import yaml

        data = yaml.safe_load(
            (_PIPE_RESULTS / "pipeline.yaml").read_text(encoding="utf-8")
        )
        assert data.get("name") == "willy_l22_integration"

    def test_yaml_steps_count(self, pipeline_result):
        import yaml

        data = yaml.safe_load(
            (_PIPE_RESULTS / "pipeline.yaml").read_text(encoding="utf-8")
        )
        assert len(data["steps"]) == 5

    def test_yaml_reloadable_as_pipeline(self, pipeline_result):
        reloaded = Pipeline.from_yaml(_PIPE_RESULTS / "pipeline.yaml")
        assert isinstance(reloaded, Pipeline)
        assert len(reloaded) == 5
        assert reloaded.name == "willy_l22_integration"

    def test_reloaded_pipeline_runs_on_same_data(
        self, pipeline_result, willy_sites
    ):
        """A pipeline reloaded from the saved YAML must run without error."""
        reloaded = Pipeline.from_yaml(_PIPE_RESULTS / "pipeline.yaml")
        result2 = reloaded.run(
            willy_sites,
            outdir=None,  # in-memory, no disk writes
            save_plots=False,
            save_edis=False,
            save_report=False,
        )
        assert result2.ok
        assert len(result2.sites_out) == _N_STATIONS


class TestStepResultDetails:
    """Verify per-step metadata recorded in PipelineResult."""

    def test_step_indices_are_sequential(self, pipeline_result):
        indices = [sr.step_idx for sr in pipeline_result.step_results]
        assert indices == list(range(1, 6))

    def test_step_codes_in_order(self, pipeline_result):
        codes = [sr.step_code for sr in pipeline_result.step_results]
        assert codes == ["FREQ002", "FREQ001", "FREQ004", "NR001", "QC001"]

    def test_step_labels_match_registry(self, pipeline_result):
        from pycsamt.pipeline import lookup_step

        for sr in pipeline_result.step_results:
            spec = lookup_step(sr.step_code)
            assert sr.step_label == spec.label

    def test_all_steps_report_elapsed(self, pipeline_result):
        for sr in pipeline_result.step_results:
            assert sr.elapsed_sec >= 0

    def test_sites_counts_non_zero(self, pipeline_result):
        for sr in pipeline_result.step_results:
            assert sr.n_sites_in > 0
            assert sr.n_sites_out > 0

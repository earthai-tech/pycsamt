"""Integration tests: method-aware presets against real MT/AMT/CSAMT data.

Confirms the "smart" gating in the mt_qc/amt_qc/csamt_qc/csumt_qc presets
actually differs by dataset — not just by name — using three real,
already-vendored datasets:

* ``data/MT/kap03lmt_edis``       — real long-period MT, has tipper.
* ``data/AMT/WILLY_DATA/L22PLT``  — real AMT, no tipper.
* ``data/CSAMT``                  — real Tongkeng CSAMT line (published
  dataset), multi-component, no resolvable source_offset in its headers.

Each is skipped individually if its data directory is not present (e.g. a
fresh checkout without the data submodule), matching
test_pipeline_integration.py's convention.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest

from pycsamt.api import read_edis
from pycsamt.emtools._core import _get_z_block, _iter_items
from pycsamt.emtools.source_effects import correct_near_field
from pycsamt.pipeline import Pipeline, configure_pipe, reset_pipe

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_DATA_MT = _PROJECT_ROOT / "data" / "MT" / "kap03lmt_edis"
_DATA_AMT = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
_DATA_CSAMT = _PROJECT_ROOT / "data" / "CSAMT"
_OUTROOT = _PROJECT_ROOT / "pipe_results" / "smart_presets_test"

_N_STATIONS = 5


def _load(data_dir: Path, n: int = _N_STATIONS):
    paths = sorted(data_dir.glob("*.edi"))[:n]
    assert paths, f"No EDI files found in {data_dir}"
    return read_edis(str(data_dir), strict=False).to_collection()


def _run(preset_name: str, data_dir: Path, outdir: Path):
    configure_pipe(show_progress=False, plot_dpi=60, plot_fmt="png")
    if outdir.exists():
        shutil.rmtree(outdir)
    sites = _load(data_dir)
    pipe = Pipeline.from_preset(preset_name)
    result = pipe.run(
        sites, outdir=outdir, save_plots=True, save_edis=False, save_report=False
    )
    reset_pipe()
    return result


def _plot_names(result) -> set[str]:
    return {p.stem for p in result.plots}


# ─────────────────────────────────────────────────────────────────────────────
# mt_qc against real MT data (has tipper)
# ─────────────────────────────────────────────────────────────────────────────


# Each class below has a class-scoped fixture that writes to a fixed
# shared directory; pytest-xdist's default load-balanced scheduling does
# not keep a class's tests on one worker, so xdist_group pins each class
# to a single worker to avoid two workers racing on the same output dir.
@pytest.mark.skipif(not _DATA_MT.exists(), reason=f"MT data not found at {_DATA_MT}")
@pytest.mark.xdist_group(name="presets_mt_qc")
class TestMtQcRealData:
    @pytest.fixture(scope="class")
    def result(self):
        return _run("mt_qc", _DATA_MT, _OUTROOT / "mt_qc")

    def test_run_ok(self, result):
        assert result.ok, "\n".join(
            f"{sr.step_code}: {sr.error}" for sr in result.step_results if not sr.ok
        )

    def test_raw_and_processed_preview_present(self, result):
        names = _plot_names(result)
        assert "plot_raw_preview" in names
        assert "plot_processed_preview" in names

    def test_phase_tensor_smart_fires_for_multicomponent_mt(self, result):
        # Real MT data is full-tensor -- the gated plot should fire.
        assert "phase_tensor_smart" in _plot_names(result)

    def test_tipper_smart_plots_fire_for_real_mt_tipper(self, result):
        # kap03lmt is real long-period MT data with a tipper channel.
        names = _plot_names(result)
        for fn_name in (
            "induction_multiperiod_map_smart",
            "induction_section_smart",
            "response_tipper_smart",
            "tipper_components_smart",
            "tipper_hodograms_smart",
        ):
            assert fn_name in names, f"{fn_name} did not fire for tipper-bearing MT data"

    def test_strike_qc_always_fires(self, result):
        names = _plot_names(result)
        assert "plot_strike_analysis" in names
        assert "plot_strike_rose" in names


# ─────────────────────────────────────────────────────────────────────────────
# amt_qc against real AMT data (no tipper)
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _DATA_AMT.exists(), reason=f"AMT data not found at {_DATA_AMT}"
)
class TestAmtQcRealData:
    @pytest.fixture(scope="class")
    def result(self):
        return _run("amt_qc", _DATA_AMT, _OUTROOT / "amt_qc")

    def test_run_ok(self, result):
        assert result.ok, "\n".join(
            f"{sr.step_code}: {sr.error}" for sr in result.step_results if not sr.ok
        )

    def test_phase_tensor_smart_fires_for_multicomponent_amt(self, result):
        assert "phase_tensor_smart" in _plot_names(result)

    def test_tipper_smart_plots_absent_for_no_tipper_amt(self, result):
        # WILLY AMT data carries no tipper channel -- every tipper-gated
        # plot must be silently skipped, not drawn as an empty placeholder.
        names = _plot_names(result)
        for fn_name in (
            "induction_multiperiod_map_smart",
            "induction_section_smart",
            "response_tipper_smart",
            "tipper_components_smart",
            "tipper_hodograms_smart",
        ):
            assert fn_name not in names, f"{fn_name} fired despite no tipper"

    def test_strike_qc_still_fires_without_tipper(self, result):
        # plot_strike_analysis is internally tipper-aware (2-panel instead
        # of 3) but must still run.
        names = _plot_names(result)
        assert "plot_strike_analysis" in names
        assert "plot_strike_rose" in names


# ─────────────────────────────────────────────────────────────────────────────
# csamt_qc / csumt_qc against real CSAMT data
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _DATA_CSAMT.exists(), reason=f"CSAMT data not found at {_DATA_CSAMT}"
)
class TestCsamtQcRealData:
    @pytest.fixture(scope="class")
    def result(self):
        return _run("csamt_qc", _DATA_CSAMT, _OUTROOT / "csamt_qc")

    def test_run_ok_even_without_resolvable_offset(self, result):
        # data/CSAMT's real headers expose no source_offset -- the preset's
        # source_offset=None default must warn-and-passthrough, never fail
        # the run.
        assert result.ok, "\n".join(
            f"{sr.step_code}: {sr.error}" for sr in result.step_results if not sr.ok
        )

    def test_field_zone_pseudosection_present(self, result):
        assert "plot_field_zones" in _plot_names(result)

    def test_no_depth_section_in_csamt_qc(self, result):
        # Bostick depth section is CSUMT-specific, not part of csamt_qc.
        assert "plot_depth_section" not in _plot_names(result)


@pytest.mark.skipif(
    not _DATA_CSAMT.exists(), reason=f"CSAMT data not found at {_DATA_CSAMT}"
)
class TestCsumtQcRealData:
    @pytest.fixture(scope="class")
    def result(self):
        return _run("csumt_qc", _DATA_CSAMT, _OUTROOT / "csumt_qc")

    def test_run_ok(self, result):
        assert result.ok, "\n".join(
            f"{sr.step_code}: {sr.error}" for sr in result.step_results if not sr.ok
        )

    def test_depth_section_present(self, result):
        assert "plot_depth_section" in _plot_names(result)


@pytest.mark.skipif(
    not _DATA_CSAMT.exists(), reason=f"CSAMT data not found at {_DATA_CSAMT}"
)
class TestNearFieldCorrectionRealData:
    """Prove the correction path itself works (not just the graceful skip)."""

    def test_explicit_offset_actually_changes_z(self):
        sites = _load(_DATA_CSAMT)
        corrected = correct_near_field(sites, source_offset=300.0, verbose=0)

        raw_z = next(iter(_get_z_block(ed)[1] for ed in _iter_items(sites)))
        cor_z = next(iter(_get_z_block(ed)[1] for ed in _iter_items(corrected)))
        assert not np.allclose(raw_z, cor_z, equal_nan=True)

    def test_none_offset_warns_and_passes_through(self):
        sites = _load(_DATA_CSAMT)
        with pytest.warns(UserWarning, match="no source offset"):
            corrected = correct_near_field(sites, source_offset=None, verbose=1)

        raw_z = next(iter(_get_z_block(ed)[1] for ed in _iter_items(sites)))
        cor_z = next(iter(_get_z_block(ed)[1] for ed in _iter_items(corrected)))
        assert np.allclose(raw_z, cor_z, equal_nan=True)

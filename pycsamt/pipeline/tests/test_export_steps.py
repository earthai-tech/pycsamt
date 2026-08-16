"""Real end-to-end tests for the built-in EXPORT001/002/003 pipeline steps.

Public entry point: :mod:`pycsamt.pipeline` (``Step("EXPORT001")`` etc.)
Implementation:     wrapper functions in :mod:`pycsamt.pipeline._registry`

These steps wrap the already-validated ``pycsamt.models.{modem,occam2d,mare2dem}``
input writers (see the doctested tutorials at
``docs/source/tutorials/prepare_{modem,occam2d,mare2dem}_inversion.rst`` for
the same writers exercised against the same real WILLY data with asserted
real numbers). This file proves the *pipeline step* wrapping — that ``Step``
plumbs ``sites``/params through correctly and that ``sites`` passes through
unchanged — not the underlying writers' own correctness, which the tutorials
already cover.

Skip condition: skipped automatically when ``data/AMT/WILLY_DATA`` is absent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.emtools._core import ensure_sites
from pycsamt.pipeline import Step, lookup_step

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_DATA_DIR = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"

pytestmark = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"WILLY EDI data not found at {_DATA_DIR}",
)

_N_STATIONS = 5


@pytest.fixture(scope="module")
def willy_sites():
    paths = sorted(_DATA_DIR.glob("*.edi"))[:_N_STATIONS]
    assert paths, f"No EDI files found in {_DATA_DIR}"
    return ensure_sites([str(p) for p in paths])


class TestExportStepRegistration:
    def test_export001_is_builtin(self):
        spec = lookup_step("EXPORT001")
        assert spec.name == "export_modem"
        assert spec.category == "export"
        assert spec.origin == "builtin"
        assert spec.returns_sites is True

    def test_export002_is_builtin(self):
        spec = lookup_step("EXPORT002")
        assert spec.name == "export_occam2d"
        assert spec.category == "export"
        assert spec.origin == "builtin"

    def test_export003_is_builtin(self):
        spec = lookup_step("EXPORT003")
        assert spec.name == "export_mare2dem"
        assert spec.category == "export"
        assert spec.origin == "builtin"


class TestExportModem:
    def test_writes_real_modem_files(self, willy_sites, tmp_path):
        step = Step("EXPORT001", workdir=str(tmp_path / "modem"))
        result = step.transform(willy_sites)

        assert result is willy_sites  # sites pass through unchanged
        # Real filenames written by InputBuilder.build() with default config
        # (ModEmConfig() defaults to a 3-D run, hence covariance.cov too).
        written = {p.name for p in (tmp_path / "modem").iterdir()}
        assert written == {"data.dat", "m0.ws", "covariance.cov", "control.inv"}

    def test_creates_workdir_if_missing(self, willy_sites, tmp_path):
        workdir = tmp_path / "nested" / "modem_out"
        assert not workdir.exists()
        Step("EXPORT001", workdir=str(workdir)).transform(willy_sites)
        assert workdir.exists()


class TestExportOccam2D:
    def test_writes_real_occam2d_files(self, willy_sites, tmp_path):
        step = Step("EXPORT002", workdir=str(tmp_path / "occam2d"))
        result = step.transform(willy_sites)

        assert result is willy_sites
        written = {p.name for p in (tmp_path / "occam2d").iterdir()}
        assert written == {
            "OccamDataFile.dat",
            "Occam2DMesh",
            "Occam2DModel",
            "Startup",
        }


class TestExportMare2dem:
    def test_writes_real_emdata_file(self, willy_sites, tmp_path):
        step = Step("EXPORT003", workdir=str(tmp_path / "mare2dem"))
        result = step.transform(willy_sites)

        assert result is willy_sites
        assert (tmp_path / "mare2dem" / "line.emdata").exists()

    def test_custom_data_filename(self, willy_sites, tmp_path):
        step = Step(
            "EXPORT003",
            workdir=str(tmp_path / "mare2dem"),
            data_filename="custom.emdata",
        )
        step.transform(willy_sites)
        assert (tmp_path / "mare2dem" / "custom.emdata").exists()

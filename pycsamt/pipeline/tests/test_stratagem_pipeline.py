"""Tests for pycsamt.pipeline.stratagem.

Covers both API levels described in the module docstring:

* Level 1 -- :class:`StratagemPipeline` (a :class:`~._pipeline.Pipeline`
  subclass) and its helper functions ``_coerce_to_sites``,
  ``_inject_coordinates``, ``_apply_hardware_mask``, ``_rename_processed``.
* Level 2 -- :class:`StratagemPreset` / :data:`STRATAGEM_PRESETS` and
  :func:`run_stratagem_preset`.

Synthetic EDI / coordinate-CSV / raw-hardware fixtures follow the patterns
already established in ``pycsamt/stratagem/tests`` (test_gis_correct.py,
test_io.py, test_rename_survey.py) rather than inventing new ones.
"""

from __future__ import annotations

import os
import textwrap
import warnings
from pathlib import Path
from unittest.mock import MagicMock

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

from pycsamt.pipeline._pipeline import Pipeline, PipelineResult
from pycsamt.pipeline._steps import Step
from pycsamt.pipeline.stratagem import (
    STRATAGEM_PRESETS,
    StratagemPipeline,
    StratagemPreset,
    _apply_hardware_mask,
    _coerce_to_sites,
    _inject_coordinates,
    _rename_processed,
    get_stratagem_preset,
    list_stratagem_presets,
    run_stratagem_preset,
    stratagem_preset_catalogue,
)
from pycsamt.seg.edi import EDIFile

# ---------------------------------------------------------------------------
# Shared EDI / CSV factories (mirrors pycsamt/stratagem/tests conventions)
# ---------------------------------------------------------------------------

_EDI_TMPL = textwrap.dedent(
    """\
>HEAD
  DATAID="S{sid:02d}"
  ACQBY=AFEDIG
  LAT={lat:.6f}
  LONG={lon:.6f}
  ELEV=200.0

>INFO
  MAXINFO=999

>=DEFINEMEAS
  MAXCHAN=5
  UNITS=M
  REFTYPE=CART
  REFLAT={lat:.6f}
  REFLONG={lon:.6f}
  REFELEV=200.0

>=MTSECT
  SECTID=S{sid:02d}
  NFREQ=4
  HX=1.001
  HY=2.001

>!****FREQUENCIES****!
>FREQ  //4
   1.000000E+02   1.000000E+01   1.000000E+00   1.000000E-01

>ZROT  //4
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

>ZXX  //4
   0.0  0.0  0.0  0.0

>ZXY  //4
   1.0  2.0  3.0  4.0

>ZYX  //4
  -1.0 -2.0 -3.0 -4.0

>ZYY  //4
   0.0  0.0  0.0  0.0

>ZXX.VAR  //4
   0.01 0.01 0.01 0.01

>ZXY.VAR  //4
   0.01 0.01 0.01 0.01

>ZYX.VAR  //4
   0.01 0.01 0.01 0.01

>ZYY.VAR  //4
   0.01 0.01 0.01 0.01

>END
"""
)

_RAW_19COL = textwrap.dedent(
    """\
 1.130e+001 2.930e+000 2.400e+001  3.728e+001  2.152e+002 -5.546e-001 -1.695e+002  1.336e+002  3.336e+004 -2.382e+002 -3.418e+003  3.016e+001  1.338e+002  2.845e+001 -1.052e+002 -2.512e+002 -1.384e+004 -2.318e+002  1.866e+004
 1.250e+001 2.930e+000 2.100e+001  2.248e+001  3.369e+002 -5.041e-002 -8.960e+001  1.700e+002  2.217e+004 -3.005e+002  5.244e+002  2.076e+001  1.842e+002  2.233e+001 -9.249e+001 -3.353e+002 -1.111e+004 -3.027e+002  1.656e+004
"""
)


def _make_edi_dir(tmp_path: Path, n: int = 3, subdir: str = "edis") -> Path:
    d = tmp_path / subdir
    d.mkdir(exist_ok=True)
    for i in range(n):
        lat = 25.77
        lon = 109.60 + i * 0.02
        (d / f"Z2HX{i + 1:03d}.edi").write_text(
            _EDI_TMPL.format(sid=i, lat=lat, lon=lon), encoding="utf-8"
        )
    return d


def _make_coord_csv(tmp_path: Path, n: int = 3, name: str = "coords.csv") -> Path:
    rows = []
    for i in range(n):
        easting = 362_589.0 + i * 20.0
        northing = 2_850_835.0 - i * 20.0
        rows.append(
            {
                "stations": f"K0+{i * 20:06.4f}",
                "longitude": northing,
                "latitude": easting,
                "elev": 261.94 + i,
                "step": i * 20.0,
            }
        )
    csv_path = tmp_path / name
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return csv_path


def _make_raw_dir(tmp_path: Path, n: int = 3, subdir: str = "raw") -> Path:
    d = tmp_path / subdir
    d.mkdir(exist_ok=True)
    for i in range(n):
        (d / f"X2HX.{i + 1:03d}").write_text(_RAW_19COL, encoding="utf-8")
    return d


def _load_edi_objects(tmp_path: Path, n: int = 3, subdir: str = "edis"):
    d = _make_edi_dir(tmp_path, n=n, subdir=subdir)
    return [EDIFile(p) for p in sorted(d.glob("*.edi"))]


# ---------------------------------------------------------------------------
# StratagemPreset / STRATAGEM_PRESETS
# ---------------------------------------------------------------------------


class TestStratagemPresetStructure:
    def test_repr_shows_name_and_steps(self):
        preset = STRATAGEM_PRESETS["basic"]
        r = repr(preset)
        assert "StratagemPreset('basic'" in r
        assert "remove_static_shift" in r
        assert "export" in r
        assert "rename" in r

    def test_three_presets_registered(self):
        assert set(STRATAGEM_PRESETS) == {
            "basic",
            "full_processing",
            "publication_ready",
        }

    @pytest.mark.parametrize(
        "name", ["basic", "full_processing", "publication_ready"]
    )
    def test_each_preset_is_stratagem_preset(self, name):
        preset = STRATAGEM_PRESETS[name]
        assert isinstance(preset, StratagemPreset)
        assert preset.name == name
        assert isinstance(preset.description, str) and preset.description

    @pytest.mark.parametrize(
        "name", ["basic", "full_processing", "publication_ready"]
    )
    def test_each_preset_has_survey_defaults(self, name):
        preset = STRATAGEM_PRESETS[name]
        assert preset.survey_defaults == {"epsg": 32649, "utm_zone": "49N"}

    @pytest.mark.parametrize(
        "name", ["basic", "full_processing", "publication_ready"]
    )
    def test_each_preset_steps_end_with_export_rename(self, name):
        steps = STRATAGEM_PRESETS[name].steps
        assert steps[-2][0] == "export"
        assert steps[-1][0] == "rename"
        for method_name, kw in steps:
            assert isinstance(method_name, str)
            assert isinstance(kw, dict)

    def test_full_processing_includes_run_qc(self):
        codes = [n for n, _ in STRATAGEM_PRESETS["full_processing"].steps]
        assert "run_qc" in codes

    def test_publication_ready_has_tighter_band_than_basic(self):
        def _fmin(name):
            for n, kw in STRATAGEM_PRESETS[name].steps:
                if n == "drop_frequencies":
                    return kw["fmin"]
            raise AssertionError("drop_frequencies step missing")

        assert _fmin("publication_ready") > _fmin("basic")


class TestGetListCatalogue:
    def test_get_valid_preset(self):
        preset = get_stratagem_preset("basic")
        assert preset is STRATAGEM_PRESETS["basic"]

    def test_get_invalid_preset_raises_key_error(self):
        with pytest.raises(KeyError, match="Unknown Stratagem preset"):
            get_stratagem_preset("not_a_real_preset")

    def test_list_stratagem_presets_returns_all_three(self):
        presets = list_stratagem_presets()
        assert len(presets) == 3
        assert all(isinstance(p, StratagemPreset) for p in presets)

    def test_catalogue_is_string_with_all_names_and_steps(self):
        cat = stratagem_preset_catalogue()
        assert isinstance(cat, str)
        for name in ("basic", "full_processing", "publication_ready"):
            assert name in cat
        assert "remove_static_shift" in cat
        assert "Stratagem pipeline presets" in cat


# ---------------------------------------------------------------------------
# _coerce_to_sites
# ---------------------------------------------------------------------------


class TestCoerceToSites:
    def test_sites_object_passthrough(self, tmp_path):
        from pycsamt.site.base import Sites

        edis = _load_edi_objects(tmp_path, n=3)
        sites_in = Sites(edis)
        out = _coerce_to_sites(sites_in)
        assert isinstance(out, Sites)
        assert len(out) == 3

    def test_dir_path(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        out = _coerce_to_sites(str(d))
        assert len(out) == 3

    def test_dir_path_as_pathlib(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        out = _coerce_to_sites(d)
        assert len(out) == 3

    def test_single_edi_file_path(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=1)
        edi_path = next(d.glob("*.edi"))
        out = _coerce_to_sites(str(edi_path))
        assert len(out) == 1

    def test_single_edifile_object(self, tmp_path):
        edis = _load_edi_objects(tmp_path, n=1)
        out = _coerce_to_sites(edis[0])
        assert len(out) == 1

    def test_list_of_edifile_objects(self, tmp_path):
        edis = _load_edi_objects(tmp_path, n=3)
        out = _coerce_to_sites(edis)
        assert len(out) == 3


# ---------------------------------------------------------------------------
# _inject_coordinates
# ---------------------------------------------------------------------------


class TestInjectCoordinates:
    # The real success path (CoordinateInjector.fit -> pyproj) is exercised
    # via the mocked test below rather than a real pyproj call: on this
    # environment's pyproj/PROJ build, the legacy `+init=EPSG:` codepath in
    # gis.utils.project_point_utm2ll (no GDAL available) raises a spurious
    # ProjVersion error specifically when a coverage tracer is attached --
    # reproducible even against the pre-existing pycsamt/stratagem/tests
    # suite, so it is a pyproj/coverage.py interaction bug, not a bug in
    # this file's code.
    def test_row_count_mismatch_warns_and_returns_original(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        csv = _make_coord_csv(tmp_path, n=5)  # mismatched row count
        sites_in = _coerce_to_sites(str(d))

        with pytest.warns(UserWarning, match="coordinate injection skipped"):
            out = _inject_coordinates(
                sites_in, csv, epsg=32649, utm_zone="49N", order="auto"
            )

        assert out is sites_in

    def test_mocked_injector_returns_ensure_sites_of_edi_objects(
        self, tmp_path, monkeypatch
    ):
        d = _make_edi_dir(tmp_path, n=2)
        sites_in = _coerce_to_sites(str(d))
        edi_objs = _load_edi_objects(tmp_path, n=2, subdir="edis")

        fake_injector = MagicMock()
        fake_injector.fit.return_value = fake_injector
        fake_injector.edi_objects_ = edi_objs

        monkeypatch.setattr(
            "pycsamt.stratagem.gis_correct.CoordinateInjector",
            MagicMock(return_value=fake_injector),
        )

        out = _inject_coordinates(
            sites_in, "unused.csv", epsg=32649, utm_zone="49N", order="auto"
        )
        assert len(out) == 2
        fake_injector.fit.assert_called_once()


# ---------------------------------------------------------------------------
# _apply_hardware_mask
# ---------------------------------------------------------------------------


class TestApplyHardwareMask:
    def test_real_small_dataset_runs_and_returns_same_count(self, tmp_path):
        edi_dir = _make_edi_dir(tmp_path, n=3)
        raw_dir = _make_raw_dir(tmp_path, n=3)
        sites_in = _coerce_to_sites(str(edi_dir))

        out = _apply_hardware_mask(sites_in, raw_dir)
        assert len(out) == 3

    def test_mocked_reader_and_filter(self, tmp_path, monkeypatch):
        d = _make_edi_dir(tmp_path, n=2)
        sites_in = _coerce_to_sites(str(d))
        edi_objs = _load_edi_objects(tmp_path, n=2, subdir="edis")

        fake_reader = MagicMock()
        fake_reader.fit.return_value = fake_reader
        monkeypatch.setattr(
            "pycsamt.stratagem.io.StratagemRawReader",
            MagicMock(return_value=fake_reader),
        )

        fake_filter = MagicMock()
        fake_filter.fit.return_value = fake_filter
        fake_filter.edi_objects_ = edi_objs
        monkeypatch.setattr(
            "pycsamt.stratagem.qc.FrequencyFilter",
            MagicMock(return_value=fake_filter),
        )

        out = _apply_hardware_mask(sites_in, "unused_raw_dir")
        assert len(out) == 2
        fake_filter.fit.assert_called_once()
        _, call_kwargs = fake_filter.fit.call_args
        assert call_kwargs["raw_reader"] is fake_reader


# ---------------------------------------------------------------------------
# _rename_processed
# ---------------------------------------------------------------------------


class TestRenameProcessed:
    def test_nonexistent_processed_dir_returns_empty_list(self, tmp_path):
        out = _rename_processed(
            tmp_path / "does_not_exist", tmp_path / "renamed", "T."
        )
        assert out == []

    def test_real_rename_creates_files(self, tmp_path):
        processed = _make_edi_dir(tmp_path, n=3, subdir="processed")
        rename_dir = tmp_path / "renamed"

        out = _rename_processed(processed, rename_dir, "T.", overwrite=False)

        assert len(out) == 3
        assert (rename_dir / "T.000.edi").exists()
        assert (rename_dir / "T.002.edi").exists()

    def test_overwrite_flag_forwarded(self, tmp_path):
        processed = _make_edi_dir(tmp_path, n=2, subdir="processed")
        rename_dir = tmp_path / "renamed"

        _rename_processed(processed, rename_dir, "S.", overwrite=False)
        second = _rename_processed(processed, rename_dir, "S.", overwrite=True)
        assert len(second) == 2


# ---------------------------------------------------------------------------
# StratagemPipeline.__init__
# ---------------------------------------------------------------------------


class TestStratagemPipelineInit:
    def test_attribute_wiring(self):
        pipe = StratagemPipeline(
            [("qc", Step("QC001"))],
            coord_file="coords.csv",
            raw_dir="raw/",
            epsg=1234,
            utm_zone="50N",
            order="reversed",
            rename_basename="T.",
            rename_dir="renamed/",
            name="my_pipe",
        )
        assert pipe.coord_file == "coords.csv"
        assert pipe.raw_dir == "raw/"
        assert pipe.epsg == 1234
        assert pipe.utm_zone == "50N"
        assert pipe.order == "reversed"
        assert pipe.rename_basename == "T."
        assert pipe.rename_dir == "renamed/"
        assert pipe.name == "my_pipe"
        assert isinstance(pipe, Pipeline)

    def test_defaults(self):
        pipe = StratagemPipeline([("qc", Step("QC001"))])
        assert pipe.coord_file is None
        assert pipe.raw_dir is None
        assert pipe.epsg == 32649
        assert pipe.utm_zone == "49N"
        assert pipe.order == "auto"
        assert pipe.rename_basename is None
        assert pipe.rename_dir is None
        assert pipe.name == "stratagem_pipeline"


# ---------------------------------------------------------------------------
# StratagemPipeline.run -- orchestration, isolated with mocks
# ---------------------------------------------------------------------------


class TestStratagemPipelineRunOrchestration:
    def _pipe(self, **kwargs):
        return StratagemPipeline([("qc", Step("QC001"))], **kwargs)

    def test_no_coord_no_raw_skips_pre_processing(self, tmp_path, monkeypatch):
        pipe = self._pipe()
        sentinel_sites = object()

        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites",
            MagicMock(return_value=sentinel_sites),
        )
        inject_mock = MagicMock()
        mask_mock = MagicMock()
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._inject_coordinates", inject_mock
        )
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._apply_hardware_mask", mask_mock
        )

        fake_result = PipelineResult(
            sites_in=sentinel_sites,
            sites_out=sentinel_sites,
            step_results=[],
            outdir=None,
            elapsed_sec=0.01,
        )
        super_run_mock = MagicMock(return_value=fake_result)
        monkeypatch.setattr(Pipeline, "run", super_run_mock)

        result = pipe.run("whatever", outdir=None)

        inject_mock.assert_not_called()
        mask_mock.assert_not_called()
        super_run_mock.assert_called_once()
        _, call_kwargs = super_run_mock.call_args
        assert call_kwargs["outdir"] is None
        assert result is fake_result

    def test_coord_and_raw_both_applied_in_order(self, tmp_path, monkeypatch):
        pipe = self._pipe(coord_file="c.csv", raw_dir="raw/", epsg=999, utm_zone="1N", order="reversed")

        s0, s1, s2 = object(), object(), object()
        coerce_mock = MagicMock(return_value=s0)
        inject_mock = MagicMock(return_value=s1)
        mask_mock = MagicMock(return_value=s2)
        monkeypatch.setattr("pycsamt.pipeline.stratagem._coerce_to_sites", coerce_mock)
        monkeypatch.setattr("pycsamt.pipeline.stratagem._inject_coordinates", inject_mock)
        monkeypatch.setattr("pycsamt.pipeline.stratagem._apply_hardware_mask", mask_mock)

        fake_result = PipelineResult(
            sites_in=s2, sites_out=s2, step_results=[], outdir=None, elapsed_sec=0.01
        )
        super_run_mock = MagicMock(return_value=fake_result)
        monkeypatch.setattr(Pipeline, "run", super_run_mock)

        pipe.run("input_sites", outdir=None)

        coerce_mock.assert_called_once_with("input_sites")
        inject_mock.assert_called_once_with(s0, "c.csv", 999, "1N", "reversed")
        mask_mock.assert_called_once_with(s1, "raw/")
        _, call_kwargs = super_run_mock.call_args
        assert call_kwargs is not None
        args, kwargs = super_run_mock.call_args
        assert args[0] is s2

    def test_rename_fires_when_basename_set_and_outdir_present(
        self, tmp_path, monkeypatch
    ):
        pipe = self._pipe(rename_basename="T.")
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites", MagicMock(return_value="S")
        )

        real_outdir = tmp_path / "out"
        real_outdir.mkdir()
        fake_result = PipelineResult(
            sites_in="S",
            sites_out="S",
            step_results=[],
            outdir=real_outdir,
            elapsed_sec=0.01,
            processed_paths=[],
        )
        monkeypatch.setattr(Pipeline, "run", MagicMock(return_value=fake_result))

        rename_mock = MagicMock(return_value=[Path("a.edi"), Path("b.edi")])
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._rename_processed", rename_mock
        )

        result = pipe.run("sites", outdir=real_outdir)

        rename_mock.assert_called_once()
        args, kwargs = rename_mock.call_args
        processed_dir_arg, dst_dir_arg, basename_arg = args[0], args[1], args[2]
        assert processed_dir_arg == real_outdir / "processed"
        assert dst_dir_arg == real_outdir / "renamed"
        assert basename_arg == "T."
        assert kwargs["overwrite"] is False
        assert result.processed_paths == [Path("a.edi"), Path("b.edi")]

    def test_rename_basename_param_overrides_self(self, tmp_path, monkeypatch):
        pipe = self._pipe(rename_basename="SELF.")
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites", MagicMock(return_value="S")
        )
        real_outdir = tmp_path / "out"
        real_outdir.mkdir()
        fake_result = PipelineResult(
            sites_in="S", sites_out="S", step_results=[], outdir=real_outdir,
            elapsed_sec=0.01,
        )
        monkeypatch.setattr(Pipeline, "run", MagicMock(return_value=fake_result))
        rename_mock = MagicMock(return_value=[])
        monkeypatch.setattr("pycsamt.pipeline.stratagem._rename_processed", rename_mock)

        pipe.run("sites", outdir=real_outdir, rename_basename="OVERRIDE.")

        args, _ = rename_mock.call_args
        assert args[2] == "OVERRIDE."

    def test_rename_dir_param_overrides_self(self, tmp_path, monkeypatch):
        pipe = self._pipe(rename_basename="T.", rename_dir=tmp_path / "self_renamed")
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites", MagicMock(return_value="S")
        )
        real_outdir = tmp_path / "out"
        real_outdir.mkdir()
        fake_result = PipelineResult(
            sites_in="S", sites_out="S", step_results=[], outdir=real_outdir,
            elapsed_sec=0.01,
        )
        monkeypatch.setattr(Pipeline, "run", MagicMock(return_value=fake_result))
        rename_mock = MagicMock(return_value=[])
        monkeypatch.setattr("pycsamt.pipeline.stratagem._rename_processed", rename_mock)

        override_dir = tmp_path / "param_renamed"
        pipe.run("sites", outdir=real_outdir, rename_dir=override_dir)

        args, _ = rename_mock.call_args
        assert args[1] == override_dir.expanduser().resolve()

    def test_overwrite_flag_forwarded_to_rename(self, tmp_path, monkeypatch):
        pipe = self._pipe(rename_basename="T.")
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites", MagicMock(return_value="S")
        )
        real_outdir = tmp_path / "out"
        real_outdir.mkdir()
        fake_result = PipelineResult(
            sites_in="S", sites_out="S", step_results=[], outdir=real_outdir,
            elapsed_sec=0.01,
        )
        monkeypatch.setattr(Pipeline, "run", MagicMock(return_value=fake_result))
        rename_mock = MagicMock(return_value=[])
        monkeypatch.setattr("pycsamt.pipeline.stratagem._rename_processed", rename_mock)

        pipe.run("sites", outdir=real_outdir, overwrite=True)

        _, kwargs = rename_mock.call_args
        assert kwargs["overwrite"] is True

    def test_no_rename_when_basename_unset(self, tmp_path, monkeypatch):
        pipe = self._pipe()
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites", MagicMock(return_value="S")
        )
        real_outdir = tmp_path / "out"
        real_outdir.mkdir()
        fake_result = PipelineResult(
            sites_in="S", sites_out="S", step_results=[], outdir=real_outdir,
            elapsed_sec=0.01,
        )
        monkeypatch.setattr(Pipeline, "run", MagicMock(return_value=fake_result))
        rename_mock = MagicMock()
        monkeypatch.setattr("pycsamt.pipeline.stratagem._rename_processed", rename_mock)

        pipe.run("sites", outdir=real_outdir)

        rename_mock.assert_not_called()

    def test_no_rename_when_outdir_none_even_with_basename(self, monkeypatch):
        pipe = self._pipe(rename_basename="T.")
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem._coerce_to_sites", MagicMock(return_value="S")
        )
        fake_result = PipelineResult(
            sites_in="S", sites_out="S", step_results=[], outdir=None, elapsed_sec=0.01
        )
        monkeypatch.setattr(Pipeline, "run", MagicMock(return_value=fake_result))
        rename_mock = MagicMock()
        monkeypatch.setattr("pycsamt.pipeline.stratagem._rename_processed", rename_mock)

        pipe.run("sites", outdir=None)

        rename_mock.assert_not_called()


# ---------------------------------------------------------------------------
# StratagemPipeline.run -- real, lightweight end-to-end
# ---------------------------------------------------------------------------


class TestStratagemPipelineRunRealLight:
    def test_real_run_no_preprocessing_no_outdir(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        pipe = StratagemPipeline([("qc", Step("QC001"))])

        result = pipe.run(str(d), outdir=None, save_edis=False, save_report=False)

        assert isinstance(result, PipelineResult)
        assert result.ok
        assert len(result.sites_out) == 3

    def test_real_run_with_outdir_and_rename(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        outdir = tmp_path / "pipe_out"
        pipe = StratagemPipeline(
            [("qc", Step("QC001"))], rename_basename="T."
        )

        result = pipe.run(str(d), outdir=outdir, save_plots=False, save_report=False)

        assert result.ok
        assert (outdir / "renamed" / "T.000.edi").exists()
        assert len(result.processed_paths) >= 3

    # A real (non-mismatched) coord_file run is intentionally not covered
    # here -- see the note on TestInjectCoordinates for why: it would call
    # the same pyproj codepath that misbehaves under a coverage tracer on
    # this environment. The mismatch path below stays real since it raises
    # ValidationError before ever reaching pyproj, and the orchestration
    # tests above already prove StratagemPipeline.run wires coord_file
    # through to _inject_coordinates correctly.
    def test_real_run_with_coord_mismatch_warns_but_completes(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        csv = _make_coord_csv(tmp_path, n=7)
        pipe = StratagemPipeline(
            [("qc", Step("QC001"))], coord_file=csv,
        )

        with pytest.warns(UserWarning, match="coordinate injection skipped"):
            result = pipe.run(str(d), outdir=None, save_edis=False, save_report=False)

        assert result.ok
        assert len(result.sites_out) == 3

    def test_real_run_with_hardware_mask(self, tmp_path):
        d = _make_edi_dir(tmp_path, n=3)
        raw_dir = _make_raw_dir(tmp_path, n=3)
        pipe = StratagemPipeline([("qc", Step("QC001"))], raw_dir=raw_dir)

        result = pipe.run(str(d), outdir=None, save_edis=False, save_report=False)

        assert result.ok
        assert len(result.sites_out) == 3


# ---------------------------------------------------------------------------
# StratagemPipeline.from_preset
# ---------------------------------------------------------------------------


class TestFromPreset:
    def test_from_preset_real_stratagem_mt(self):
        pipe = StratagemPipeline.from_preset(
            "stratagem_mt", coord_file="2.csv", epsg=32649, rename_basename="T."
        )
        assert isinstance(pipe, StratagemPipeline)
        assert pipe.coord_file == "2.csv"
        assert pipe.epsg == 32649
        assert pipe.rename_basename == "T."
        assert pipe.name == "stratagem_mt"
        assert len(pipe) > 0

    def test_pipeline_name_override(self):
        pipe = StratagemPipeline.from_preset(
            "stratagem_mt", pipeline_name="custom_name"
        )
        assert pipe.name == "custom_name"

    def test_from_preset_forwards_all_kwargs(self, monkeypatch):
        from pycsamt.pipeline._presets import Preset

        fake_preset = Preset(
            name="fake_preset",
            description="fake",
            steps=[("qc", Step("QC001"))],
        )
        monkeypatch.setattr(
            "pycsamt.pipeline._presets.get_preset",
            MagicMock(return_value=fake_preset),
        )

        pipe = StratagemPipeline.from_preset(
            "fake_preset",
            coord_file="c.csv",
            raw_dir="r/",
            epsg=1,
            utm_zone="2N",
            order="mapping",
            rename_basename="B.",
            rename_dir="d/",
        )
        assert pipe.coord_file == "c.csv"
        assert pipe.raw_dir == "r/"
        assert pipe.epsg == 1
        assert pipe.utm_zone == "2N"
        assert pipe.order == "mapping"
        assert pipe.rename_basename == "B."
        assert pipe.rename_dir == "d/"
        assert pipe.name == "fake_preset"
        assert len(pipe) == 1

    def test_from_preset_invalid_name_raises(self):
        with pytest.raises(KeyError):
            StratagemPipeline.from_preset("not_a_real_preset")


# ---------------------------------------------------------------------------
# run_stratagem_preset
# ---------------------------------------------------------------------------


def _spec_mock_survey():
    from pycsamt.stratagem.survey import StratagemSurvey

    inst = MagicMock(spec=StratagemSurvey)
    inst.fit.return_value = inst
    return inst


class TestRunStratagemPreset:
    def test_basic_preset_calls_expected_methods_and_defaults(
        self, tmp_path, monkeypatch
    ):
        inst = _spec_mock_survey()
        cls_mock = MagicMock(return_value=inst)
        monkeypatch.setattr("pycsamt.stratagem.survey.StratagemSurvey", cls_mock)

        outdir = tmp_path / "out"
        sv = run_stratagem_preset(
            "basic",
            edi_dir="edi/",
            coord_file="c.csv",
            outdir=outdir,
        )

        assert sv is inst
        cls_mock.assert_called_once_with(
            edi_dir="edi/",
            coord_file="c.csv",
            raw_dir=None,
            verbose=0,
            epsg=32649,
            utm_zone="49N",
        )
        inst.fit.assert_called_once()
        inst.remove_static_shift.assert_called_once_with(half_window=3, weights="tri")
        inst.drop_frequencies.assert_called_once_with(
            fmin=10.0, snr_thresh=2.5, min_frac=0.4, use_hardware_mask=True
        )
        inst.remove_noises.assert_called_once_with(
            mains_hz=50.0, n_harm=30, hampel_win=3, smooth=False
        )
        inst.export.assert_called_once_with(outdir / "corrected", overwrite=False)
        inst.rename.assert_called_once_with(
            basename="S", dst_path=outdir / "renamed", overwrite=False
        )

    def test_explicit_epsg_utm_and_raw_dir_forwarded(self, tmp_path, monkeypatch):
        inst = _spec_mock_survey()
        cls_mock = MagicMock(return_value=inst)
        monkeypatch.setattr("pycsamt.stratagem.survey.StratagemSurvey", cls_mock)

        run_stratagem_preset(
            "basic",
            edi_dir="edi/",
            coord_file="c.csv",
            outdir=tmp_path / "out",
            raw_dir="raw/",
            epsg=1111,
            utm_zone="9N",
        )
        cls_mock.assert_called_once_with(
            edi_dir="edi/",
            coord_file="c.csv",
            raw_dir="raw/",
            verbose=0,
            epsg=1111,
            utm_zone="9N",
        )

    def test_rename_basename_and_dir_and_overwrite_forwarded(
        self, tmp_path, monkeypatch
    ):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )
        custom_rename_dir = tmp_path / "custom_renamed"
        run_stratagem_preset(
            "basic",
            edi_dir="edi/",
            coord_file="c.csv",
            outdir=tmp_path / "out",
            rename_basename="T.",
            rename_dir=custom_rename_dir,
            overwrite=True,
        )
        inst.export.assert_called_once_with(
            (tmp_path / "out" / "corrected"), overwrite=True
        )
        inst.rename.assert_called_once_with(
            basename="T.",
            dst_path=custom_rename_dir.expanduser().resolve(),
            overwrite=True,
        )

    def test_rename_dir_defaults_to_outdir_renamed(self, tmp_path, monkeypatch):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )
        outdir = tmp_path / "out"
        run_stratagem_preset(
            "basic", edi_dir="edi/", coord_file="c.csv", outdir=outdir
        )
        _, kwargs = inst.rename.call_args
        assert kwargs["dst_path"] == outdir.expanduser().resolve() / "renamed"

    def test_step_overrides_merge_with_preset_defaults(self, tmp_path, monkeypatch):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )
        run_stratagem_preset(
            "basic",
            edi_dir="edi/",
            coord_file="c.csv",
            outdir=tmp_path / "out",
            step_overrides={"remove_static_shift": {"half_window": 9}},
        )
        inst.remove_static_shift.assert_called_once_with(half_window=9, weights="tri")

    def test_full_processing_preset_calls_run_qc(self, tmp_path, monkeypatch):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )
        run_stratagem_preset(
            "full_processing",
            edi_dir="edi/",
            coord_file="c.csv",
            outdir=tmp_path / "out",
        )
        inst.run_qc.assert_called_once_with(
            include_skew=True, min_frac_ok=0.6, min_snr_med=2.0, max_skew_med=6.0
        )

    def test_unknown_preset_raises_key_error(self, tmp_path):
        with pytest.raises(KeyError):
            run_stratagem_preset(
                "not_a_preset", edi_dir="edi/", coord_file="c.csv",
                outdir=tmp_path / "out",
            )

    def test_unknown_step_name_skipped_and_verbose_prints(
        self, tmp_path, monkeypatch, capsys
    ):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )

        bogus_preset = StratagemPreset(
            name="bogus",
            description="has an unknown step",
            survey_defaults={"epsg": 32649, "utm_zone": "49N"},
            steps=[("totally_bogus_method", {})],
        )
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem.get_stratagem_preset",
            MagicMock(return_value=bogus_preset),
        )

        run_stratagem_preset(
            "bogus", edi_dir="edi/", coord_file="c.csv", outdir=tmp_path / "out",
            verbose=1,
        )
        captured = capsys.readouterr()
        assert "unknown step" in captured.out
        assert "totally_bogus_method" in captured.out

    def test_unknown_step_name_silent_without_verbose(
        self, tmp_path, monkeypatch, capsys
    ):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )
        bogus_preset = StratagemPreset(
            name="bogus2",
            description="has an unknown step",
            survey_defaults={},
            steps=[("totally_bogus_method", {})],
        )
        monkeypatch.setattr(
            "pycsamt.pipeline.stratagem.get_stratagem_preset",
            MagicMock(return_value=bogus_preset),
        )

        run_stratagem_preset(
            "bogus2", edi_dir="edi/", coord_file="c.csv", outdir=tmp_path / "out",
            verbose=0,
        )
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_returns_the_survey_instance(self, tmp_path, monkeypatch):
        inst = _spec_mock_survey()
        monkeypatch.setattr(
            "pycsamt.stratagem.survey.StratagemSurvey", MagicMock(return_value=inst)
        )
        sv = run_stratagem_preset(
            "basic", edi_dir="edi/", coord_file="c.csv", outdir=tmp_path / "out"
        )
        assert sv is inst

    # A fully-real (no-mock) end-to-end run is intentionally not exercised
    # here: StratagemSurvey.fit() always calls CoordinateInjector, which
    # hits the same pyproj/coverage.py interaction bug noted on
    # TestInjectCoordinates above. The mocked tests in this class already
    # give full statement/branch coverage of run_stratagem_preset's own
    # orchestration logic; the real StratagemSurvey behavior itself is
    # covered by pycsamt/stratagem/tests/test_rename_survey.py.

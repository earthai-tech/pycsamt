"""Tests for stratagem.rename and stratagem.survey."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# shared EDI factory (same as test_qc_process.py)
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
  NFREQ=8
  HX=1.001
  HY=2.001

>!****FREQUENCIES****!
>FREQ  //8
   1.000000E+04   1.000000E+03   1.000000E+02   1.000000E+01
   1.000000E+00   1.000000E-01   1.000000E-02   1.000000E-03

>ZROT  //8
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00
   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00

>ZXX  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZXY  //8
   1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0

>ZYX  //8
  -1.0 -2.0 -3.0 -4.0 -5.0 -6.0 -7.0 -8.0

>ZYY  //8
   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

>ZXX.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>ZXY.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>ZYX.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>ZYY.VAR  //8
   0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01

>END
"""
)

_COORD_CSV_TMPL = "stations,longitude,latitude,elev,step\n"


def _make_edi_dir(tmp_path: Path, n: int = 4) -> Path:
    d = tmp_path / "edis"
    d.mkdir(exist_ok=True)
    for i in range(n):
        lat = 25.77
        lon = 109.60 + i * 0.02
        (d / f"Z2HX{i + 1:03d}.edi").write_text(
            _EDI_TMPL.format(sid=i, lat=lat, lon=lon), encoding="utf-8"
        )
    return d


def _make_coord_csv(tmp_path: Path, n: int = 4) -> Path:
    """Synthetic CSV in the Stratagem survey column format."""
    import pandas as pd

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
    csv = tmp_path / "coords.csv"
    pd.DataFrame(rows).to_csv(csv, index=False)
    return csv


def _load_edis(tmp_path: Path, n: int = 4):
    from pycsamt.stratagem.io import EDIBatch

    _make_edi_dir(tmp_path, n=n)
    return EDIBatch(tmp_path / "edis").fit().edi_objects_


# ---------------------------------------------------------------------------
# EDIRenamer
# ---------------------------------------------------------------------------


class TestEDIRenamer:
    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        _make_edi_dir(tmp_path, n=3)
        rn = EDIRenamer(basename="T.", zero_pad=3)
        assert rn.fit(tmp_path / "edis", tmp_path / "renamed") is rn

    def test_creates_correct_files(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        _make_edi_dir(tmp_path, n=3)
        out = tmp_path / "renamed"
        rn = EDIRenamer(basename="T.", zero_pad=3).fit(tmp_path / "edis", out)
        assert rn.n_renamed_ == 3
        expected = ["T.000.edi", "T.001.edi", "T.002.edi"]
        for name in expected:
            assert (out / name).exists(), f"{name} not found"

    def test_dataid_updated(self, tmp_path):
        from pycsamt.seg.edi import EDIFile
        from pycsamt.stratagem.rename import EDIRenamer

        _make_edi_dir(tmp_path, n=2)
        out = tmp_path / "renamed"
        EDIRenamer(basename="X.", zero_pad=3).fit(tmp_path / "edis", out)
        edi = EDIFile(out / "X.000.edi")
        assert edi.station == "X.000"

    def test_fit_with_edi_objects(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        edis = _load_edis(tmp_path, n=3)
        out = tmp_path / "renamed2"
        rn = EDIRenamer(basename="S", zero_pad=3).fit(edis, out)
        assert rn.n_renamed_ == 3

    def test_dst_paths_returns_list(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        _make_edi_dir(tmp_path, n=2)
        out = tmp_path / "renamed"
        rn = EDIRenamer(basename="A").fit(tmp_path / "edis", out)
        paths = rn.dst_paths()
        assert isinstance(paths, list) and len(paths) == 2

    def test_no_overwrite_skips(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        _make_edi_dir(tmp_path, n=2)
        out = tmp_path / "renamed"
        EDIRenamer(basename="T.").fit(tmp_path / "edis", out)
        rn2 = EDIRenamer(basename="T.", overwrite=False).fit(
            tmp_path / "edis", out
        )
        assert len(rn2.skipped_) == 2

    def test_overwrite_replaces(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        _make_edi_dir(tmp_path, n=2)
        out = tmp_path / "renamed"
        EDIRenamer(basename="T.").fit(tmp_path / "edis", out)
        rn = EDIRenamer(basename="T.", overwrite=True).fit(
            tmp_path / "edis", out
        )
        assert len(rn.skipped_) == 0

    def test_missing_source_raises(self, tmp_path):
        from pycsamt.stratagem.rename import EDIRenamer

        with pytest.raises(Exception):
            EDIRenamer().fit("/does/not/exist", tmp_path / "out")

    def test_dst_paths_before_fit_raises(self):
        from pycsamt.stratagem.rename import EDIRenamer

        with pytest.raises(Exception):
            EDIRenamer().dst_paths()


# ---------------------------------------------------------------------------
# EDIWriter
# ---------------------------------------------------------------------------


class TestEDIWriter:
    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.rename import EDIWriter

        edis = _load_edis(tmp_path, n=3)
        wr = EDIWriter()
        assert wr.fit(edis, tmp_path / "out") is wr

    def test_writes_correct_count(self, tmp_path):
        from pycsamt.stratagem.rename import EDIWriter

        edis = _load_edis(tmp_path, n=3)
        out = tmp_path / "out"
        wr = EDIWriter().fit(edis, out)
        assert wr.n_written_ == 3

    def test_dataid_prefix_applied(self, tmp_path):
        from pycsamt.seg.edi import EDIFile
        from pycsamt.stratagem.rename import EDIWriter

        edis = _load_edis(tmp_path, n=2)
        out = tmp_path / "out"
        EDIWriter(dataid_prefix="X", zero_pad=3).fit(edis, out)
        # files written as X000.edi, X001.edi
        edi = EDIFile(out / "X000.edi")
        assert edi.station == "X000"

    def test_output_files_exist(self, tmp_path):
        from pycsamt.stratagem.rename import EDIWriter

        edis = _load_edis(tmp_path, n=3)
        out = tmp_path / "out"
        wr = EDIWriter().fit(edis, out)
        assert all(p.exists() for p in wr.written_)


# ---------------------------------------------------------------------------
# StratagemSurvey
# ---------------------------------------------------------------------------


class TestStratagemSurvey:
    def _setup(self, tmp_path: Path, n: int = 5):
        _make_edi_dir(tmp_path, n=n)
        csv = _make_coord_csv(tmp_path, n=n)
        return tmp_path / "edis", csv

    def test_fit_returns_self(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649)
        assert sv.fit() is sv

    def test_fit_populates_batch_and_injector(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        assert sv.batch_ is not None
        assert sv.injector_ is not None
        assert sv.edi_objects_ is not None
        assert len(sv.edi_objects_) == 5

    def test_n_stations(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path, n=4)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        assert sv.n_stations_ == 4

    def test_run_qc_populates_qc(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = (
            StratagemSurvey(edi_dir, csv, epsg=32649)
            .fit()
            .run_qc(include_skew=False)
        )
        assert sv.qc_ is not None
        assert len(sv.qc_.report_) == 5

    def test_remove_static_shift_runs(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        assert sv.remove_static_shift() is sv
        assert hasattr(sv, "_ss_corrector_")

    def test_drop_frequencies_runs(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        assert sv.drop_frequencies(fmin=10.0) is sv
        assert hasattr(sv, "_freq_filter_")

    def test_remove_noises_runs(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        assert sv.remove_noises() is sv
        assert hasattr(sv, "_noise_remover_")

    def test_export_writes_files(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        out = tmp_path / "export"
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit().export(out)
        assert hasattr(sv, "_writer_")
        assert sv._writer_.n_written_ == 5
        assert len(list(out.glob("*.edi"))) == 5

    def test_rename_writes_files(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        out_exp = tmp_path / "exported"
        out_ren = tmp_path / "renamed"
        sv = (
            StratagemSurvey(edi_dir, csv, epsg=32649)
            .fit()
            .export(out_exp)
            .rename(basename="T.", dst_path=out_ren)
        )
        assert hasattr(sv, "_renamer_")
        assert sv._renamer_.n_renamed_ == 5
        assert (out_ren / "T.000.edi").exists()
        assert (out_ren / "T.004.edi").exists()

    def test_fluent_full_pipeline(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        out = tmp_path / "final"
        ren = tmp_path / "renamed"
        sv = (
            StratagemSurvey(edi_dir, csv, epsg=32649)
            .fit()
            .remove_static_shift()
            .drop_frequencies(fmin=1.0)
            .remove_noises()
            .export(out)
            .rename(basename="T.", dst_path=ren)
        )
        assert sv._renamer_.n_renamed_ == 5

    def test_summary_is_string(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        s = sv.summary()
        assert isinstance(s, str)
        assert "StratagemSurvey" in s

    def test_coordinate_frame(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        frame = sv.coordinate_frame
        assert frame.shape == (5, 4)

    def test_processing_before_fit_raises(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649)
        with pytest.raises(Exception):
            sv.remove_static_shift()

    def test_accepts_list_of_edifile(self, tmp_path):
        """edi_dir may be a plain list of EDIFile, not just a directory."""
        from pycsamt.stratagem.survey import StratagemSurvey

        edis = _load_edis(tmp_path, n=5)
        csv = _make_coord_csv(tmp_path, n=5)
        sv = StratagemSurvey(edis, csv, epsg=32649).fit()
        assert sv.batch_ is None  # no directory was involved
        assert sv.n_stations_ == 5
        assert len(sv.edi_objects_) == 5

    def test_accepts_sites_object(self, tmp_path):
        """edi_dir may be an already-built Sites, per pycsamt.emtools convention."""
        from pycsamt.site.base import Sites
        from pycsamt.stratagem.survey import StratagemSurvey

        edis = _load_edis(tmp_path, n=5)
        csv = _make_coord_csv(tmp_path, n=5)
        sv = StratagemSurvey(Sites(edis), csv, epsg=32649).fit()
        assert sv.batch_ is None
        assert sv.n_stations_ == 5

    def test_sites_property_wraps_current_state(self, tmp_path):
        from pycsamt.site.base import Sites
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        sites = sv.sites_
        assert isinstance(sites, Sites)
        assert len(sites) == len(sv.edi_objects_)

    def test_sites_property_writes_via_site_api(self, tmp_path):
        edi_dir, csv = self._setup(tmp_path)
        from pycsamt.stratagem.survey import StratagemSurvey

        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        out = tmp_path / "sites_out"
        paths = sv.sites_.write(out, exist_ok=True)
        assert len(paths) == len(sv.edi_objects_)
        assert all(p.exists() for p in paths)

    def test_raw_dir_not_required(self, tmp_path):
        from pycsamt.stratagem.survey import StratagemSurvey

        edi_dir, csv = self._setup(tmp_path)
        sv = StratagemSurvey(edi_dir, csv, epsg=32649).fit()
        assert sv.raw_reader_ is None

"""Tests for stratagem.io (StratagemRawReader, EDIBatch)."""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pytest

from pycsamt.stratagem.io import (
    EDIBatch,
    StratagemRawReader,
    _edi_sort_key,
    _read_19col,
    _station_number,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_HEADER = dedent(
    """\
>HEAD
  DATAID="S{sid:02d}"
  ACQBY=AFEDIG
  LAT=0:00:00.00
  LONG=0:00:00.00
  ELEV=0

>INFO
  MAXINFO=999

>=DEFINEMEAS
  MAXCHAN=5
  UNITS=M
  REFTYPE=CART

>=MTSECT
  SECTID=S{sid:02d}
  NFREQ=2
  HX=1.001
  HY=2.001

>!****FREQUENCIES****!
>FREQ  //2
   1.000000E+02   1.000000E+01

>ZROT  //2
   0.000000E+00   0.000000E+00

>ZXX.VAR  //2
   1.0E+32  1.0E+32

>ZXY.VAR  //2
   1.0E+32  1.0E+32

>ZYX.VAR  //2
   1.0E+32  1.0E+32

>ZYY.VAR  //2
   1.0E+32  1.0E+32

>ZXX  //2
   0.0  0.0

>ZXY  //2
   1.0  2.0

>ZYX  //2
  -1.0 -2.0

>ZYY  //2
   0.0  0.0

>END
"""
)

_RAW_19COL = dedent(
    """\
 1.130e+001 2.930e+000 2.400e+001  3.728e+001  2.152e+002 -5.546e-001 -1.695e+002  1.336e+002  3.336e+004 -2.382e+002 -3.418e+003  3.016e+001  1.338e+002  2.845e+001 -1.052e+002 -2.512e+002 -1.384e+004 -2.318e+002  1.866e+004
 1.250e+001 2.930e+000 0.000e+000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.500e+001 2.930e+000 2.100e+001  2.248e+001  3.369e+002 -5.041e-002 -8.960e+001  1.700e+002  2.217e+004 -3.005e+002  5.244e+002  2.076e+001  1.842e+002  2.233e+001 -9.249e+001 -3.353e+002 -1.111e+004 -3.027e+002  1.656e+004
"""
)


def _make_edi_dir(tmp_dir: Path, n: int = 3) -> Path:
    d = tmp_dir / "edis"
    d.mkdir()
    for i in range(n):
        (d / f"Z2HX{i + 1:03d}.edi").write_text(_HEADER.format(sid=i), encoding="utf-8")
    return d


def _make_raw_dir(tmp_dir: Path, n: int = 3, n_rows: int = 3) -> Path:
    d = tmp_dir / "raw"
    d.mkdir()
    for i in range(n):
        (d / f"X2HX.{i + 1:03d}").write_text(_RAW_19COL, encoding="utf-8")
    return d


# ---------------------------------------------------------------------------
# _station_number
# ---------------------------------------------------------------------------


class TestStationNumber:
    def test_three_digit_suffix(self):
        assert _station_number("X2HX.001") == 1

    def test_stem_with_no_digits(self):
        assert _station_number("SENSORS") == 0

    def test_stem_with_middle_digits(self):
        assert _station_number("X2HX.087") == 87


# ---------------------------------------------------------------------------
# _edi_sort_key
# ---------------------------------------------------------------------------


class TestEdiSortKey:
    def test_orders_by_number(self):
        paths = [
            Path("Z2HX010.edi"),
            Path("Z2HX002.edi"),
            Path("Z2HX001.edi"),
        ]
        sorted_paths = sorted(paths, key=_edi_sort_key)
        assert [p.stem for p in sorted_paths] == [
            "Z2HX001",
            "Z2HX002",
            "Z2HX010",
        ]


# ---------------------------------------------------------------------------
# _read_19col
# ---------------------------------------------------------------------------


class TestRead19Col:
    def test_parses_valid_rows(self, tmp_path):
        f = tmp_path / "X.001"
        f.write_text(_RAW_19COL, encoding="utf-8")
        mat = _read_19col(f)
        assert mat.shape == (3, 19)
        assert mat[0, 0] == pytest.approx(11.3)
        assert mat[1, 2] == pytest.approx(0.0)  # zero stack row

    def test_zero_stack_detected(self, tmp_path):
        f = tmp_path / "X.001"
        f.write_text(_RAW_19COL, encoding="utf-8")
        mat = _read_19col(f)
        stacks = mat[:, 2].astype(int)
        assert stacks[1] == 0  # second row has stack = 0
        assert stacks[0] == 24  # first row: 24 stacks

    def test_empty_file_returns_zeros(self, tmp_path):
        f = tmp_path / "X.001"
        f.write_text("", encoding="utf-8")
        mat = _read_19col(f)
        assert mat.shape == (0, 19)


# ---------------------------------------------------------------------------
# StratagemRawReader
# ---------------------------------------------------------------------------


class TestStratagemRawReader:
    def test_fit_returns_self(self, tmp_path):
        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw")
        assert rdr.fit() is rdr

    def test_stations_and_shapes(self, tmp_path):
        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        assert rdr.n_stations_ == 3
        assert rdr.n_freqs_ == 3
        assert rdr.snr_mask_.shape == (3, 3)
        assert rdr.stack_counts_.shape == (3, 3)

    def test_snr_mask_correct(self, tmp_path):
        _make_raw_dir(tmp_path, n=2)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        # row 1 (index 1) has stack=0 → snr_mask False
        assert not rdr.snr_mask_[0, 1]  # second row of file has stack=0
        assert rdr.snr_mask_[0, 0]  # first row has stack=24

    def test_missing_dir_raises(self):
        with pytest.raises(Exception):
            StratagemRawReader("/does/not/exist").fit()

    def test_no_files_raises(self, tmp_path):
        with pytest.raises(Exception):
            StratagemRawReader(tmp_path).fit()

    def test_station_coverage_range(self, tmp_path):
        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        cov = rdr.station_coverage()
        assert 0.0 <= cov <= 1.0

    def test_usable_freq_counts_shape(self, tmp_path):
        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        counts = rdr.usable_freq_counts()
        assert counts.shape == (3,)


# ---------------------------------------------------------------------------
# EDIBatch
# ---------------------------------------------------------------------------


class TestEDIBatch:
    def test_fit_returns_self(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis")
        assert batch.fit() is batch

    def test_n_stations(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        assert batch.n_stations_ == 3

    def test_len_protocol(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        assert len(batch) == 3

    def test_iter_protocol(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        assert len(list(batch)) == 3

    def test_getitem(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        edi = batch[0]
        assert edi is not None

    def test_station_names(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        names = batch.station_names()
        assert len(names) == 3

    def test_natural_sort_order(self, tmp_path):
        d = tmp_path / "edis"
        d.mkdir()
        for sid in [10, 2, 1]:
            (d / f"Z2HX{sid:03d}.edi").write_text(
                _HEADER.format(sid=sid), encoding="utf-8"
            )
        batch = EDIBatch(d).fit()
        # should be sorted: 001, 002, 010
        stems = [p.stem for p in batch.edi_paths_]
        assert stems == sorted(stems, key=lambda s: int(s[-3:]))

    def test_missing_dir_raises(self):
        with pytest.raises(Exception):
            EDIBatch("/does/not/exist").fit()

    def test_empty_dir_raises(self, tmp_path):
        (tmp_path / "empty").mkdir()
        with pytest.raises(Exception):
            EDIBatch(tmp_path / "empty").fit()

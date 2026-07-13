"""Phase 4 tests: StratagemRawReader enhancements and station-matching fix."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# shared test fixtures
# ---------------------------------------------------------------------------

_RAW_19COL = textwrap.dedent(
    """\
 1.130e+001 2.930e+000 2.400e+001  3.728e+001  2.152e+002 -5.546e-001 -1.695e+002  1.336e+002  3.336e+004 -2.382e+002 -3.418e+003  3.016e+001  1.338e+002  2.845e+001 -1.052e+002 -2.512e+002 -1.384e+004 -2.318e+002  1.866e+004
 1.250e+001 2.930e+000 0.000e+000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.500e+001 2.930e+000 2.100e+001  2.248e+001  3.369e+002 -5.041e-002 -8.960e+001  1.700e+002  2.217e+004 -3.005e+002  5.244e+002  2.076e+001  1.842e+002  2.233e+001 -9.249e+001 -3.353e+002 -1.111e+004 -3.027e+002  1.656e+004
"""
)

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
>=MTSECT
  SECTID=S{sid:02d}
  NFREQ=3
  HX=1.001
  HY=2.001
>!****FREQUENCIES****!
>FREQ  //3
   1.130e+001   1.250e+001   1.500e+001
>ZROT  //3
   0.  0.  0.
>ZXX  //3
   0.  0.  0.
>ZXY  //3
   1.0  2.0  3.0
>ZYX  //3
  -1.0 -2.0 -3.0
>ZYY  //3
   0.  0.  0.
>ZXX.VAR  //3
   0.01 0.01 0.01
>ZXY.VAR  //3
   0.01 0.01 0.01
>ZYX.VAR  //3
   0.01 0.01 0.01
>ZYY.VAR  //3
   0.01 0.01 0.01
>END
"""
)


def _make_raw_dir(tmp_path: Path, n: int = 4, prefix: str = "X2HX") -> Path:
    """Create n raw station files starting from station 1 (X2HX.001, ...)."""
    d = tmp_path / "raw"
    d.mkdir(exist_ok=True)
    for i in range(n):
        (d / f"{prefix}.{i + 1:03d}").write_text(_RAW_19COL, encoding="utf-8")
    return d


def _make_edi_dir(tmp_path: Path, n: int = 3, start: int = 2) -> Path:
    """Create n EDI files starting from station `start` (e.g., Z2HX002.edi)."""
    d = tmp_path / "edis"
    d.mkdir(exist_ok=True)
    for i in range(n):
        sid = start + i
        lat = 25.77 + i * 0.01
        lon = 109.60 + i * 0.02
        (d / f"Z2HX{sid:03d}.edi").write_text(
            _EDI_TMPL.format(sid=sid, lat=lat, lon=lon), encoding="utf-8"
        )
    return d


# ---------------------------------------------------------------------------
# StratagemRawReader: stations_ and station_numbers_
# ---------------------------------------------------------------------------


class TestStationsAttribute:
    def test_stations_stores_file_names(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        # must be file names like 'X2HX.001', not stems like 'X2HX'
        assert all("." in s for s in rdr.stations_)
        assert rdr.stations_[0] == "X2HX.001"
        assert rdr.stations_[2] == "X2HX.003"

    def test_station_numbers_correct(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        np.testing.assert_array_equal(rdr.station_numbers_, [1, 2, 3])

    def test_station_numbers_dtype(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        assert rdr.station_numbers_.dtype == np.int32

    def test_all_stations_unique(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=4)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        assert len(set(rdr.stations_)) == 4


# ---------------------------------------------------------------------------
# StratagemRawReader: match_to_edis
# ---------------------------------------------------------------------------


class TestMatchToEdis:
    def test_exact_match(self, tmp_path):
        """Raw 1-4, EDI starting at 2: EDI[0]=station2 → raw[1]."""
        from pycsamt.stratagem.io import (
            EDIBatch,
            StratagemRawReader,
        )

        _make_raw_dir(tmp_path, n=4)
        _make_edi_dir(tmp_path, n=3, start=2)  # Z2HX002..Z2HX004
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        batch = EDIBatch(tmp_path / "edis").fit()

        mapping = rdr.match_to_edis(batch.edi_objects_)
        # EDI[0] = Z2HX002 = station 2 → raw index 1 (X2HX.002)
        assert mapping[0] == 1
        # EDI[1] = Z2HX003 = station 3 → raw index 2 (X2HX.003)
        assert mapping[1] == 2

    def test_no_match_missing_station(self, tmp_path):
        """EDI station 5 has no raw file when raw only goes to 4."""
        from pycsamt.stratagem.io import (
            EDIBatch,
            StratagemRawReader,
        )

        _make_raw_dir(tmp_path, n=4)
        _make_edi_dir(tmp_path, n=1, start=5)  # Z2HX005.edi — no raw file
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        batch = EDIBatch(tmp_path / "edis").fit()

        mapping = rdr.match_to_edis(batch.edi_objects_)
        assert 0 not in mapping  # no match for station 5

    def test_returns_dict(self, tmp_path):
        from pycsamt.stratagem.io import (
            EDIBatch,
            StratagemRawReader,
        )

        _make_raw_dir(tmp_path, n=3)
        _make_edi_dir(tmp_path, n=2, start=2)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        batch = EDIBatch(tmp_path / "edis").fit()
        m = rdr.match_to_edis(batch.edi_objects_)
        assert isinstance(m, dict)

    def test_before_fit_raises(self):
        from pycsamt.stratagem.io import StratagemRawReader

        with pytest.raises(Exception):
            StratagemRawReader().match_to_edis([])


# ---------------------------------------------------------------------------
# StratagemRawReader: station_frame
# ---------------------------------------------------------------------------


class TestStationFrame:
    def test_returns_dataframe(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        df = rdr.station_frame()
        assert isinstance(df, pd.DataFrame)

    def test_row_count(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=4)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        assert len(rdr.station_frame()) == 4

    def test_required_columns(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        df = StratagemRawReader(tmp_path / "raw").fit().station_frame()
        for col in ("station", "station_number", "usable_freqs", "coverage"):
            assert col in df.columns

    def test_coverage_in_range(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        df = StratagemRawReader(tmp_path / "raw").fit().station_frame()
        assert df["coverage"].between(0, 1).all()

    def test_before_fit_raises(self):
        from pycsamt.stratagem.io import StratagemRawReader

        with pytest.raises(Exception):
            StratagemRawReader().station_frame()


# ---------------------------------------------------------------------------
# StratagemRawReader: freq_frame
# ---------------------------------------------------------------------------


class TestFreqFrame:
    def test_returns_dataframe(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        df = StratagemRawReader(tmp_path / "raw").fit().freq_frame()
        assert isinstance(df, pd.DataFrame)

    def test_row_count_matches_n_freqs(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        assert len(rdr.freq_frame()) == rdr.n_freqs_

    def test_required_columns(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        df = StratagemRawReader(tmp_path / "raw").fit().freq_frame()
        for col in ("freq_hz", "stations_ok", "frac_ok"):
            assert col in df.columns

    def test_before_fit_raises(self):
        from pycsamt.stratagem.io import StratagemRawReader

        with pytest.raises(Exception):
            StratagemRawReader().freq_frame()


# ---------------------------------------------------------------------------
# StratagemRawReader: stack_audit
# ---------------------------------------------------------------------------


class TestStackAudit:
    def test_returns_dataframe(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=3)
        df = StratagemRawReader(tmp_path / "raw").fit().stack_audit()
        assert isinstance(df, pd.DataFrame)

    def test_shape(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=3)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        df = rdr.stack_audit()
        assert df.shape == (rdr.n_stations_, rdr.n_freqs_)

    def test_index_is_station_names(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        rdr = StratagemRawReader(tmp_path / "raw").fit()
        df = rdr.stack_audit()
        assert list(df.index) == rdr.stations_

    def test_values_non_negative(self, tmp_path):
        from pycsamt.stratagem.io import StratagemRawReader

        _make_raw_dir(tmp_path, n=2)
        df = StratagemRawReader(tmp_path / "raw").fit().stack_audit()
        assert (df.values >= 0).all()

    def test_before_fit_raises(self):
        from pycsamt.stratagem.io import StratagemRawReader

        with pytest.raises(Exception):
            StratagemRawReader().stack_audit()


# ---------------------------------------------------------------------------
# FrequencyFilter: correct station-number alignment
# ---------------------------------------------------------------------------


class TestFrequencyFilterAlignment:
    """Verify the hardware mask is applied to the correct EDI stations."""

    def test_alignment_by_station_number(self, tmp_path):
        """Raw stations 1-4, EDIs from station 2.
        The hardware mask for EDI[0] (station 2) should come from raw[1], NOT raw[0].
        """
        from pycsamt.stratagem.io import (
            EDIBatch,
            StratagemRawReader,
        )
        from pycsamt.stratagem.qc import FrequencyFilter

        # raw: 4 stations (1-4), station 1 has ALL bad freqs
        d = tmp_path / "raw"
        d.mkdir()
        all_bad = textwrap.dedent(
            """\
 1.130e+001 2.930e+000 0.000e+000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
 1.250e+001 2.930e+000 0.000e+000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
 1.500e+001 2.930e+000 0.000e+000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
"""
        )
        (d / "X2HX.001").write_text(
            all_bad, encoding="utf-8"
        )  # station 1 - all bad
        (d / "X2HX.002").write_text(
            _RAW_19COL, encoding="utf-8"
        )  # station 2 - has valid data
        (d / "X2HX.003").write_text(_RAW_19COL, encoding="utf-8")
        (d / "X2HX.004").write_text(_RAW_19COL, encoding="utf-8")

        _make_edi_dir(tmp_path, n=3, start=2)  # Z2HX002..Z2HX004

        rdr = StratagemRawReader(tmp_path / "raw").fit()
        batch = EDIBatch(tmp_path / "edis").fit()
        edis = batch.edi_objects_

        # Confirm station-number alignment:
        # EDI[0] = Z2HX002 (station 2) → raw[1] = X2HX.002 (has valid data)
        # NOT raw[0] = X2HX.001 (all bad)
        mapping = rdr.match_to_edis(edis)
        assert mapping[0] == 1, (
            f"EDI[0] should map to raw[1], got {mapping[0]}"
        )

        # FrequencyFilter should apply the raw[1] mask to EDI[0].
        # Disable the incoherence step (snr_thresh/min_frac=0) so the test
        # only exercises the hardware-mask alignment logic; it avoids false
        # failures from all-zero Z values in the minimal test EDI template
        # (zero Z → SNR=0 → incoherence filter would mask everything).
        ff = FrequencyFilter(
            use_hardware_mask=True,
            snr_thresh=0.0,
            min_frac=0.0,
        ).fit(edis, raw_reader=rdr)
        z0 = ff.edi_objects_[0].Z.z
        if z0 is not None:
            # EDI[0] maps to raw station 2 (partial valid: snr=[T,F,T]).
            # After hardware masking, rows 0 and 2 are kept (finite = 0+0j),
            # row 1 is NaN.  If alignment were wrong (station 1 = all-bad),
            # all three rows would be NaN and .any() would be False.
            assert np.isfinite(z0).any(), (
                "All Z values are NaN — hardware mask alignment used station 1 "
                "(all-bad) instead of station 2 (partial valid)."
            )

"""Tests for stratagem.gis_correct (CoordinateInjector, StationLocator)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.stratagem.gis_correct import (
    CoordinateInjector,
    StationLocator,
    _detect_coord_cols,
)
from pycsamt.stratagem.io import EDIBatch

# ---------------------------------------------------------------------------
# minimal EDI template (reused from test_io helpers)
# ---------------------------------------------------------------------------

_EDI_TMPL = textwrap.dedent(
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

>ZXY  //2
   1.0  2.0

>ZYX  //2
  -1.0 -2.0

>ZXX  //2
   0.0  0.0

>ZYY  //2
   0.0  0.0

>ZXX.VAR  //2
   1.0E+32  1.0E+32

>ZXY.VAR  //2
   1.0E+32  1.0E+32

>ZYX.VAR  //2
   1.0E+32  1.0E+32

>ZYY.VAR  //2
   1.0E+32  1.0E+32

>END
"""
)


def _make_edi_dir(tmp_path: Path, n: int = 3) -> Path:
    d = tmp_path / "edis"
    d.mkdir(exist_ok=True)
    for i in range(n):
        (d / f"Z2HX{i + 1:03d}.edi").write_text(
            _EDI_TMPL.format(sid=i), encoding="utf-8"
        )
    return d


def _make_csv(tmp_path: Path, n: int = 3) -> Path:
    """Create a minimal coordinate CSV matching the Stratagem survey format.

    Uses the misleading ``longitude`` / ``latitude`` column names that
    actually store UTM Northing / Easting values (as in the real survey).
    Coordinates are in WGS84 UTM Zone 49N so that project_point_utm2ll
    can round-trip them.
    """
    rows = []
    for i in range(n):
        # synthetic UTM 49N coords inside China (roughly Guangxi)
        easting = 362_589.0 + i * 20.0
        northing = 2_850_835.0 - i * 20.0
        rows.append(
            {
                "stations": f"K0+{i * 20:06.4f}",
                "longitude": northing,  # intentionally mislabeled
                "latitude": easting,  # intentionally mislabeled
                "elev": 261.94 + i,
                "step": i * 20.0,
            }
        )
    df = pd.DataFrame(rows)
    csv_path = tmp_path / "coords.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


# ---------------------------------------------------------------------------
# _detect_coord_cols
# ---------------------------------------------------------------------------


class TestDetectCoordCols:
    def test_explicit_override(self):
        df = pd.DataFrame({"E": [1.0], "N": [2.0], "elev": [0.0]})
        e, n = _detect_coord_cols(df, "E", "N", exclude={"elev"})
        assert e == "E" and n == "N"

    def test_auto_detect_by_magnitude(self):
        df = pd.DataFrame(
            {
                "longitude": [2_850_835.0, 2_850_815.0],  # northing (big)
                "latitude": [362_589.0, 362_609.0],  # easting (small)
                "elev": [261.9, 262.0],
            }
        )
        e, n = _detect_coord_cols(df, None, None, exclude={"elev", "step"})
        assert e == "latitude"  # smaller median → easting
        assert n == "longitude"  # larger median  → northing

    def test_raises_when_too_few_cols(self):
        df = pd.DataFrame({"elev": [1.0, 2.0]})
        with pytest.raises(Exception):
            _detect_coord_cols(df, None, None, exclude={"elev"})

    def test_raises_when_too_many_cols(self):
        """>2 numeric candidates is ambiguous — must raise, not guess.

        Regression: a station-index column (small median) used to get
        silently picked as 'easting' instead of the real easting column.
        """
        df = pd.DataFrame(
            {
                "station_num": [1, 2, 3],  # small median — not a coordinate
                "longitude": [2_850_835.0, 2_850_815.0, 2_850_795.0],
                "latitude": [362_589.0, 362_609.0, 362_629.0],
                "elev": [261.9, 262.0, 262.1],
            }
        )
        with pytest.raises(Exception):
            _detect_coord_cols(df, None, None, exclude={"elev"})

    def test_explicit_cols_bypass_ambiguity(self):
        """Explicit easting_col/northing_col still work with >2 candidates."""
        df = pd.DataFrame(
            {
                "station_num": [1, 2, 3],
                "longitude": [2_850_835.0, 2_850_815.0, 2_850_795.0],
                "latitude": [362_589.0, 362_609.0, 362_629.0],
                "elev": [261.9, 262.0, 262.1],
            }
        )
        e, n = _detect_coord_cols(
            df, "latitude", "longitude", exclude={"elev"}
        )
        assert e == "latitude" and n == "longitude"


# ---------------------------------------------------------------------------
# StationLocator
# ---------------------------------------------------------------------------


class TestStationLocator:
    def _dummy_edis(self, n: int):
        # minimal objects with no path; StationLocator only needs len(edi_objects)
        class _E:
            path = None
            station = None

        return [_E() for _ in range(n)]

    def test_forward_order(self):
        n = 5
        lats = np.linspace(25.0, 25.5, n)
        lons = np.linspace(109.5, 110.0, n)
        loc = StationLocator(order="forward").fit(
            self._dummy_edis(n), lats, lons
        )
        assert loc.index_map_ == list(range(n))
        assert loc.reversed_ is False

    def test_reversed_order(self):
        n = 5
        lats = np.linspace(25.0, 25.5, n)
        lons = np.linspace(109.5, 110.0, n)
        loc = StationLocator(order="reversed").fit(
            self._dummy_edis(n), lats, lons
        )
        assert loc.index_map_ == list(range(n - 1, -1, -1))
        assert loc.reversed_ is True

    def test_mapping_order(self):
        n = 4
        lats = np.zeros(n)
        lons = np.zeros(n)
        mapping = [3, 2, 1, 0]
        loc = StationLocator(order="mapping", mapping=mapping).fit(
            self._dummy_edis(n), lats, lons
        )
        assert loc.index_map_ == mapping

    def test_mapping_wrong_length_raises(self):
        n = 4
        with pytest.raises(Exception):
            StationLocator(order="mapping", mapping=[0, 1]).fit(
                self._dummy_edis(n),
                np.zeros(n),
                np.zeros(n),
            )

    def test_auto_returns_valid_map(self):
        n = 5
        lats = np.linspace(25.0, 25.5, n)
        lons = np.zeros(n)
        loc = StationLocator(order="auto").fit(
            self._dummy_edis(n), lats, lons
        )
        assert len(loc.index_map_) == n
        assert set(loc.index_map_) == set(range(n))


# ---------------------------------------------------------------------------
# CoordinateInjector
# ---------------------------------------------------------------------------


class TestCoordinateInjector:
    def test_fit_returns_self(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=3)
        inj = CoordinateInjector(epsg=32649, utm_zone="49N")
        assert inj.fit(batch, csv) is inj

    def test_latitudes_shape(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=3)
        inj = CoordinateInjector(epsg=32649, utm_zone="49N").fit(batch, csv)
        assert inj.latitudes_.shape == (3,)
        assert inj.longitudes_.shape == (3,)
        assert inj.elevations_.shape == (3,)

    def test_coordinates_in_valid_range(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=3)
        inj = CoordinateInjector(epsg=32649, utm_zone="49N").fit(batch, csv)
        # should be somewhere in southern China
        assert np.all(inj.latitudes_ > 20.0)
        assert np.all(inj.latitudes_ < 35.0)
        assert np.all(inj.longitudes_ > 100.0)
        assert np.all(inj.longitudes_ < 120.0)

    def test_edi_head_updated(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=3)
        inj = CoordinateInjector(epsg=32649, utm_zone="49N").fit(batch, csv)
        head = inj.edi_objects_[0].get_section("head")
        assert head is not None
        assert head.lat is not None and head.lat != 0.0
        assert head.lon is not None and head.lon != 0.0

    def test_mismatch_raises(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=5)  # 5 rows, 3 EDIs → mismatch
        with pytest.raises(Exception):
            CoordinateInjector(epsg=32649, utm_zone="49N").fit(batch, csv)

    def test_coordinate_frame_shape(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=3)
        inj = CoordinateInjector(epsg=32649, utm_zone="49N").fit(batch, csv)
        frame = inj.coordinate_frame()
        assert frame.shape == (3, 4)
        assert list(frame.columns) == [
            "station",
            "latitude",
            "longitude",
            "elevation",
        ]

    def test_export_creates_files(self, tmp_path):
        _make_edi_dir(tmp_path, n=3)
        batch = EDIBatch(tmp_path / "edis").fit()
        csv = _make_csv(tmp_path, n=3)
        inj = CoordinateInjector(epsg=32649, utm_zone="49N").fit(batch, csv)
        out = tmp_path / "out"
        paths = inj.export(out)
        assert len(paths) == 3
        assert all(p.exists() for p in paths)

    def test_export_before_fit_raises(self):
        inj = CoordinateInjector()
        with pytest.raises(Exception):
            inj.export("/tmp/nowhere")

    def test_explicit_column_names(self, tmp_path):
        """Explicit easting_col / northing_col bypass auto-detection."""
        _make_edi_dir(tmp_path, n=2)
        batch = EDIBatch(tmp_path / "edis").fit()

        df = pd.DataFrame(
            {
                "stations": ["A", "B"],
                "E": [362_589.0, 362_609.0],
                "N": [2_850_835.0, 2_850_815.0],
                "elev": [261.9, 262.0],
            }
        )
        csv = tmp_path / "custom.csv"
        df.to_csv(csv, index=False)

        inj = CoordinateInjector(epsg=32649, utm_zone="49N").fit(
            batch, csv, easting_col="E", northing_col="N"
        )
        assert inj.latitudes_.shape == (2,)

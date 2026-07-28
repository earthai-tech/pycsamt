# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Format-parser tests for :mod:`pycsamt.tdem.io` and the
:class:`~pycsamt.tdem.reader.TEMReader` auto-detection front end.

Each synthetic fixture is the minimal document shown in the parser's
own docstring; the bundled ``data/TEMAVG/JIANGSU`` files provide the
real-file smoke paths.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.tdem import io as tio
from pycsamt.tdem.reader import TEMReader, _detect_format

_DATA = Path(__file__).resolve().parents[3] / "data" / "TEMAVG" / "JIANGSU"


# ── synthetic documents ──────────────────────────────────────────────────────

GEOSOFT_DAT = """\
/ Geosoft TEM export
/ Current: 8.0
/ LoopSide: 200.0
/ GateTimes(ms): 0.021 0.037 0.065 0.090
/ DataUnit: nV/Am2
X        Y        ELEV    G1       G2       G3       G4
L100
1000.0   2000.0   10.0    1.24e4   8.91e3   5.23e3   3.12e3
2000.0   2000.0   12.0    9.87e3   6.54e3   4.32e3   2.61e3
"""

AMIRA_TEM = """\
; EMIT TEM data - site survey
TRANSMITTER_AREA   40000
TRANSMITTER_CURRENT  8.0
TRANSMITTER_TURNS    1
RECEIVER_AREA  1.0
RECEIVER_TURNS 1
GATE_TIMES_UNIT ms
GATE_TIMES  0.021 0.037 0.065 0.090 0.135 0.200
DATA_UNIT   nV/Am2

DATA
; StnName  X         Y         Elev    D1        D2 ...
S001  1000.0  2000.0    0.0  1.24e+04  8.91e+03  5.23e+03  3.12e+03  2.0e+03  1.1e+03
S002  2000.0  2000.0    5.0  9.87e+03  6.54e+03  4.32e+03  2.61e+03  1.7e+03  0.9e+03
END DATA
"""

ZONGE_GDP = """\
* Zonge GDP-32 TEM sounding
* Station: S001
* X: 1000.0   Y: 2000.0   Elev: 0.0
* Current: 8.00 A
* TxArea: 40000 m^2
* GateTimes(ms): 0.021 0.037 0.065 0.090
* DataUnit: nV/Am2
Win  Time(ms)    Hz(nV/Am2)   Error
  1    0.021     1.2400e+04   5.0e+02
  2    0.037     8.9100e+03   4.0e+02
  3    0.065     5.2300e+03   2.5e+02
  4    0.090     3.1200e+03   1.8e+02
* Station: S002
* X: 1200.0   Y: 2000.0   Elev: 2.0
* Current: 8.00 A
* TxArea: 40000 m^2
* GateTimes(ms): 0.021 0.037 0.065 0.090
* DataUnit: nV/Am2
Win  Time(ms)    Hz(nV/Am2)   Error
  1    0.021     1.1000e+04   5.0e+02
  2    0.037     8.0000e+03   4.0e+02
  3    0.065     5.0000e+03   2.5e+02
  4    0.090     3.0000e+03   1.8e+02
"""

WALKTEM = """\
/ WalkTEM data - Aarhus Workbench export
/
GENERALHEADER
WalkTEM System
END_GENERALHEADER

TXOPERATION
LoopSize    40.000
Current      6.60
TurnsNumber  1
Waveform     Full_Waveform
END_TXOPERATION

RXOPERATION
CoilMoment  0.01
END_RXOPERATION

GATESET
NumberOfGates 5
GateTimes  0.008 0.015 0.025 0.043 0.072
GateSigns  1 1 1 1 1
END_GATESET

DATA
! Line  Fid   X        Y       Elev   d1      d2      d3      d4      d5
100  1000  1000.0  2000.0   0.0  1.5e-7  9.2e-8  5.1e-8  2.9e-8  1.4e-8
100  2000  1200.0  2000.0   0.0  1.3e-7  8.1e-8  4.6e-8  2.5e-8  1.2e-8
END_DATA
"""

XYZ_SIMPLE = """\
# time(ms)   dBdt(nT/s)   err
0.021   1.24e2   5.0
0.037   0.89e2   4.0
0.065   0.52e2   2.5
0.090   0.31e2   1.8
"""

XYZ_MULTI = """\
# site  time(s)  data(SI)
A1  2.1e-5  1.2e-7
A1  3.7e-5  0.9e-7
A2  2.1e-5  1.4e-7
A2  3.7e-5  1.0e-7
"""


@pytest.fixture()
def fmt_files(tmp_path: Path) -> dict[str, Path]:
    files = {
        "geosoft": tmp_path / "survey.dat",
        "amira": tmp_path / "survey_amira.tem",
        "zonge": tmp_path / "gdp.avg",
        "walkttem": tmp_path / "walk.tem",
        "xyz": tmp_path / "sounding.txt",
        "multi": tmp_path / "multi.txt",
    }
    files["geosoft"].write_text(GEOSOFT_DAT, encoding="utf-8")
    files["amira"].write_text(AMIRA_TEM, encoding="utf-8")
    files["zonge"].write_text(ZONGE_GDP, encoding="utf-8")
    files["walkttem"].write_text(WALKTEM, encoding="utf-8")
    files["xyz"].write_text(XYZ_SIMPLE, encoding="utf-8")
    files["multi"].write_text(XYZ_MULTI, encoding="utf-8")
    return files


# ── read_xyz / read_xyz_multisite ────────────────────────────────────────────


def test_read_xyz_units_and_error_column(fmt_files):
    snd = tio.read_xyz(
        fmt_files["xyz"],
        current=8.0,
        loop_side=200.0,
        time_col=0,
        data_col=1,
        error_col=2,
        time_unit="ms",
        data_unit="nT/s",
        station_name="S1",
    )
    assert snd.time_gates.shape == (4,)
    # 0.021 ms -> 2.1e-5 s
    assert snd.time_gates[0] == pytest.approx(2.1e-5)
    # 124 nT/s -> 1.24e-7 T/s
    assert snd.data[0] == pytest.approx(1.24e-7)
    assert snd.error is not None
    assert snd.current == pytest.approx(8.0)
    assert snd.tx_area == pytest.approx(200.0**2)


def test_read_xyz_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        tio.read_xyz(tmp_path / "nope.txt", current=1.0, tx_area=1.0)


def test_read_xyz_multisite_groups_by_site(fmt_files):
    out = tio.read_xyz_multisite(
        fmt_files["multi"],
        site_col=0,
        time_col=1,
        data_col=2,
        current=5.0,
        tx_area=1600.0,
    )
    assert set(out) == {"A1", "A2"}
    assert out["A1"].time_gates.shape == (2,)
    assert out["A2"].station_name == "A2"


# ── geosoft ──────────────────────────────────────────────────────────────────


def test_read_geosoft_dat_parses_header_metadata(fmt_files):
    soundings = tio.read_geosoft_dat(fmt_files["geosoft"])
    assert len(soundings) == 2
    s0 = soundings[0]
    assert s0.time_gates.shape == (4,)
    assert s0.time_gates[0] == pytest.approx(2.1e-5)  # ms -> s
    assert s0.current == pytest.approx(8.0)
    assert s0.tx_area == pytest.approx(200.0**2)
    # nV/Am2 -> SI conversion happened (values far below raw 1.24e4)
    assert s0.data[0] < 1.0
    assert np.all(np.diff(s0.data) < 0)  # decaying transient


def test_read_geosoft_dat_missing_gate_times(tmp_path):
    p = tmp_path / "nogates.dat"
    p.write_text(
        "/ Current: 8.0\n/ LoopSide: 200.0\nX Y ELEV G1 G2\n1 2 0 1e4 5e3\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        tio.read_geosoft_dat(p)


# ── amira ────────────────────────────────────────────────────────────────────


def test_read_amira_keyword_block(fmt_files):
    soundings = tio.read_amira(fmt_files["amira"])
    assert len(soundings) == 2
    s0 = soundings[0]
    assert s0.station_name == "S001"
    assert s0.current == pytest.approx(8.0)
    assert s0.tx_area == pytest.approx(40000.0)
    assert s0.time_gates.shape == (6,)
    assert s0.time_gates[-1] == pytest.approx(2.0e-4)  # 0.2 ms


def test_read_amira_kwarg_overrides(fmt_files):
    soundings = tio.read_amira(fmt_files["amira"], current=4.0)
    assert soundings[0].current == pytest.approx(4.0)


# ── zonge GDP ────────────────────────────────────────────────────────────────


def test_read_zonge_multi_station_blocks(fmt_files):
    soundings = tio.read_zonge(fmt_files["zonge"])
    assert len(soundings) == 2
    names = [s.station_name for s in soundings]
    assert names == ["S001", "S002"]
    s0 = soundings[0]
    assert s0.time_gates.shape == (4,)
    assert s0.current == pytest.approx(8.0)
    assert s0.error is not None


# ── walkTEM ──────────────────────────────────────────────────────────────────


def test_read_walkttem_sections(fmt_files):
    soundings = tio.read_walkttem(fmt_files["walkttem"])
    assert len(soundings) == 2
    s0 = soundings[0]
    assert s0.current == pytest.approx(6.6)
    assert s0.tx_area == pytest.approx(40.0**2)
    assert s0.time_gates.shape == (5,)


# ── real bundled JIANGSU files ───────────────────────────────────────────────


@pytest.fixture(scope="module")
def jiangsu() -> Path:
    if not _DATA.exists():
        pytest.skip(f"TEMAVG dataset not found: {_DATA}")
    return _DATA


def test_read_temavg_real_file(jiangsu):
    obj = tio.read_temavg(jiangsu / "TEM100.AVG")
    assert obj is not None
    df = getattr(obj, "df", None)
    if df is not None:
        assert len(df) > 0


def test_read_tem_z_real_file(jiangsu):
    obj = tio.read_tem_z(jiangsu / "TEM100.Z")
    assert obj is not None


def test_read_tem_log_real_file(jiangsu):
    obj = tio.read_tem_log(jiangsu / "TEM100.LOG")
    assert obj is not None


# ── TEMReader auto-detection ─────────────────────────────────────────────────


def test_detect_format_disambiguation(fmt_files, jiangsu):
    assert _detect_format(fmt_files["geosoft"]) == "geosoft"
    assert _detect_format(fmt_files["amira"]) == "amira"
    assert _detect_format(fmt_files["walkttem"]) == "walkttem"
    assert _detect_format(fmt_files["zonge"]) == "zonge"
    assert _detect_format(jiangsu / "TEM100.AVG") == "temavg"
    assert _detect_format(fmt_files["xyz"]) == "xyz"


def test_temreader_reads_all_synthetic_formats(fmt_files):
    reader = TEMReader()
    geo = reader.read(fmt_files["geosoft"])
    assert len(geo) == 2
    am = reader.read(fmt_files["amira"])
    assert len(am) == 2
    zg = reader.read(fmt_files["zonge"])
    assert len(zg) == 2
    wt = reader.read(fmt_files["walkttem"])
    assert len(wt) == 2


def test_temreader_explicit_fmt_and_kwargs(fmt_files):
    reader = TEMReader()
    snd = reader.read(
        fmt_files["xyz"],
        fmt="xyz",
        current=8.0,
        tx_area=100.0,
        time_unit="ms",
        data_unit="nT/s",
        error_col=2,
    )
    assert snd.time_gates.shape == (4,)


def test_temreader_unknown_format_raises(fmt_files):
    reader = TEMReader()
    with pytest.raises((KeyError, ValueError)):
        reader.read(fmt_files["xyz"], fmt="not-a-format")


def test_temreader_lists_formats():
    fmts = TEMReader().formats
    assert "geosoft" in fmts
    assert "walkttem" in fmts


def test_temreader_named_methods_and_cache(fmt_files):
    reader = TEMReader(store=True)
    am = reader.read_amira(fmt_files["amira"])
    assert len(am) == 2
    geo = reader.read_geosoft_dat(fmt_files["geosoft"])
    assert len(geo) == 2
    # both results cached under their paths
    assert len(reader.results) == 2
    reader.clear()
    assert not reader.results


def test_temreader_named_zonge_and_walkttem(fmt_files):
    reader = TEMReader()
    assert len(reader.read_zonge(fmt_files["zonge"])) == 2
    assert len(reader.read_walkttem(fmt_files["walkttem"])) == 2


def test_temreader_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        TEMReader().read("does/not/exist.avg")

# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from math import isclose
from pathlib import Path

from pycsamt.seg.meas import (
    DefineMeas,
    DefineMeasMixin,
    EMeasMixin,
    Emeasurement,
    HMeasMixin,
    Hmeasurement,
    MeasMixin,
)


def _write_sample_edi(p: Path) -> None:
    head = (
        ">HEAD\n\n"
        '  DATAID="DEMO"\n'
        "  ACQDATE=01/01/24\n"
        "  FILEDATE=01/01/24\n"
        "  LAT=0\n"
        "  LONG=0\n"
        "  EMPTY=1.0E+32\n\n"
    )
    info = ">INFO\n  MAXINFO=9\n\n"
    definemeas = (
        ">=DEFINEMEAS\n"
        "  MAXCHAN=7\n"
        "  MAXRUN=999\n"
        "  MAXMEAS=7\n"
        "  UNITS=M\n"
        "  REFTYPE=CART\n"
        "  REFLAT=-22:22:14.9\n"
        "  REFLONG=139:11:19.1\n"
        "  REFELEV=200\n"
        "  >HMEAS ID=251.025 CHTYPE=HX X=8.5 Y=8.5 AZM=0 "
        "ACQCHAN=CH3\n"
        "  >HMEAS ID=252.025 CHTYPE=HY X=-8.5 Y=8.5 AZM=90 "
        "ACQCHAN=CH4\n"
        "  >HMEAS ID=253.025 CHTYPE=HZ X=21.2 Y=-21.2 AZM=0 "
        "ACQCHAN=CH5\n"
        "  >EMEAS ID=254.025 CHTYPE=EX X=-50.0 Y=0.0 X2=50.0 "
        "Y2=0.0 ACQCHAN=CH1\n"
        "  >EMEAS ID=255.025 CHTYPE=EY X=22.4 Y=-44.7 X2=-22.4 "
        "Y2=44.7 ACQCHAN=CH2\n\n"
    )
    spectra = '>=SPECTRASECT\n  SECTID="demo"\n  NCHAN=5\n\n'
    tail = ">END\n"

    p.write_text(
        head + info + definemeas + spectra + tail,
        encoding="utf-8",
    )


def test_define_meas_from_file_parses_sample(tmp_path: Path) -> None:
    p = tmp_path / "sample.edi"
    _write_sample_edi(p)

    dm = DefineMeas.from_file(p)
    assert dm is not None
    assert dm.maxchan == 7
    assert dm.units == "M"
    assert dm.reftype == "CART"

    # -22 - 22/60 - 14.9/3600 = -22.370805555...
    assert dm.reflat is not None
    assert isclose(dm.reflat, -22.3708055, rel_tol=0.0, abs_tol=1e-6)

    # 139 + 11/60 + 19.1/3600 = 139.18863888...
    assert dm.reflong is not None
    assert isclose(dm.reflong, 139.1886388, rel_tol=0.0, abs_tol=1e-6)

    assert dm.refelev == 200

    # counts
    assert len(dm.hmeas) == 3
    assert len(dm.emeas) == 2

    # spot-check H
    h0 = dm.hmeas[0]
    assert isinstance(h0, Hmeasurement)
    assert h0.id == "251.025"
    assert h0.chtype == "HX"
    assert h0.x == 8.5
    assert h0.y == 8.5
    assert h0.azm == 0.0
    assert h0.acqchan == "CH3"

    # spot-check E
    e1 = dm.emeas[1]
    assert isinstance(e1, Emeasurement)
    assert e1.chtype == "EY"
    assert e1.x == 22.4
    assert e1.y == -44.7
    assert e1.x2 == -22.4
    assert e1.y2 == 44.7
    assert e1.acqchan == "CH2"


def test_define_meas_write_roundtrip(tmp_path: Path) -> None:
    dm = DefineMeas()
    dm.maxchan = 5
    dm.units = "M"
    dm.reftype = "CART"
    dm.reflat = -22.3708055
    dm.reflong = 139.1886388
    dm.refelev = 200

    dm.hmeas.append(
        Hmeasurement(
            id="251.025",
            chtype="hx",
            x=8.5,
            y=8.5,
            azm=0,
            acqchan="CH3",
        )
    )
    dm.emeas.append(
        Emeasurement(
            id="254.025",
            chtype="EX",
            x=-50.0,
            y=0.0,
            x2=50.0,
            y2=0.0,
            acqchan="CH1",
        )
    )

    # compose a minimal full-EDI for from_file
    text = (
        ">HEAD\n\n"
        + "".join(dm.write())
        + ">=SPECTRASECT\n"
        + "  NCHAN=5\n\n"
        + ">END\n"
    )
    p = tmp_path / "roundtrip.edi"
    p.write_text(text, encoding="utf-8")

    dm2 = DefineMeas.from_file(p)
    assert dm2.maxchan == 5
    assert dm2.units == "M"
    assert dm2.reftype == "CART"
    assert isclose(dm2.reflat, dm.reflat, abs_tol=1e-6)
    assert isclose(dm2.reflong, dm.reflong, abs_tol=1e-6)
    assert dm2.refelev == 200
    assert len(dm2.hmeas) == 1
    assert len(dm2.emeas) == 1
    assert dm2.hmeas[0].chtype == "HX"
    assert dm2.emeas[0].chtype == "EX"


def test_mixins_return_expected(tmp_path: Path) -> None:
    p = tmp_path / "sample_mixins.edi"
    _write_sample_edi(p)

    class HostAll(MeasMixin):
        pass

    class HostD(DefineMeasMixin):
        pass

    class HostE(EMeasMixin):
        pass

    class HostH(HMeasMixin):
        pass

    dm = HostAll.from_file(p)
    assert isinstance(dm, DefineMeas)

    dm2 = HostD.from_file(p)
    assert isinstance(dm2, DefineMeas)

    e_list = HostE.from_file(p)
    h_list = HostH.from_file(p)

    assert all(isinstance(e, Emeasurement) for e in e_list)
    assert all(isinstance(h, Hmeasurement) for h in h_list)
    assert len(e_list) == 2
    assert len(h_list) == 3
